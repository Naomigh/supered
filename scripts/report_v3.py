"""Audit the timed code and identical corpora, then write the third-round report."""
from pathlib import Path
import csv,hashlib,json,subprocess
root=Path(__file__).resolve().parents[1];run=root/'results/run_v3_20261022'
assert (run/'ANALYSIS_AUDIT_PASS').exists()
data=json.loads((run/'analysis.json').read_text());rows=data['summary']
lookup={(r['workload'],r['K'],r['algorithm']):r for r in rows}
inventory=data['inventory'];v1=root.parent/'results/run_20261022';v2=root.parent/'work_v2/results/run_v2_20261022'
for earlier in [v1,v2]:
    previous={r['group']:r for r in csv.DictReader((earlier/'input_inventory.csv').open())}
    for r in inventory:assert r['input_sha256']==previous[r['group']]['input_sha256']
manifest=json.loads((run/'source_sha256.json').read_text());bad=[]
for name,h in manifest.items():
    if name=='CMakeLists.txt' or name.startswith(('include/','src/','tests/','bench/','generated/','binary/')):
        p=run/(name.removeprefix('binary/')+'.executable') if name.startswith('binary/') else root/name
        if hashlib.sha256(p.read_bytes()).hexdigest()!=h:bad.append(name)
assert not bad,bad
for repo in ['edlib_origin','leap_origin']:
    assert not subprocess.check_output(['git','-C',str(root/repo),'status','--porcelain'],text=True)
new=[r for r in rows if r['algorithm']=='myers_banded64x4']
ratios=[lookup[(r['workload'],r['K'],'myers_2x64x4_v2')]['ns_per_pair']/r['ns_per_pair'] for r in new]
profiling=json.loads((root/'results/profile_final/summary.json').read_text())
text=[]
def add(s=''):text.append(s)
add('# 第三轮：相同阈值语义下的 64-bit 移动带优化')
add()
add(f"八组正式 cache 内拟合中，`myers_banded64x4<K>` 为 **{min(r['ns_per_pair'] for r in new):.3f}–{max(r['ns_per_pair'] for r in new):.3f} ns/pair**，"
    f"即单核 **{min(r['pair_per_second'] for r in new)/1e6:.3f}–{max(r['pair_per_second'] for r in new)/1e6:.3f} million pair/s**。"
    f"相对同轮、同二进制中的 v2，速度为 **{min(ratios):.3f}–{max(ratios):.3f} 倍**。")
add('这里加速的是给定 K 的有界全局距离接口。它利用等长序列的阈值约束缩小 DP 带宽，仍返回精确距离或 -1。不设上界的完整距离并未得到这个加速数字；K>=64 使用保留的完整双 block。')
add()
add('## 冻结与代码位置')
add()
add('此前约 98 ns/pair 的 v2 已先冻结在 `../../../versions/v2_20261022/`，代码和二进制共 42 个文件有 SHA256 并设为只读；`original_results` 引用完整 v2 原始实验。开发在独立 `work_v3/` 中完成。两个原始 git 仓库仍保持 clean，工作区通过符号链接引用它们。')
add('新入口为 `include/strong/banded.hpp` 中的 `myers_banded64x4<K>`，默认四列展开；旧的 `myers_2x64x4<K>` 和完整 128-bit 状态源码逐字保持 v2 版本。v2 参考头文件只改变 namespace 与限定 wrapper 调用，避免 ADL 歧义。')
add('`DESIGN_V3.md` 给出带宽证明、逐步递推与最终 score 公式。生产代码没有加入长度/字符/空指针检查、相等特判、CPU dispatch 或 bad-case fallback。K 的编译期选择属于算法所需的带宽选择。')
add()
add('## 正式性能')
add()
add('单位为 ns/pair，± 为五轮独立斜率的 t(4) 95% CI 半宽。所有方法在同一个 Clang 22.1.2 / C++20 / `-O3 -march=native` / LTO 二进制中测量；固定 guest CPU 24。K=50 不纳入 LEAP。')
add()
algs=[('edlib_original_hot','原始 hot'),('edlib_original_stack','原始 stack'),('myers_2x64x4_v1','v1'),('myers_2x64x4_v2','v2'),('myers_banded64x4','新移动带'),('myers_128x4','完整 128-bit'),('leap_corrected_avx2','修正 LEAP')]
add('| 负载 | K | '+' | '.join(n for _,n in algs)+' |')
add('|---|---:|'+'---:|'*len(algs))
for w in ['near','cross']:
    for k in [0,2,5,50]:
        cells=[]
        for a,_ in algs:
            r=lookup.get((w,k,a));cells.append(f"{r['ns_per_pair']:.3f} ± {r['ci95_half_ns']:.3f}" if r else '—')
        add(f'| {w} | {k} | '+' | '.join(cells)+' |')
add()
add('原始 hot 是预建 Peq/target codes 的原始 NW 核心，保留内部两 block 分配；stack 参考仅更换 block 存储方式。没有将 `edlibAlign()` 的编码/API 开销混入原始核心测速。LEAP 源码与路由未作修改，仍测量其 pair-dependent mask/reset/run。')
add('`near` 真实距离均不超过 K；`cross` 是四条固定 query 与候选的全交叉，包含正常超阈值结果。K=0 near 在单窗内相同，但跨窗固定 query 不同。K=2 near 的性能数据只有替换，长 indel 的正确性由另外的独立测试覆盖，不把该 LEAP 速度外推到所有 K=2 分布。')
add()
add('## 为什么能进一步减少计算')
add()
add('等长 NW 路径的 insertion 与 deletion 数相等。路径若偏离主对角线 d 格，为到达该偏移再返回终点，至少需 2|d| 次 indel；因此所有距离 <=K 的路径都有 `|i-j|<=floor(K/2)`。K<=63 时最大偏移 31，均可容纳在 64-bit 移动带中。K=64 已可能需要偏移 32，因此使用完整双 block。')
add('v2 每列计算两个 64-bit block，靠跨列交错减少边界依赖。新版本每列只计算一个 64-bit 带，不再传递两个 block 的 Ph/Mh；新增 query 位窗口的每列移动，仍只使用原来的 32-byte 输入，全部窗口构造与移动均计时。所有 pair 都执行 128 列，没有数据相关的提前退出。')
add('将普通 Myers 的 Ph/Mh 左移与带窗口右移合并后，令 `V=Xv>>1`，可写成 `PV_next=Mh|~(V|Ph)`、`MV_next=V&Ph`。最高 PV 位为 1、最高 MV 位为 0 的不变量保证新最高 PV 位天然为 1，省去初版每列的一次 OR。')
add('移动 origin 的 score 不能照搬完整列公式。初始 origin=-32，S=32；每列 `S_next=S+1-(Xv&1)`。用 VPSHRDQ 收集这些 bit，每 64 列 popcount 一次；最后 origin=96，目标 row128 距它 32 行：`D=160-sum(hits)+popcount(PV&0xffffffff)-popcount(MV&0xffffffff)`。')
add('负行的初始代价与合法 insertion 的投影保证不会产生低于真实 NW 的值；窗口边界只限制路径，所有距离 <=K 的路径完整保留。因此阈值内结果精确，超阈值返回 -1，不是近似距离或启发式过滤。详细证明见 `DESIGN_V3.md`。')
add()
add('## Profiling 与针对性对照')
add()
add('根据 v2 的既有 profiling，并用本轮重测确认，load 等待很少，主要成本是每列两个 block 的递推与执行依赖。先保留该基线，再测试移动带、消去最高位 OR、展开程度。下表是同一热点循环每次 12,800,000 pairs、三次 perf stat 的整进程 PMU 归一化值；包含少量启动成本，正式吞吐仍采用上面的斜率拟合。')
add()
add('| 版本 | cycles/pair | instructions/pair | IPC | execution wait/cycles | load wait/cycles |')
add('|---|---:|---:|---:|---:|---:|')
for r in profiling:
    add(f"| {r['algorithm']} | {r['cycles_per_pair']:.2f} | {r['instructions_per_pair']:.2f} | {r['ipc']:.3f} | {100*r['execution_wait_fraction']:.2f}% | {100*r['load_wait_fraction']:.3f}% |")
add()
add('初版移动带约 262 cycles/pair，去掉每列 OR 后约 231，四列展开约 228；v2 约 359。指令数约从 1071 降到 691，表明本次收益来自工作量和依赖链减少，没有以不同 cache 状态解释结果。PMU 等待事件不能单独给出某条指令的精确延迟，结论同时依靠源码等价推导和独立对照。')
add('小样本 CPU 时间拟合曾出现计时器固定成本的波动，全部原始结果保留为 pilot，没有当作正式成绩或用最低斜率汇报。最终按原有大样本、五轮、五个 N 的规程重新测量；cycle 与 wall time 也全部保留。')
add()
add('## 计时与缓存审计')
add()
add('计时程序的 trial 函数、PerfGate/RDPMC 实现、PackedDNA128 布局和 Sequence 预处理函数与冻结 v2 逐字核对一致。仅增加候选算法及基线。八组输入分别与 v1、v2 的完整 SHA256 核对一致。')
add('每个窗口固定一条或四条 query，流入 4/8/12/16/20 个候选，与全部固定 query 比较。加载新批、编码、Peq 和 cache touch 均在计时外；load/转置、window 形成、Eq、fresh 状态、全部 128 列、score/threshold 和 checksum 在计时内。固定窗口数改变 N，拟合 `T=a+bN`，用 b 表示 ns/pair。')
add(f"每点至少 {min(r['min_N'] for r in rows):,} pairs，最大 {max(r['max_N'] for r in rows):,}；每组包含 {min(r['unique_DNA_sequences'] for r in inventory):,}–{max(r['unique_DNA_sequences'] for r in inventory):,} 条不同 DNA。全部拟合最小 R²={min(r['r2_min'] for r in rows):.6f}。")
add(f"本轮接受 {data['accepted_windows']:,} 个窗口，观测到的 DRAM/IO demand fills 全部为 0；拒绝并重测 {data['rejected_windows']:,} 次尝试，不跳过输入。保留 L2 miss，允许 L2/L3 命中。每 trial PMU running/enabled>=0.99。正/负缓存对照、原始时间与 PMU 数据均保存。")
add('审计重新求和所有窗口二进制记录，逐行核对 rejected CSV 的次数/CPU 时间/DRAM fills，核对 checksum、输入大小和 PMU 调度率。这里的 cache 条件不等于独占物理核心、所有访问命中 L1 或没有硬件预取；环境仍是 VMware。')
add()
add('## 正确性范围')
add()
add('主验证：131,072 pairs，K=0/2/5/50/128，比较新移动带、完整双 block、完整 128-bit、原始 Edlib API/核心和适用的 LEAP；2,551 pairs 再对照独立标量 DP。另有 40,000 组原生 128-bit 原语检查。')
add('新增 229,376 pairs：在位置 0、60、120 上分别枚举全部二元 8-base 字符串组合，共 196,608；再生成 32,768 个最大 indel excursion 64 的长位移案例。K=0/2/5/50/62/63/64/128，展开 1/2/4/8 均逐 pair 对照原始 Edlib，33,534 pairs 另与独立 DP 对照。原有 16,384 个跨 63/64 位的 K=2 indel 正例也改为检查新候选并通过。')
add(f"正式输入另外逐项用原始 API 验证 {data['benchmark_oracle_pairs']:,} 个 pair；计时后每窗 checksum 都一致。该结果是以上测试范围内的一致性，不声称穷举了全部 128 bp 输入。")
add()
add('## 复现和原始文件')
add()
add('```bash\ncmake -S . -B build_clang -DCMAKE_BUILD_TYPE=Release \\\n  -DCMAKE_C_COMPILER=clang-22 -DCMAKE_CXX_COMPILER=clang++-22 \\\n  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=/usr/bin/ld.lld-22\ncmake --build build_clang -j4\nflock /tmp/edlib_strong_measure.lock python3 scripts/run_experiment.py \\\n  results/new_run --build-dir build_clang --input-root ../results/run_20261022\npython3 scripts/analyze.py results/new_run\n```')
add()
add('本目录保存 source_snapshot、实际可执行文件、编译参数、输入与环境；每组 `trials.csv`、`*.windows.bin`、`rejected_windows.csv` 和 `inputs.bin` 是原始数据。`summary.csv`、`fits_per_repeat.csv`、`fit_points.csv` 为拟合与残差，`FINAL_AUDIT.json` 为最终审计。`../profile_final/` 有本轮 perf stat/data/annotate 与汇编；`../../../versions/v2_20261022/` 为冻结基线。')
(run/'REPORT.md').write_text('\n'.join(text)+'\n')
audit=dict(groups=len(inventory),timing_points=len(data['points']),algorithm_groups=len(rows),accepted_windows=data['accepted_windows'],rejected_windows=data['rejected_windows'],all_inputs_match_v1_and_v2=True,timed_source_mismatches=bad,baseline_audit=json.loads((root/'results/profile_final/baseline_audit.json').read_text()),speedup_vs_v2_min=min(ratios),speedup_vs_v2_max=max(ratios),new_ns_min=min(r['ns_per_pair'] for r in new),new_ns_max=max(r['ns_per_pair'] for r in new),min_r2=min(r['r2_min'] for r in rows))
(run/'FINAL_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
