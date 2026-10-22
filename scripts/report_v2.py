"""Produce the final reproducibility audit and Chinese report after all groups pass."""
from pathlib import Path
import csv,hashlib,json,statistics,subprocess
root=Path(__file__).resolve().parents[1]
run=root/'results/run_v2_20261022'
assert (run/'ANALYSIS_AUDIT_PASS').exists()
data=json.loads((run/'analysis.json').read_text())
rows=data['summary']; lookup={(r['workload'],r['K'],r['algorithm']):r for r in rows}
old=root.parent/'results/run_20261022'
oldrows=list(csv.DictReader((old/'summary.csv').open()))
oldlookup={(r['workload'],int(r['K']),r['algorithm']):r for r in oldrows}
inventory=data['inventory']; previous={r['group']:r for r in csv.DictReader((old/'input_inventory.csv').open())}
for r in inventory:assert r['input_sha256']==previous[r['group']]['input_sha256']
manifest=json.loads((run/'source_sha256.json').read_text())
bad=[]
for name,expected in manifest.items():
    path=run/(name.removeprefix('binary/')+'.executable') if name.startswith('binary/') else root/name
    # Analysis/report scripts can evolve after launch; timed C++ and executables must not.
    if name.startswith(('include/','src/','tests/','bench/','generated/','binary/')):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:bad.append(name)
assert not bad,bad
for origin in ['edlib_origin','leap_origin']:
    assert not subprocess.check_output(['git','-C',str(root/origin),'status','--porcelain'],text=True)
new=[r for r in rows if r['algorithm']=='myers_2x64x4']
profile=json.loads((root/'results/profile_final/summary.json').read_text())
names=[('edlib_original_hot','Edlib hot'),('edlib_original_stack','Edlib stack'),
       ('myers_2x64x4_v1','旧双 block'),('myers_2x64x4','新双 block'),
       ('myers_128x4','完整 128-bit'),('leap_corrected_avx2','修正 LEAP')]
def fmt(r):return f"{r['ns_per_pair']:.3f} ± {r['ci95_half_ns']:.3f}"
text=[]
def add(s=''):text.append(s)
add('# 双 block 四 pair Myers 第二轮优化')
add()
add(f"新内核八组正式拟合为 **{min(r['ns_per_pair'] for r in new):.3f}–{max(r['ns_per_pair'] for r in new):.3f} ns/pair**，"
    f"对应单核 **{min(r['pair_per_second'] for r in new)/1e6:.3f}–{max(r['pair_per_second'] for r in new)/1e6:.3f} million pair/s**。"
    f"最大均值加 95% CI 半宽为 {max(r['ns_per_pair']+r['ci95_half_ns'] for r in new):.3f} ns/pair。")
add('这是本机 VMware guest CPU 24 的测量结论；生产入口仍为 `myers_2x64x4<K>(Quad,Quad)`，统一每条 32-byte 输入。没有增加输入防御、bad-case fallback 或 CPU dispatch。')
add()
add('## 版本与输入冻结')
add()
add('优化前实际代码与二进制冻结在 `../../../versions/v1_20261022/`，`MANIFEST.json` 覆盖 34 个文件；第一轮原始结果保持在 `../../../results/run_20261022/`，冻结目录通过 `original_results` 引用它。上游两个 git 仓库通过符号链接共享，源码保持 clean。')
add('本轮在独立 `work_v2/` 中修改。正式 run 保存 self-contained `source_snapshot/`、全部被测可执行文件、编译命令、源代码 SHA256、执行命令、环境、输入与逐窗口结果。')
add(f"八个输入 SHA256 与第一轮完全一致；重新使用未修改 Edlib API 验证 {data['benchmark_oracle_pairs']:,} 个完整距离。每组不同 DNA sequence 为 {min(r['unique_DNA_sequences'] for r in inventory):,}–{max(r['unique_DNA_sequences'] for r in inventory):,} 条。")
add()
add('## 正式结果')
add()
add('下表全部为同一个 Clang 22.1.2、`-O3 -march=native`、C++20、LTO 二进制中的比较；单位 ns/pair，± 为五轮斜率均值的 t(4) 95% CI 半宽。K=50 按要求不测 LEAP。')
add()
add('| 负载 | max_edit | '+' | '.join(n for _,n in names)+' |')
add('|---|---:|'+'---:|'*len(names))
for workload in ['near','cross']:
    for k in [0,2,5,50]:
        cells=[fmt(lookup[(workload,k,a)]) if (workload,k,a) in lookup else '—' for a,_ in names]
        add(f"| {workload} | {k} | "+' | '.join(cells)+' |')
add()
add('`near` 的 pair 均在阈值内；`cross` 为四条 query 与候选全交叉，包含正常超阈值结果。K=0 near 的窗口内序列相同，但跨窗口 query 不同。K=2 near 只含替换，其 indel 正例由独立测试覆盖，不将该 LEAP 速度外推到所有 K=2 indel 分布。')
add()
add('| 负载 | K | 冻结 GCC 旧双 block | 同轮 Clang 旧双 block | Clang 新双 block | 同编译器加速比 |')
add('|---|---:|---:|---:|---:|---:|')
for workload in ['near','cross']:
    for k in [0,2,5,50]:
        v1=lookup[(workload,k,'myers_2x64x4_v1')]['ns_per_pair'];v2=lookup[(workload,k,'myers_2x64x4')]['ns_per_pair']
        gcc=float(oldlookup[(workload,k,'myers_2x64x4')]['ns_per_pair'])
        add(f'| {workload} | {k} | {gcc:.3f} | {v1:.3f} | {v2:.3f} | {v1/v2:.3f}× |')
add()
add('因此不能把从冻结 GCC 旧版到当前 Clang 新版的全部差异归为算法/代码收益。原始 Edlib hot 保留两 block 内部分配，stack 参考仅将其 block 存储放到栈上；两者都使用预建 Peq/target codes，没有用 `edlibAlign()` API 开销充当原始核心耗时。LEAP 的 pair-dependent masks/reset/run 仍计时。')
add()
add('## 先 profiling 后修改的依据')
add()
add('初始 GCC 热循环约 138 ns/pair，反汇编在每列重载 query 位平面，并反复将 Ph/Mh 边界编码成 signed hin 再解码。L2 miss 很少，没有把瓶颈假设为 DRAM 带宽。初始 stat 的六事件组发生复用，只作探索；最终四事件固定组不复用，所有正式吞吐来自独立 cache 门控拟合。')
add('`../profile_baseline/` 保留修改前 perf.data、annotate、完整二进制和汇编；`../profile_final/` 保存顺序运行的三次 stat 与采样。下表是长热点循环每轮 7,680,000 pairs 的整进程 PMU 归一化值，含很小的启动成本，用于定位而非代替正式拟合。')
add()
add('| 内核/对照 | cycles/pair | instructions/pair | IPC | 执行未完成等待/cycles | load 未完成等待/cycles |')
add('|---|---:|---:|---:|---:|---:|')
for r in profile:
    add(f"| {r['algorithm']} | {r['cycles_per_pair']:.2f} | {r['instructions_per_pair']:.2f} | {r['ipc']:.3f} | {100*r['execution_wait_fraction']:.2f}% | {100*r['load_wait_fraction']:.3f}% |")
add()
add('## 实际采用的修改')
add()
add('1. `equality4` 用 XOR + ternary truth table 0x09 直接生成匹配位，减少旧 Eq 表达式的补码步骤。')
add('2. 跨 block 保留原始 Ph/Mh，而不是 `(Ph>>63)-(Mh>>63)` 的 signed hin；高 block 的 Eq 只读取 Mh 的最高位，移位时用 `VPSHLDQ` 直接注入低 block 的最高位。')
add('3. 先剥离低 block 第 0 列，循环同时推进高 block 的前一列与低 block 的当前列，最后排空高 block 第 127 列。二者的递推在每个循环体中独立，边界来自先前保存的 Ph/Mh。')
add('4. 删除每列 score 维护。在最终列，PV/MV 分别表示纵向差分 +1/-1，列顶端 D(0,128)=128，故最终距离为 `128 + popcount(PV_lo)+popcount(PV_hi) - popcount(MV_lo)-popcount(MV_hi)`。最后才应用 K 的精确/超阈值语义。')
add('5. 四条 32-byte 输入用四次 YMM load 后 unpack/permute 转置。拼装仍计时，没有把 pair 相关工作移到预处理。循环两列展开，在 profiling 后选择 Clang 的寄存器分配/调度。')
add('所有配置选择均为编译期常量；内核仍完整计算 128 列，没有靠输入检查、提前相等判断、band 或 early rejection 得到这些成绩。试验模板及未采纳的反向遍历、ZMM block 交错、score 位累积、展开方案保存在 `bench/exploration.hpp` 和原始 `profile_*.time` 中。')
add()
add('## 当前完整 128-bit 状态为什么更慢')
add()
add('ZMM 中是四个 pair 的 `[lo,hi]`，并不是八个独立 pair。AVX512 的整数加法和移位在这里仍以 64-bit qword 为单位；128-bit carry/shift 需要显式连接，因此换成一个 ZMM 没有消除递推依赖。')
add('在当前 Clang 汇编中，完整加法的关键链为 `VPADDQ → unsigned VPCMP → KSHIFTLB → VPMOVM2Q → VPSUBQ → Xh/Ph/Mh`。GCC 则将后两步表达为 masked add。随后 Ph/Mh 各需邻接低 word 的边界注入才能形成下一列 PV/MV。当前 Clang 已把原始 permutation 优化为 `VPSLLDQ 8` 与 `VPSHLDQ 1`，所以单纯将源码改成 lane-local shuffle 没有收益。')
add('完整状态退休指令比旧双 block 更少，却耗费约 1.83 倍 cycles，IPC 从约 3.21 降为 1.14；执行未完成等待约 44.7%，load 等待约 0.03%。结合汇编与对照实验，这支持长算术/掩码依赖与执行调度瓶颈，不支持 DRAM 容量/带宽是主要原因。PMU 等待事件本身不能单独证明某一条指令的精确延迟。')
add('独立对照把 carry 留在向量寄存器，用 `(a&b)|((a|b)&~sum)` 的最高位计算 carry，cycles 从约 852 降为 788/pair，改善约 7.5%；仍远慢于旧双 block 的 466 和新双 block 的 359。只改 shift 的对照几乎不变，两者同时修改也约 788。三个诊断版本各用 16,384 pair 对照原始 Edlib，均正确。')
add('新双 block 通过列错位将两个 64-bit recurrence 暴露为可并行调度的工作，并移除了 score 与 signed 边界转换。这里的结论针对已保存的完整状态实现，不声称所有可能的 128-bit/ZMM 方案都必然更慢，也没有证据把差异归因为 AVX512 降频。LLVM-MCA 的 Znver4 模型明显低估完整状态实测周期，只作定性参考，不用其预测当作测量值。')
add()
add('## 计时、cache 与拟合审计')
add()
add('每窗固定一条或四条 query，流入 4/8/12/16/20 个候选，与全部固定 query 比较。输入加载、解码、Peq/编码预处理和 cache touch 全在计时外；内核 load、转置、状态初始化、Eq、全部递推、最终 score、threshold 与结果 checksum 在计时内。')
add('每条拟合固定窗口数，只改变每窗 pair 数；五个 N 分别拟合 `T=a+bN`，斜率 b 为 ns/pair。每组五轮，随机化算法与 N 顺序。固定 timer/RDPMC/window 调用开销进入截距，每 pair 的必要成本保留在斜率中。')
add(f"所有拟合点至少 {min(r['min_N'] for r in rows):,} pairs；最大点 {max(r['max_N'] for r in rows):,} pairs。全部拟合最小 R²={min(r['r2_min'] for r in rows):.6f}。")
add(f"本轮接受 {data['accepted_windows']:,} 个窗口，全部 DRAM/IO demand fill 增量为 0；拒绝 {data['rejected_windows']:,} 次尝试并重测相同输入，未跳过 pair。保留全部 L2 miss；允许 L2/L3 命中，不将其称为纯算术或独占 L1。PMU running/enabled 均至少 0.99。")
add('hot-16KiB 对照五轮 demand fills 均为 0；cold-clflush-64MiB 对照约 1.05 million fills/轮。固定 guest CPU 24、单个测量进程顺序运行；VM 无物理核心独占和真实共享拓扑保证。原始 wall time 同时保存，主结论使用 thread CPU 时间。')
add('`analyze.py` 重新求和每个 `<QQQII` 接受窗口记录，核对重试次数、checksum、PMU 调度率、输入大小、不同 DNA 数量，并逐行审计 rejected CSV 的条数/CPU ns/DRAM fills。输入 SHA 和冻结源代码另做交叉审计。')
add()
add('## 正确性与复现')
add()
add('独立 `validate`：131,072 pairs、K=0/2/5/50/128，完整及 bounded 距离对照未修改 Edlib API，含随机/重复/同一序列、替换、平衡 indel、不同 query lanes、阈值边界；2,551 pairs 额外对照标量 DP，40,000 组原生 128-bit 加法/移位原语检查。`validate_indels` 再检查 16,384 个跨 63/64 位的 K=2 indel 正例，每个都与标量 DP 对照。')
add(f"正式输入还独立验证 {data['benchmark_oracle_pairs']:,} 个 pair；所有计时窗口结果 checksum 与原始距离一致。旧双 block 头文件仅改 namespace，与冻结头文件逐字核对一致。")
add()
add('```bash\ncmake -S . -B build_clang -DCMAKE_BUILD_TYPE=Release \\\n  -DCMAKE_C_COMPILER=clang-22 -DCMAKE_CXX_COMPILER=clang++-22 \\\n  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=/usr/bin/ld.lld-22\ncmake --build build_clang -j4\nflock /tmp/edlib_strong_measure.lock python3 scripts/run_experiment.py \\\n  results/new_run --build-dir build_clang --input-root ../results/run_20261022\npython3 scripts/analyze.py results/new_run\n```')
add()
add('正式路径 `run_v2_20261022` 已有结果，复现应使用新的输出目录。`REPORT.md` 为报告，`summary.csv`/`fits_per_repeat.csv`/`fit_points.csv` 为估计与残差，`trials.csv`、`*.windows.bin`、`rejected_windows.csv`、`inputs.bin` 是未经汇总的原始证据。Profiling 的 stat/采样与正式计时目录分开；标为 overlap exploration 的早期试验不用于结论。')
(run/'REPORT.md').write_text('\n'.join(text)+'\n')
audit=dict(groups=len(inventory),timing_points=len(data['points']),algorithm_groups=len(rows),accepted_windows=data['accepted_windows'],rejected_windows=data['rejected_windows'],input_sha_all_match_v1=True,timed_source_mismatches=bad,all_new_mean_le100=all(r['ns_per_pair']<=100 for r in new),all_new_ci95_upper_le100=all(r['ns_per_pair']+r['ci95_half_ns']<=100 for r in new),min_r2=min(r['r2_min'] for r in rows))
(run/'FINAL_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
