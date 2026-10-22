from pathlib import Path
import argparse
import json

NAMES = {'edlib_original_hot':'原始 Edlib 核心', 'edlib_original_stack':'原始 Edlib 栈存储',
         'myers_2x64x4':'SIMD 2×64 ×4', 'myers_128x4':'SIMD 128 ×4',
         'leap_corrected_avx2':'修复版 LEAP AVX2'}

TEMPLATE = r'''
import {BarChart, Divider, Grid, H1, H2, LineChart, Row, Select, Stack, Stat, Table, Text, useCanvasState, useHostTheme} from 'cursor/canvas';
const report = __DATA__;
const names: Record<string,string> = __NAMES__;
const order = ['edlib_original_hot','edlib_original_stack','myers_2x64x4','myers_128x4','leap_corrected_avx2'];
export default function EdlibBenchmark() {
  const {tokens:t} = useHostTheme();
  const [mode,setMode] = useCanvasState('edlib-workload','near');
  const [threshold,setThreshold] = useCanvasState('edlib-K','5');
  const [chosen,setChosen] = useCanvasState('edlib-algorithm','myers_2x64x4');
  const rows = report.summary.filter(r=>r.workload===mode && r.K===Number(threshold)).sort((a,b)=>order.indexOf(a.algorithm)-order.indexOf(b.algorithm));
  const fast = rows.find(r=>r.algorithm==='myers_2x64x4')!;
  const algorithm = rows.some(r=>r.algorithm===chosen) ? chosen : 'myers_2x64x4';
  const current = rows.find(r=>r.algorithm===algorithm)!;
  const points = report.points.filter(r=>r.workload===mode && r.K===Number(threshold) && r.algorithm===algorithm);
  const sizes = Array.from(new Set(points.map(r=>r.N))).sort((a,b)=>a-b);
  const mean = (n:number,key:'cpu_ns'|'fitted_cpu_ns')=>{
    const rr=points.filter(r=>r.N===n);return rr.reduce((s,r)=>s+r[key],0)/rr.length/1e6;
  };
  const inventory=report.inventory.find(r=>r.group===`k${threshold}_${mode}`)!;
  return <Stack gap={20}>
    <H1>四 pair Myers：两个 64-bit block 更快</H1>
    <Text>128 bp DNA，统一 32-byte 位平面。单个 guest CPU 24，cache 内计算，预处理与数据加载排除。结果为 T=a+bN 的斜率。</Text>
    <Row gap={12} wrap>
      <Select value={mode} onChange={setMode} options={[{value:'near',label:'阈值内相近 pair'},{value:'cross',label:'四 query 全交叉，含超阈值'}]}/>
      <Select value={threshold} onChange={setThreshold} options={['0','2','5','50'].map(k=>({value:k,label:`max_edit = ${k}`}))}/>
    </Row>
    <Grid columns={3} gap={18}>
      <Stat value={`${fast.ns_per_pair.toFixed(2)} ns`} label='双 block SIMD / pair'/>
      <Stat value={`${(fast.pair_per_second/1e6).toFixed(2)} M`} label='pair / CPU-second'/>
      <Stat value={`${fast.speedup_vs_original_stack.toFixed(2)}×`} label='相对原始栈存储核心'/>
    </Grid>
    <Text weight='medium'>通用 Myers 尚未达到 100 ns/pair；完整 128-bit 状态在本轮实测中更慢。</Text>
    <Table headers={['实现','ns/pair ± 95% CI','M pair/s','相对原始核心','最小 R²']}
      rows={rows.map(r=>[names[r.algorithm],`${r.ns_per_pair.toFixed(2)} ± ${r.ci95_half_ns.toFixed(2)}`,
        (r.pair_per_second/1e6).toFixed(2),`${r.speedup_vs_original_hot.toFixed(2)}×`,r.r2_min.toFixed(6)])}
      columnAlign={['left','right','right','right','right']}/>
    <BarChart categories={rows.map(r=>names[r.algorithm])} series={[{name:'ns/pair，越低越快',data:rows.map(r=>r.ns_per_pair),color:t.accent.primary}]} horizontal height={240} valueSuffix=' ns'/>
    <Divider/>
    <H2>检查线性拟合</H2>
    <Select value={algorithm} onChange={setChosen} options={rows.map(r=>({value:r.algorithm,label:names[r.algorithm]}))}/>
    <LineChart categories={sizes.map(n=>`${(n/1e6).toFixed(3)} M pair`)} series={[
      {name:'五轮实际 CPU 时间均值',data:sizes.map(n=>mean(n,'cpu_ns')),color:t.accent.primary},
      {name:'五条拟合线的均值',data:sizes.map(n=>mean(n,'fitted_cpu_ns')),color:t.text.tertiary}
    ]} height={280} valueSuffix=' ms'/>
    <Text tone='secondary'>b={current.ns_per_pair.toFixed(3)} ns/pair；平均截距={(current.intercept_ns_mean/1e6).toFixed(3)} ms。每轮固定窗口数，改变每窗 pair 数；5 个 N、5 轮，95% CI 使用五条斜率和 t(4)。</Text>
    <Grid columns={2} gap={24}>
      <Stack gap={8}><H2>正确性与输入</H2>
        <Text>独立随机/结构验证 131,072 pair；另有 16,384 个跨 word 的 K=2 插入删除正例，全部与原始 API 和相应独立 DP 一致。</Text>
        <Text>正式数据集共逐 pair 验证 {report.benchmark_oracle_pairs.toLocaleString()} 个记录；当前组 {inventory.unique_DNA_sequences.toLocaleString()} 条不同 DNA，最大点 {inventory.max_point_pairs.toLocaleString()} pair。</Text>
        <Text tone='secondary'>near 的 K=2 性能输入只有替换；K=5/50 在编辑预算达到 4 时加入 indel。cross 的一半候选随机、一半来自固定 query 的邻域，与四条 query 全交叉。</Text>
      </Stack>
      <Stack gap={8}><H2>缓存与计时边界</H2>
        <Text>{report.accepted_windows.toLocaleString()} 个接受窗口的 demand DRAM/IO fill 增量均为 0；{report.rejected_windows.toLocaleString()} 个失败尝试保留并对同一输入重测。允许 L2/L3 命中。</Text>
        <Text>包含 fresh 状态、SIMD 拼装、pair-dependent mask 和结果累计；排除编码、Peq 构建、候选加载和 cache touch。原始核心保留内部两 block 分配；另列栈存储对照。</Text>
        <Text tone='secondary'>合成 DNA，VMware guest；未控制宿主机核心共享或物理频率。K=50 按要求不运行 LEAP。</Text>
      </Stack>
    </Grid>
  </Stack>;
}
'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('--canvas',type=Path)
    args=ap.parse_args();root=args.run.resolve()
    a=json.loads((root/'analysis.json').read_text())
    assert (root/'ANALYSIS_AUDIT_PASS').exists()
    lines=['# 四 pair SIMD Myers 实测结果（2026-10-22）','',
        '两个 64-bit block 的四 pair 内核约为 135 ns/pair；完整 128-bit 状态约为 259 ns/pair。',
        '通用 Myers 约 7.4M pair/CPU-second，尚未达到 10M 目标。所有数值为线性拟合斜率，非单次计时比值。','',
        '以下单位均为 ns/pair，顺序为原始 NW 核心、原始栈存储核心、SIMD 双 block、SIMD 完整 128-bit、修复版 LEAP。','']
    lookup={(r['workload'],r['K'],r['algorithm']):r for r in a['summary']}
    order=list(NAMES)
    for mode,title in [('near','阈值内相近 pair'),('cross','四 query 全交叉，含正常超阈值结果')]:
        lines+=['## '+title,'']
        for k in [0,2,5,50]:
            values=[f"{lookup[(mode,k,name)]['ns_per_pair']:.2f}" if (mode,k,name) in lookup else '未测' for name in order]
            lines.append(f"- K={k}："+' / '.join(values)+'。')
        lines.append('')
    lines+=['## 验证与测量范围','',
        '- 初始独立验证 131,072 pair，包含 40,000 组 128-bit 加法/移位测试和 2,551 个独立 DP 对照。',
        '- 额外 16,384 个跨 63/64 位边界的 K=2 插入删除正例，与两个内核、原始 API、独立 DP 和修复版 LEAP 一致。',
        f"- 正式数据逐 pair 验证 {a['benchmark_oracle_pairs']:,} 个记录；每组不同 DNA 数量均至少 114,688。输入 SHA256 见 input_inventory.csv。",
        f"- 共 {a['accepted_windows']:,} 个接受窗口，观察到的 demand DRAM/IO fill 为 0；{a['rejected_windows']:,} 个失败尝试已记录并重测相同输入。",
        '- 950 个计时点。每组五个 N、五轮独立拟合；所有拟合 R² 至少 0.99718，置信区间及残差完整保留。',
        '- near 的 K=2 数据只有替换；K=5/50 的部分候选包含平衡插入删除。这是合成数据性能结果，不外推为真实 mapping 工作流吞吐。',
        '- 固定 guest CPU 24；VMware 环境，不代表独占宿主机物理核心。主结果使用线程 CPU 时间，wall 时间同样保留。',
        '- 每条拟合线固定窗口数、改变每窗比较数，固定窗口成本进入截距。编码、预建 Peq、流入数据准备、cache touch 全在计时外。',
        '- SIMD 输入装载/拼装、fresh 状态初始化、LEAP pair mask、递推、结果累计在计时内。两个新 Myers 内核均计算完整矩阵后应用阈值。',
        '- 原始核心保留其内部 new/delete；栈存储对照仅修改 block 存储。二者的差异也可能包含编译器优化变化，不能全部归因于 malloc 延迟。',
        '- 缓存事件来自 Linux perf 的 [AMD Zen4 定义](https://raw.githubusercontent.com/torvalds/linux/master/tools/perf/pmu-events/arch/x86/amdzen4/cache.json)；正负对照及所有窗口记录可审计。','',
        '## 文件','',
        '- summary.csv：全部斜率、95% CI、pair/s、R²、相对两个原始基线的加速比。',
        '- fits_per_repeat.csv / fit_points.csv：每轮拟合、每点拟合值及残差。',
        '- k*/trials.csv / *.windows.bin / rejected_windows.csv：逐点及逐窗口原始值。',
        '- k*/inputs.bin / distance_histogram.csv：统一编码输入和 oracle 真实距离。',
        '- source_snapshot / source_sha256.json / benchmark.executable / assembly.txt：测量对应源码、版本、二进制和指令证据。',
        '- 主目录 README.md：完整计时边界、输入格式与复现命令。','']
    (root/'RESULTS.md').write_text('\n'.join(lines))
    if args.canvas:
        content=TEMPLATE.replace('__DATA__',json.dumps(a,ensure_ascii=False,separators=(',',':'))).replace('__NAMES__',json.dumps(NAMES,ensure_ascii=False))
        args.canvas.write_text(content)
        print('Canvas:',args.canvas)
    print('Report:',root/'RESULTS.md')

if __name__=='__main__':main()
