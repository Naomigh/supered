"""Check immutable inputs/code, then summarize the measured public interfaces."""
from pathlib import Path
import argparse
import hashlib
import json
import statistics as st

ROOT=Path(__file__).resolve().parents[1]
KS=[0,2,5,10,20,30,40,50,60,70,80,90]
METHODS=['supered_single','supered_batch4','edlib_original_hot','edlib_stack_reference']

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);args=ap.parse_args()
    run=args.run.resolve();assert (run/'ANALYSIS_AUDIT_PASS').exists()
    data=json.loads((run/'analysis.json').read_text());rows=data['summary']
    assert len(rows)==96
    lookup={(r['workload'],r['K'],r['algorithm']):r for r in rows}
    manifest=json.loads((run/'source_sha256.json').read_text())
    for name,h in manifest.items():
        p=ROOT/'build_clang'/name.removeprefix('binary/') if name.startswith('binary/') else ROOT/name
        assert hashlib.sha256(p.read_bytes()).hexdigest()==h,name
    previous=json.loads((ROOT/'results/integration_20261022/source_sha256.json').read_text())
    production=[name for name in previous if name.startswith(('include/','src/','generated/'))]
    for name in production:
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==previous[name],name
    for k in KS:
        for w in ['near','mixed']:
            status=(run/f'k{k}_{w}.process_status').read_text()
            assert 'Threads:\t1\n' in status
            assert 'Cpus_allowed_list:\t24\n' in status
            for method in METHODS:
                r=lookup[(w,k,method)]
                assert r['min_N']==131072 and r['max_N']==655360
                assert r['accepted_dram_fills']==0
    ratio={}
    for w in ['near','mixed']:
        ratio[w]={str(k):lookup[(w,k,'supered_single')]['ns_per_pair']/lookup[(w,k,'supered_batch4')]['ns_per_pair'] for k in KS}
    audit=dict(timed_groups=24,methods=4,independent_fits=480,windows_per_point=4096,
        pair_counts=[131072,262144,393216,524288,655360],
        oracle_pairs=data['oracle_pairs'],independent_DP_pairs=data['independent_DP_pairs'],
        accepted_windows=data['accepted_windows'],rejected_windows=data['rejected_windows'],
        accepted_dram_fills=0,min_fit_r2=min(r['r2_min'] for r in rows),
        verified_manifest_files=len(manifest),unchanged_production_files=len(production),
        all_processes_single_thread_cpu24=True,single_divided_by_batch=ratio)
    (run/'FINAL_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    text=[]
    def add(line=''):text.append(line)
    add('supered 单 pair / 四 pair / 原始 Edlib：128 bp cache 内测速')
    add('单位 ns/pair；数值为五次独立斜率的均值，± 是 t(4) 95% CI 半宽。')
    add()
    add('每个点 N=131072,262144,393216,524288,655360；每窗32..160个独立pair，4096窗。')
    add('原始 Edlib 是预建 Peq/target codes 后的原始 NW 计算核心；stack 仅更换 block 存储。')
    add('supered 直接调用公开 distance()/distance4()，编译器正常优化，没有额外的 pair 级 noinline。')
    add('本轮四 pair 使用四个不同 query；先前约62ns的数据包含共享 query，不能把跨实验差异全归因于接口。')
    add('K=0/2 使用 LEAP；K=5..60 使用64-bit移动带；K=70..90 使用完整双block。')
    for w in ['near','mixed']:
        add();add('近似输入，真实距离<=K' if w=='near' else '混合输入，一半近似、一半独立随机')
        add('K        supered_single       supered_batch4        original_hot        stack_reference')
        for k in KS:
            cells=[f"{lookup[(w,k,m)]['ns_per_pair']:9.3f} ± {lookup[(w,k,m)]['ci95_half_ns']:6.3f}" for m in METHODS]
            add(f'{k:<3}  '+'    '.join(cells))
    add();add('测量边界与限制')
    add('编码、Peq、上下文分配、候选批拷入cache和预热在计时外；输入读取、转置、Eq、LEAP masks/reset、DP、阈值与checksum计时。')
    add('单pair结果是循环调用公开接口的摊销成本；没有用时钟包住孤立一次调用。')
    add('near 的K=2数据只有替换；更高K的近似数据交替包含替换和平衡indel。随机目标在较大K时也可能通过；histogram与passing_fraction保存实际分布。')
    add('所有方法逐项使用同一输入。每组至少10^5条不同DNA；输入32-byte bit-plane与精确距离原始文件均保存。')
    add('每个pair对照原始Edlib完整NW距离；每窗另取一个pair对照独立标量DP；每窗checksum验证。')
    add('cache gate只拒绝非零DRAM/IO demand fills；同一输入重测。保留L2 miss，不按耗时筛选。')
    add('本机VMware guest，单计算线程绑定CPU24；未声称独占物理核心或所有访问命中L1。')
    add(f"验证pair={data['oracle_pairs']:,}；独立DP={data['independent_DP_pairs']:,}；接受窗={data['accepted_windows']:,}；拒绝尝试={data['rejected_windows']:,}；最小R²={audit['min_fit_r2']:.8f}。")
    add('热/冷cache正负对照见cache_control.log；完整原始数据、拟合点/残差和480个独立拟合均保留。')
    add('两次准备阶段试跑保留在single_pair_large_window_pilot_20261022与single_pair_reencode_pilot_20261022，没有合并入本次正式结果。')
    add();add('主要文件：summary.csv / fits_per_repeat.csv / fit_points.csv / input_inventory.csv / FINAL_AUDIT.json')
    add('各k*_near,mixed目录：inputs.bin / exact_distances.i16 / trials.csv / *.windows.bin / rejected_windows.csv / dataset.json / distance_histogram.csv')
    add('source_snapshot、benchmark_single_pair.artifact、libsupered.a.artifact、compile_commands.json、source_sha256.json保存当次构建。')
    add('disassembly.txt与compute_k5/k70_m0/m1.asm保存单pair广播及四pair执行路径。')
    (run/'REPORT.txt').write_text('\n'.join(text)+'\n')
    (run/'visualization_data.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
