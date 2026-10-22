"""Audit all raw windows and fit each process repeat independently."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import statistics as st
import struct
from analyze import fit,write_csv,audit_windows,audit_rejections,T95_DF4

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);args=ap.parse_args()
    root=args.run.resolve();assert (root/'ALL_GROUPS_COMPLETE').exists()
    summaries=[];fits=[];points=[];inventory=[]
    for group in sorted(root.glob('k*_*')):
        if not group.is_dir():continue
        assert (group/'COMPLETE').exists()
        meta=json.loads((group/'dataset.json').read_text())
        raw=(group/'inputs.bin').read_bytes();distances=(group/'exact_distances.i16').read_bytes()
        assert len(raw)==meta['oracle_pairs']*64
        assert len(distances)==meta['oracle_pairs']*2
        unique={raw[i:i+32] for i in range(0,len(raw),32)}
        assert len(unique)>=100000,(group,len(unique))
        exact=[d[0] for d in struct.iter_unpack('<h',distances)]
        assert sum(exact)==meta['distance_checksum']
        inventory.append(dict(group=group.name,pairs=meta['oracle_pairs'],unique_DNA_sequences=len(unique),
            input_sha256=hashlib.sha256(raw).hexdigest(),distance_sha256=hashlib.sha256(distances).hexdigest(),
            independent_DP_pairs=meta['independent_DP_pairs']))
        rows=list(csv.DictReader((group/'trials.csv').open()));audit_rejections(group,rows)
        for row in rows:
            stem=f"{row['algorithm']}_r{row['repeat']}_m{row['multiplier']}"
            audit_windows(group/(stem+'.windows.bin'),row)
            count=int(row['pairs_per_window']);expected=0;passing=0
            for w in range(meta['windows']):
                ds=exact[w*meta['max_pairs']:w*meta['max_pairs']+count]
                passing+=sum(d<=meta['K'] for d in ds)
                expected+=sum(d+1 for d in ds if d<=meta['K'])
            assert passing==int(row['passing_pairs'])
            assert expected==int(row['expected_checksum'])
        for algorithm in sorted({r['algorithm'] for r in rows}):
            selected=[r for r in rows if r['algorithm']==algorithm];repeat_fits=[]
            for repeat in range(meta['repeats']):
                pr=sorted([r for r in selected if int(r['repeat'])==repeat],key=lambda r:int(r['N']))
                assert len(pr)==5
                x=[int(r['N']) for r in pr]
                f=fit(x,[int(r['cpu_ns']) for r in pr])
                wf=fit(x,[int(r['wall_ns']) for r in pr]);cf=fit(x,[int(r['user_cycles']) for r in pr])
                assert f['slope']>0
                repeat_fits.append(f)
                fits.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,repeat=repeat,
                    slope_ns_per_pair=f['slope'],intercept_ns=f['intercept'],r2=f['r2'],
                    within_fit_stderr=f['stderr'],max_abs_residual_ns=f['max_abs_residual'],
                    wall_slope_ns_per_pair=wf['slope'],cycle_slope_per_pair=cf['slope']))
                for r in pr:
                    predicted=f['intercept']+f['slope']*int(r['N'])
                    points.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,
                        repeat=repeat,N=int(r['N']),cpu_ns=int(r['cpu_ns']),fitted_cpu_ns=predicted,
                        residual_ns=int(r['cpu_ns'])-predicted))
            slopes=[f['slope'] for f in repeat_fits];mean=st.mean(slopes)
            summaries.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,
                ns_per_pair=mean,ci95_half_ns=T95_DF4*st.stdev(slopes)/math.sqrt(len(slopes)),
                pair_per_second=1e9/mean,r2_min=min(f['r2'] for f in repeat_fits),
                cycle_per_pair=st.mean(f['cycle_slope_per_pair'] for f in fits
                    if f['workload']==meta['workload'] and f['K']==meta['K'] and f['algorithm']==algorithm),
                min_N=min(int(r['N']) for r in selected),max_N=max(int(r['N']) for r in selected),
                passing_fraction_min=min(int(r['passing_pairs'])/int(r['N']) for r in selected),
                passing_fraction_max=max(int(r['passing_pairs'])/int(r['N']) for r in selected),
                accepted_windows=sum(int(r['windows']) for r in selected),
                rejected_windows=sum(int(r['rejected_windows']) for r in selected),
                accepted_l2_misses=sum(int(r['accepted_l2_misses']) for r in selected),accepted_dram_fills=0))
    lookup={(r['workload'],r['K'],r['algorithm']):r for r in summaries}
    for r in summaries:
        r['speedup_vs_original_hot']=lookup[(r['workload'],r['K'],'edlib_original_hot')]['ns_per_pair']/r['ns_per_pair']
        r['speedup_vs_original_stack']=lookup[(r['workload'],r['K'],'edlib_stack_reference')]['ns_per_pair']/r['ns_per_pair']
    write_csv(root/'summary.csv',summaries);write_csv(root/'fits_per_repeat.csv',fits)
    write_csv(root/'fit_points.csv',points);write_csv(root/'input_inventory.csv',inventory)
    payload=dict(summary=summaries,inventory=inventory,
        oracle_pairs=sum(r['pairs'] for r in inventory),
        independent_DP_pairs=sum(r['independent_DP_pairs'] for r in inventory),
        accepted_windows=sum(r['accepted_windows'] for r in summaries),
        rejected_windows=sum(r['rejected_windows'] for r in summaries))
    (root/'analysis.json').write_text(json.dumps(payload,indent=2)+'\n')
    (root/'ANALYSIS_AUDIT_PASS').write_text('Every binary window, rejected attempt, exact checksum, PMU ratio and input size audited.\n')
    for workload in ['near','mixed']:
        for k in [0,2,5,10,20,30,40,50,60,70,80,90]:
            print(workload,k,' '.join(f"{a}={lookup[(workload,k,a)]['ns_per_pair']:.3f}" for a in
                ['supered_single','supered_batch4','edlib_original_hot','edlib_stack_reference']))

if __name__=='__main__':main()
