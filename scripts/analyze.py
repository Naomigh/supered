"""Audit raw accepted windows and estimate T=a+b*N, independently per repeat."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import statistics as st
import struct
from collections import defaultdict

T95_DF4 = 2.7764451051977987

def fit(x, y):
    mx, my = st.mean(x), st.mean(y)
    sxx = sum((v-mx)**2 for v in x)
    slope = sum((a-mx)*(b-my) for a,b in zip(x,y))/sxx
    intercept = my - slope*mx
    residual = [b-(intercept+slope*a) for a,b in zip(x,y)]
    sse = sum(v*v for v in residual)
    sst = sum((v-my)**2 for v in y)
    return dict(slope=slope, intercept=intercept,
                r2=1-sse/sst, stderr=math.sqrt(sse/(len(x)-2)/sxx),
                max_abs_residual=max(map(abs,residual)))

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def audit_windows(path, row):
    n=int(row['windows'])
    data=path.read_bytes()
    assert len(data)==n*32, (path,len(data),n)
    values=list(struct.iter_unpack('<QQQII',data))
    for field,index in [('cpu_ns',0),('wall_ns',1),('user_cycles',2),('accepted_l2_misses',3)]:
        assert sum(v[index] for v in values)==int(row[field]),(path,field)
    assert sum(v[4]-1 for v in values)==int(row['rejected_windows'])
    assert int(row['accepted_dram_fills'])==0
    assert row['checksum']==row['expected_checksum']
    assert int(row['running_ns'])/int(row['enabled_ns'])>=.99

def audit_rejections(group,rows):
    totals=defaultdict(lambda:[0,0,0])
    with (group/'rejected_windows.csv').open() as f:
        for r in csv.DictReader(f):
            assert int(r['dram_fills'])>0
            v=totals[(r['algorithm'],r['repeat'],r['multiplier'])]
            v[0]+=1;v[1]+=int(r['cpu_ns']);v[2]+=int(r['dram_fills'])
    for r in rows:
        v=totals.pop((r['algorithm'],r['repeat'],r['multiplier']),[0,0,0])
        assert v==[int(r['rejected_windows']),int(r['rejected_cpu_ns']),int(r['rejected_dram_fills'])]
    assert not totals

def input_inventory(group, meta):
    data=(group/'inputs.bin').read_bytes()
    stride=meta['raw_tile_bytes']; nq=meta['fixed_per_tile']; tiles=meta['tiles']
    assert len(data)==stride*tiles
    unique=set(); checksum=0
    for offset in range(0,len(data),stride):
        for i in range(nq): unique.add(data[offset+32*i:offset+32*(i+1)])
        for i in range(20): unique.add(data[offset+128+32*i:offset+160+32*i])
        expected=struct.unpack_from('<80h',data,offset+768)
        checksum+=sum(expected[q*20+t] for q in range(nq) for t in range(20))
    assert len(unique)>=100000, (group,'fewer than 100k distinct DNA sequences',len(unique))
    return dict(group=group.name,unique_DNA_sequences=len(unique),
                max_point_pairs=tiles*nq*20,exact_distance_checksum=checksum,
                input_sha256=hashlib.sha256(data).hexdigest())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);args=ap.parse_args()
    root=args.run.resolve();assert (root/'ALL_GROUPS_COMPLETE').exists()
    summaries=[];fits=[];points=[];inventory=[]
    for group in sorted(root.glob('k*_*')):
        if not group.is_dir():continue
        assert (group/'COMPLETE').exists(),group
        meta=json.loads((group/'dataset.json').read_text())
        inventory.append(input_inventory(group,meta))
        rows=list(csv.DictReader((group/'trials.csv').open()))
        audit_rejections(group,rows)
        for row in rows:
            stem=f"{row['algorithm']}_r{row['repeat']}_m{row['multiplier']}"
            audit_windows(group/(stem+'.windows.bin'),row)
        for algorithm in sorted({r['algorithm'] for r in rows}):
            ar=[r for r in rows if r['algorithm']==algorithm]; rr=[]
            for repeat in sorted({int(r['repeat']) for r in ar}):
                pr=sorted((r for r in ar if int(r['repeat'])==repeat),key=lambda r:int(r['N']))
                assert len(pr)==5
                x=[int(r['N']) for r in pr];y=[int(r['cpu_ns']) for r in pr]
                f=fit(x,y);wf=fit(x,[int(r['wall_ns']) for r in pr]);cf=fit(x,[int(r['user_cycles']) for r in pr])
                rr.append(f)
                fits.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,repeat=repeat,
                    slope_ns_per_pair=f['slope'],intercept_ns=f['intercept'],r2=f['r2'],
                    within_fit_stderr=f['stderr'],max_abs_residual_ns=f['max_abs_residual'],
                    wall_slope_ns_per_pair=wf['slope'],cycle_slope_per_pair=cf['slope']))
                for r in pr:
                    n=int(r['N']); predicted=f['intercept']+f['slope']*n
                    points.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,
                        repeat=repeat,N=n,cpu_ns=int(r['cpu_ns']),fitted_cpu_ns=predicted,
                        residual_ns=int(r['cpu_ns'])-predicted))
            slopes=[r['slope'] for r in rr];mean=st.mean(slopes)
            half=T95_DF4*st.stdev(slopes)/math.sqrt(len(slopes))
            summaries.append(dict(workload=meta['workload'],K=meta['K'],algorithm=algorithm,
                ns_per_pair=mean,ci95_half_ns=half,pair_per_second=1e9/mean,
                repeat_slope_min=min(slopes),repeat_slope_max=max(slopes),
                r2_min=min(f['r2'] for f in rr),r2_mean=st.mean(f['r2'] for f in rr),
                intercept_ns_mean=st.mean(f['intercept'] for f in rr),
                min_N=min(int(r['N']) for r in ar),max_N=max(int(r['N']) for r in ar),
                passing_fraction_min=min(int(r['passing_pairs'])/int(r['N']) for r in ar),
                passing_fraction_max=max(int(r['passing_pairs'])/int(r['N']) for r in ar),
                accepted_windows=sum(int(r['windows']) for r in ar),
                rejected_windows=sum(int(r['rejected_windows']) for r in ar),
                accepted_l2_misses=sum(int(r['accepted_l2_misses']) for r in ar),
                accepted_dram_fills=0))
    lookup={(r['workload'],r['K'],r['algorithm']):r for r in summaries}
    for r in summaries:
        key=(r['workload'],r['K'])
        for reference in ['edlib_original_hot','edlib_original_stack']:
            if key+(reference,) in lookup:
                r['speedup_vs_'+reference.removeprefix('edlib_')]=lookup[key+(reference,)]['ns_per_pair']/r['ns_per_pair']
    write_csv(root/'summary.csv',summaries);write_csv(root/'fits_per_repeat.csv',fits)
    write_csv(root/'fit_points.csv',points);write_csv(root/'input_inventory.csv',inventory)
    payload=dict(summary=summaries,points=points,inventory=inventory,
        independent_validation_pairs=131072,benchmark_oracle_pairs=sum(x['max_point_pairs'] for x in inventory),
        accepted_windows=sum(r['accepted_windows'] for r in summaries),
        rejected_windows=sum(r['rejected_windows'] for r in summaries))
    (root/'analysis.json').write_text(json.dumps(payload,indent=2)+'\n')
    (root/'ANALYSIS_AUDIT_PASS').write_text('All binary window totals, checksums, PMU ratios and input sizes audited.\n')
    for r in summaries:
        print(r['workload'],r['K'],r['algorithm'],f"{r['ns_per_pair']:.3f} +/- {r['ci95_half_ns']:.3f} ns/pair",f"R2_min={r['r2_min']:.6f}")

if __name__=='__main__':main()
