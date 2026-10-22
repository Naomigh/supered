"""Benchmark the public router on the unchanged v3 cache-resident corpora."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/integration_20261022'
BUILD = ROOT / 'build_clang'
INPUT = ROOT.parent / 'work_v3/results/run_v3_20261022'

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for folder in ['include', 'src', 'generated', 'bench', 'tests', 'scripts', 'examples']:
        for path in sorted((ROOT/folder).rglob('*')):
            if path.is_file() and '__pycache__' not in str(path):
                manifest[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest['CMakeLists.txt'] = hashlib.sha256((ROOT/'CMakeLists.txt').read_bytes()).hexdigest()
    for name in ['benchmark', 'supered_example', 'supered_validate', 'libsupered.a']:
        path=BUILD/name
        manifest['binary/'+name]=hashlib.sha256(path.read_bytes()).hexdigest()
        shutil.copy2(path,OUT/(name+'.artifact'))
    (OUT/'source_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for folder in ['include','src','generated','bench','tests','scripts','examples']:
        shutil.copytree(ROOT/folder,OUT/'source_snapshot'/folder,
            ignore=shutil.ignore_patterns('__pycache__'),dirs_exist_ok=True)
    shutil.copy2(ROOT/'CMakeLists.txt',OUT/'source_snapshot')
    shutil.copy2(BUILD/'compile_commands.json',OUT/'compile_commands.json')
    environment={}
    for name,cmd in [('cpu',['lscpu']),('compiler',['clang++-22','--version']),
                     ('processes',['ps','-eo','pid,pcpu,psr,args','--sort=-pcpu'])]:
        environment[name]=subprocess.check_output(cmd,text=True)
    for repo in ['edlib_origin','leap_origin']:
        environment[repo+'_commit']=json.loads((ROOT/repo/'UPSTREAM.json').read_text())['commit']
        shutil.copytree(ROOT/repo,OUT/'source_snapshot'/repo,
            ignore=shutil.ignore_patterns('.git','__pycache__'))
    (OUT/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    commands=[]
    for k in [0,2,5,50]:
        for workload in ['near','cross']:
            tiles=131072 if (k==0 and workload=='near') else (32768 if workload=='near' else 8192)
            cmd=[str(BUILD/'benchmark'),str(OUT/f'k{k}_{workload}'),str(k),workload,
                 str(tiles),'5',str(INPUT/f'k{k}_{workload}'/'inputs.bin')]
            commands.append(cmd)
            (OUT/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            print(datetime.datetime.now().isoformat(),'START',k,workload,flush=True)
            with (OUT/f'k{k}_{workload}.log').open('w') as log:
                subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,
                    env=dict(os.environ,STRONG_ALGORITHMS='11'))
            print(datetime.datetime.now().isoformat(),'COMPLETE',k,workload,flush=True)
    (OUT/'ALL_GROUPS_COMPLETE').write_text('Eight public-API groups complete.\n')
    subprocess.run(['python3',str(ROOT/'scripts/analyze.py'),str(OUT)],check=True)

if __name__=='__main__':
    main()
