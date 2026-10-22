"""Single-worker benchmark, sequential groups, immutable source/binary snapshots."""
from pathlib import Path
import argparse
import datetime
import fcntl
import hashlib
import json
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
KS=[0,2,5,10,20,30,40,50,60,70,80,90]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path)
    ap.add_argument('--windows',type=int,default=4096);ap.add_argument('--repeats',type=int,default=5)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    lock=open('/tmp/edlib_strong_measure.lock','a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    build=ROOT/'build_clang';manifest={}
    for folder in ['include','src','generated','bench','tests','scripts','examples']:
        for p in sorted((ROOT/folder).rglob('*')):
            if p.is_file() and '__pycache__' not in str(p):
                manifest[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        shutil.copytree(ROOT/folder,out/'source_snapshot'/folder,
            ignore=shutil.ignore_patterns('__pycache__'),dirs_exist_ok=True)
    for name in ['CMakeLists.txt','AGENTS.md','README.md','DESIGN_V3.md']:
        manifest[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        shutil.copy2(ROOT/name,out/'source_snapshot'/name)
    for name in ['benchmark_single_pair','cache_control','libsupered.a']:
        manifest['binary/'+name]=hashlib.sha256((build/name).read_bytes()).hexdigest()
        shutil.copy2(build/name,out/(name+'.artifact'))
    (out/'source_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
    shutil.copy2(build/'compile_commands.json',out/'compile_commands.json')
    environment={}
    for name,cmd in [('cpu',['lscpu']),('compiler',['clang++-22','--version']),
                     ('processes',['ps','-eo','pid,pcpu,psr,args','--sort=-pcpu'])]:
        environment[name]=subprocess.check_output(cmd,text=True)
    for repo in ['edlib_origin','leap_origin']:
        environment[repo+'_commit']=json.loads((ROOT/repo/'UPSTREAM.json').read_text())['commit']
        shutil.copytree(ROOT/repo,out/'source_snapshot'/repo,ignore=shutil.ignore_patterns('.git','__pycache__'))
    (out/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    with (out/'cache_control.log').open('w') as log:
        subprocess.run([str(build/'cache_control')],stdout=log,stderr=subprocess.STDOUT,check=True)
    commands=[]
    for k in KS:
        for workload in ['near','mixed']:
            group=f'k{k}_{workload}'
            cmd=[str(build/'benchmark_single_pair'),str(out/group),str(k),workload,str(args.windows),str(args.repeats)]
            commands.append(cmd);(out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            print(datetime.datetime.now().isoformat(),'START',group,flush=True)
            with (out/(group+'.log')).open('w') as log:
                p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
                time.sleep(.15)
                status=Path(f'/proc/{p.pid}/status')
                if status.exists():(out/(group+'.process_status')).write_text(status.read_text())
                if p.wait():raise RuntimeError('group failed: '+group)
            print(datetime.datetime.now().isoformat(),'COMPLETE',group,flush=True)
    (out/'ALL_GROUPS_COMPLETE').write_text('24 groups complete.\n')
    subprocess.run(['python3',str(ROOT/'scripts/analyze_single_pair.py'),str(out)],check=True)

if __name__=='__main__':main()
