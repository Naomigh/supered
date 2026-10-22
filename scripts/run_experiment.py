"""Reproducible single-process-at-a-time benchmark; all kernels pin CPU 24."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', type=Path)
    ap.add_argument('--repeats', type=int, default=5)
    ap.add_argument('--build-dir',type=Path,default=ROOT/'build_clang')
    ap.add_argument('--input-root',type=Path,default=ROOT.parent/'results/run_20261022')
    args = ap.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    build=args.build_dir.resolve()
    sources = {}
    for folder in ['include', 'src', 'tests', 'bench', 'scripts', 'generated']:
        for path in sorted((ROOT / folder).rglob('*')):
            if path.is_file() and '__pycache__' not in str(path):
                sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    sources['CMakeLists.txt']=hashlib.sha256((ROOT/'CMakeLists.txt').read_bytes()).hexdigest()
    for name in ['benchmark','validate','validate_indels','validate_banded','cache_control','profile','profile_v3']:
        sources['binary/'+name]=hashlib.sha256((build/name).read_bytes()).hexdigest()
    (out / 'source_sha256.json').write_text(json.dumps(sources, indent=2) + '\n')
    snapshot = out / 'source_snapshot'
    for folder in ['include', 'src', 'tests', 'bench', 'scripts', 'generated']:
        shutil.copytree(ROOT / folder, snapshot / folder, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(ROOT / 'CMakeLists.txt', snapshot)
    for name in ['benchmark','validate','validate_indels','validate_banded','cache_control','profile','profile_v3']:
        shutil.copy2(build/name,out/(name+'.executable'))
    shutil.copy2(build/'compile_commands.json',out/'compile_commands.json')
    for name in ['edlib_origin','leap_origin']:
        shutil.copytree(ROOT/name,snapshot/name,ignore=shutil.ignore_patterns('.git','__pycache__'))
    for name in ['validate','validate_indels','validate_banded','cache_control']:
        cmd=[str(build/name)]
        with (out/(name+'.log')).open('w') as f:
            subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    environment={}
    for name,cmd in [('cpu',['lscpu']),('compiler',['clang++-22','--version']),
                     ('processes',['ps','-eo','pid,pcpu,psr,args','--sort=-pcpu']),
                     ('edlib_commit',['git','-C',str(ROOT/'edlib_origin'),'rev-parse','HEAD']),
                     ('leap_commit',['git','-C',str(ROOT/'leap_origin'),'rev-parse','HEAD'])]:
        environment[name]=subprocess.check_output(cmd,text=True)
    (out/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    commands = []
    for k in [0, 2, 5, 50]:
        for workload in ['near', 'cross']:
            tiles = 131072 if (k == 0 and workload == 'near') else (32768 if workload == 'near' else 8192)
            group = out / f'k{k}_{workload}'
            cmd = [str(build/'benchmark'), str(group), str(k), workload,
                   str(tiles), str(args.repeats),str(args.input_root.resolve()/f'k{k}_{workload}'/'inputs.bin')]
            commands.append(cmd)
            (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
            print(datetime.datetime.now().isoformat(), 'START', k, workload, flush=True)
            with (out / f'k{k}_{workload}.log').open('w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            print(datetime.datetime.now().isoformat(), 'COMPLETE', k, workload, flush=True)
    (out / 'ALL_GROUPS_COMPLETE').write_text('8 groups complete.\n')

if __name__ == '__main__':
    main()
