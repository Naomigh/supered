"""Invoke under flock /tmp/edlib_strong_measure.lock; one pinned process at a time."""
from pathlib import Path
import subprocess,json,datetime,shutil
root=Path(__file__).resolve().parents[1]
out=root/'results/profile_final'
out.mkdir()
subprocess.run(['cmake','--build',str(root/'build_clang'),'--target','profile','-j4'],check=True,
               stdout=(out/'build.log').open('w'),stderr=subprocess.STDOUT)
binary=root/'build_clang/profile'
shutil.copy2(binary,out/'profile.executable')
shutil.copy2(root/'bench/profile.cpp',out/'profile.cpp')
commands=[]
for a in [0,1,16,17,18,19]:
    cmd=['perf','stat','-r','3','-e','{cycles:u,instructions:u,ex_no_retire.not_complete:u,ex_no_retire.load_not_complete:u}',
         '-o',str(out/f'a{a}.stat'),str(binary),str(a),'60000']
    commands.append(dict(start=datetime.datetime.now().isoformat(),cmd=cmd))
    (out/'commands.json').write_text(json.dumps(commands,indent=2))
    with (out/f'a{a}.time').open('w') as f:subprocess.run(cmd,stdout=f,check=True)
    commands[-1]['end']=datetime.datetime.now().isoformat()
    (out/'commands.json').write_text(json.dumps(commands,indent=2))
    print('PROFILE_COMPLETE',a,flush=True)
for a in [0,1,16]:
    cmd=['perf','record','-e','cycles:u','-c','2000000','-o',str(out/f'a{a}.perf.data'),'--',str(binary),str(a),'100000']
    with (out/f'a{a}.sample.time').open('w') as f,(out/f'a{a}.record.log').open('w') as err:
        subprocess.run(cmd,stdout=f,stderr=err,check=True)
    with (out/f'a{a}.annotate.txt').open('w') as f:
        subprocess.run(['perf','annotate','--stdio','-i',str(out/f'a{a}.perf.data')],stdout=f,check=True)
with (out/'assembly.txt').open('w') as f:
    subprocess.run(['objdump','-d','-C','-Mintel',str(binary)],stdout=f,check=True)
(out/'COMPLETE').write_text('Serial profiling completed, no benchmark overlapped.\n')
