import sys,json,subprocess,time
from pathlib import Path
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path[:0]=[str(root/'tst/shocks-mc'),str(root/'vis/python')]
from common import sha,provenance,edit_input
from bin_convert import read_binary
out=Path('/private/tmp/shock-highres128-cadence-smoke-fixed-20261014');out.mkdir(exist_ok=False)
inputfile=root/'inputs/tests/shock_mc_highres_long.athinput';binary=root/'build-release/src/athena'
variant=out/'input.athinput';variant.write_text(edit_input(inputfile.read_text(),{'time/nlim':2}))
command=[str(binary),'-i',str(variant),'time/tlim=.0004','output3/dt=.0001','output6/dt=.0001','output4/dt=-1','output5/dt=-1','shock_tracking/trace=false',f'shock_tracking/directory={out}/shock_events']
start=time.monotonic()
with (out/'run.log').open('w') as f:code=subprocess.run(command,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=90).returncode
assert code==0
cycles=[]
for path in sorted((out/'bin').glob('*hydro_w*.bin')):
 d=read_binary(path);cycles.append(d['cycle'])
assert sorted(set(cycles))==[0,1,2],cycles
report={**provenance(binary),'test':'Corrected input native timed-output smoke, two evolved steps only','exit':code,'command':command,'input_sha256':sha(inputfile),'cycles':cycles,'seconds':time.monotonic()-start,'log_sha256':sha(out/'run.log')}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
