from pathlib import Path
import json,subprocess,time,sys
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path.insert(0,str(root/'tst/shocks-mc'))
from common import edit_input,provenance,sha
work=Path('/private/tmp/shock-highres128-long-20261014');work.mkdir(exist_ok=False)
binary=root/'build-release/src/athena';base=(root/'inputs/tests/shock_mc_showcase.athinput').read_text()
changes={'job/basename':'shock_highres128_long','mesh/nx1':128,'mesh/nx2':128,'mesh/nx3':128,
 'meshblock/nx1':128,'meshblock/nx2':128,'meshblock/nx3':128,'mesh/nghost':8,
 'time/tlim':1.,'time/ndiag':50,'shock_tracking/width':8,'shock_tracking/trace':'true',
 'shock_tracking/lineage':'showcase-m5-edot50-128-long','shock_tracking/directory':str(work/'shock_events'),
 'output1/dt':.02,'output2/dcycle':50,'output3/dcycle':-1,'output3/dt':.1,
 'output5/dcycle':-1,'output5/dt':.1,'output6/dcycle':-1,'output6/dt':.1,
 'output4/file_type':'rst','output4/dt':.25}
text=edit_input(base,changes).replace('# Native-forced Mach5 Fourier box, no cooling.','# Native-forced Mach5,128^3,long t=1.0 box; no cooling; full shock traces.')
(root/'inputs/tests/shock_mc_highres_long.athinput').write_text(text)
inputfile=work/'input.athinput';inputfile.write_text(text);command=[str(binary),'-i',str(inputfile)]
report={**provenance(binary),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'command':command,'input_sha256':sha(inputfile),
 'configuration':{'cells':128,'particles':8192,'tlim':1.,'width':8,'nghost':8,'trace':True,'gamma':1.4,'dedt':50,'initial_mach':5},'status':'running','started_unix':time.time()}
(work/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
start=time.monotonic()
with (work/'run.log').open('w') as f:code=subprocess.run(command,cwd=work,stdout=f,stderr=subprocess.STDOUT).returncode
report.update(status='completed' if code==0 else 'failed',exit=code,wall_seconds=time.monotonic()-start,log_sha256=sha(work/'run.log'))
(work/'execution.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
assert code==0
