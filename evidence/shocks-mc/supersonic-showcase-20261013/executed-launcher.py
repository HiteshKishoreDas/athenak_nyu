import json,subprocess,time
from pathlib import Path
import sys
r=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path.insert(0,str(r/'tst/shocks-mc'))
from common import edit_input,provenance,sha
w=Path('/private/tmp/supersonic-showcase-20261013');w.mkdir(exist_ok=False)
binary=r/'build-release/src/athena'
changes={'job/basename':'supersonic_showcase','mesh/nx1':64,'mesh/nx2':64,'mesh/nx3':64,
 'meshblock/nx1':64,'meshblock/nx2':64,'meshblock/nx3':64,'time/tlim':.3,
 'initial_perturbations/velocity_rms':5,'tracer_seed_initial/count_per_event':8192,
 'shock_tracking/event_capacity':32768,'shock_tracking/lineage':'showcase-m5-edot50-64',
 'turb_driving/normalization':'edot','turb_driving/dedt':50,
 'output3/variable':'hydro_w','output3/dt':-1,'output3/dcycle':50,
 'output5/file_type':'bin','output5/variable':'turb_force','output5/id':'force','output5/dt':-1,'output5/dcycle':50,'output5/data_precision':'real',
 'output6/file_type':'bin','output6/variable':'hydro_u','output6/id':'cons','output6/dt':-1,'output6/dcycle':50,'output6/data_precision':'real'}
base=edit_input((r/'inputs/tests/shock_mc_forced.athinput').read_text(),changes)
base='\n'.join(line for line in base.splitlines() if not line.startswith('accel_rms ='))+'\n'
(r/'inputs/tests/shock_mc_showcase.athinput').write_text(base)
report={**provenance(binary),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=r,text=True).strip(),'runs':[],'status':'in_progress'}
for name,overrides in [('baseline',{}),('tracking_off',{'shock_tracking/enabled':'false'}),('half_cfl',{'time/cfl_number':.15,'shock_tracking/lineage':'showcase-m5-edot50-64-cfl015'}),('coarse32',{'mesh/nx1':32,'mesh/nx2':32,'mesh/nx3':32,'meshblock/nx1':32,'meshblock/nx2':32,'meshblock/nx3':32,'shock_tracking/lineage':'showcase-m5-edot50-32'})]:
 d=w/name;d.mkdir();overrides['shock_tracking/directory']=str(d/'shock_events');i=d/'input.athinput';i.write_text(edit_input(base,overrides))
 command=[str(binary),'-i',str(i)];t=time.monotonic()
 with (d/'run.log').open('w') as f:
  try:code=subprocess.run(command,cwd=d,stdout=f,stderr=subprocess.STDOUT,timeout=600).returncode
  except subprocess.TimeoutExpired:code='timeout'
 row={'name':name,'command':command,'exit':code,'wall_seconds':time.monotonic()-t,'input_sha256':sha(i),'log_sha256':sha(d/'run.log')};report['runs'].append(row)
 (w/'execution.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(row),flush=True)
 assert code==0,row
report['status']='simulations_completed';(w/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
