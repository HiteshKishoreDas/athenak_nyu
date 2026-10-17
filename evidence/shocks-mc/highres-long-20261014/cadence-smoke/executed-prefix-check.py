from pathlib import Path
import csv,json,sys,numpy as np
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path[:0]=[str(root/'tst/shocks-mc'),str(root/'vis/python')]
from common import thermo
from shock_history import parse_observation
from bin_convert import read_binary
smoke=Path('/private/tmp/shock-highres128-cadence-smoke-fixed-20261014');run=Path('/private/tmp/shock-highres128-long-20261014');obs,times=thermo(smoke);assert set(times)=={0,2}
for g in (0,2):
 with (run/'shock_events'/f'trace{g}.rank0.csv').open() as f:
  for row in csv.reader(f):
   pid=int(row[0]);o=parse_observation(row[4:]);s=obs[g,pid];assert s['time']==o['time'] and list(s['x'])==o['x'] and s['fields']['density']==o['rho'] and np.isclose(s['fields']['pressure'],o['pressure'],rtol=1e-12,atol=1e-12)
x=read_binary(next((smoke/'bin').glob('*hydro_w.00000.bin')))['mb_data'];y=read_binary(next((run/'bin').glob('*hydro_w.00000.bin')))['mb_data'];assert x.keys()==y.keys() and all(np.array_equal(x[k],y[k]) for k in x)
report={'check':'Corrected-cadence smoke vs original long-run prefix','initial_fluid_fields_exact':True,'matched_thermo_samples':len(obs),'generations':[0,2],'time_and_positions_exact':True,'density_exact':True,'pressure_rtol':1e-12,'pressure_atol':1e-12}
(smoke/'prefix-equivalence.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
