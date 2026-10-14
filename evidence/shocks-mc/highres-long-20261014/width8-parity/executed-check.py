from pathlib import Path
import json,sys,subprocess,io,csv,hashlib
import numpy as np
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path.insert(0,str(root/'tst/shocks'));from parity import REFERENCE_SHA256
refpath=root/'evidence/shocks/reference';sys.path.insert(0,str(refpath))
from shockit.fields import FluidCube
from shockit.config import ShockFinderConfig
from shockit.local_pass import analyze_cube_pass
out=Path('/private/tmp/shock-highres128-width8-20261014');out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name,h in REFERENCE_SHA256.items():assert sha(refpath/'shockit'/name)==h
binary=root/'build-release/shock_detector_test';records=[]
def check(name,rho,p,velocity,spacing,physical=None):
 gamma=1.4;energy=p/(gamma-1);pressure=(gamma-1)*energy;shape=rho.shape;file=out/(name+'.txt')
 array=np.column_stack([x.ravel() for x in (rho,*velocity,energy)])
 head=' '.join(map(str,shape))+'\n1.4 1.1 1000 1e-30 0 8 0 0 '+' '.join(map(str,spacing))
 np.savetxt(file,array,fmt='%.17g',header=head,comments='')
 proc=subprocess.run([str(binary),str(file)],capture_output=True,text=True,check=True);(out/(name+'.csv')).write_text(proc.stdout)
 header=next(l for l in proc.stdout.splitlines() if l.startswith('# i,j,k')).removeprefix('# ').strip().split(',')
 rows=np.loadtxt(io.StringIO(proc.stdout),delimiter=',',comments='#');native={k:rows[:,i].reshape(shape) for i,k in enumerate(header)}
 cfg=ShockFinderConfig(min_mach=1.1,mach_max=1000,shock_width_cells=8,grad_T_min=1e-30,div_v_threshold=0,reduce_to_centers=False)
 ref=analyze_cube_pass(FluidCube(rho,pressure,*velocity,*spacing,gamma=gamma),cfg)
 expected={'candidate':ref.shock_zone_mask,'inside':ref.full_shock_mask,'normal_x':ref.normals.normal_x,'normal_y':ref.normals.normal_y,'normal_z':ref.normals.normal_z,'rho_up':ref.jumps.density_up,'rho_down':ref.jumps.density_down,'pressure_up':ref.jumps.pressure_up,'pressure_down':ref.jumps.pressure_down,'temperature_up':ref.jumps.temperature_up,'temperature_down':ref.jumps.temperature_down,'mach_pressure':ref.mach_fields.mach_pressure,'mach_temperature':ref.mach_fields.mach_temperature,'divergence':ref.derived.div_v,'temperature':ref.derived.temperature}
 assert np.all(native['valid']==1)
 for key,value in expected.items():
  if key in ('candidate','inside'):assert np.array_equal(native[key],value),key
  else:np.testing.assert_allclose(native[key],value,rtol=1e-12,atol=1e-12,equal_nan=True)
 if physical is not None:
  selected=native['inside'].astype(bool);assert selected.any()
  np.testing.assert_allclose(native['mach_pressure'][selected],physical,rtol=1e-12,atol=1e-12)
  np.testing.assert_allclose(native['mach_temperature'][selected],physical,rtol=1e-12,atol=1e-12)
 records.append({'case':name,'status':'pass','cells':rho.size,'inside_cells':int(native['inside'].sum()),'fixture_sha256':sha(file),'native_csv_sha256':sha(out/(name+'.csv')),'physical_Mach':physical})
shape=(64,6,4);M=3.;gamma=1.4;ratio=(gamma+1)*M*M/((gamma-1)*M*M+2);q=(2*gamma*M*M-(gamma-1))/(gamma+1)
rho=np.ones(shape);p=np.ones(shape);vx=np.full(shape,M*np.sqrt(gamma));rho[32:]=ratio;p[32:]=q;vx[32:]/=ratio
check('NUM08_NUM10_width8_RH_M3',rho,p,[vx,np.zeros(shape),np.zeros(shape)],(1/128,1/128,1/128),physical=M)
rng=np.random.default_rng(72941);shape=(17,11,7);rho=np.exp(rng.normal(0,.2,shape));p=np.exp(rng.normal(0,.3,shape));velocity=[rng.normal(0,.2,shape) for _ in range(3)]
check('NUM08_width8_noncubic_seed72941',rho,p,velocity,(.017,.031,.047))
report={'status':'passed','width':8,'rtol':1e-12,'atol':1e-12,'masks':'exact','detector_sha256':sha(binary),'source_sha256':sha(root/'src/shocks/detector.hpp'),'reference_sha256':REFERENCE_SHA256,'cases':records,'scientific_signoff':'pending'}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
