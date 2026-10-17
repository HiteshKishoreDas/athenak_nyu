#!/usr/bin/env python3
"""Frozen replay of the completed CSV/conserved-volume audit; read-only evidence."""
import collections,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'vis/python')
from bin_convert import read_binary
packet=Path('/private/tmp/shock-highres128-analysis-20261014');run=Path('/private/tmp/shock-highres128-long-20261014');report=json.loads((packet/'report.json').read_text());g=report['through_generation'];horizon=report['end_time'];n=8192
last=[-1]*n;complete=[0]*n;outcomes=collections.Counter();flags=collections.Counter();pressure=[];starts03=0;bounds=0;open_ids=set();bits={'left_censored':1,'invalid_observation':2,'merged_flicker':4,'no_plateau':8,'cadence_limited':16,'temperature_invalid':32,'right_censored':64,'domain_exit':128,'direction_ambiguous':256,'downstream_reentry':512}
with (packet/'episodes.csv').open() as f:
 for e in csv.DictReader(f):
  pid,seq,fl=int(e['particle_id']),int(e['sequence']),int(e['flags']);assert 0<=pid<n and seq==last[pid]+1;last[pid]=seq
  assert fl & ~1023==0;assert int(e['emission_generation'])<=g
  first,end=float(e['first_inside']),float(e['last_inside']);assert 0<=first<=end<=horizon
  starts03+=first<=.3;outcome=e['outcome'];outcomes[outcome]+=1
  for name,bit in bits.items():flags[name]+=bool(fl&bit)
  if e['pre_time']:assert 0<=float(e['pre_time'])<=first
  if e['post_time']:assert end<=float(e['post_time'])<=horizon
  if outcome=='complete':
   assert not fl&(1|64|128) and e['pre_time'] and e['post_time'];complete[pid]+=1
   lower=end-first;upper=float(e['post_time'])-float(e['pre_time']);assert 0<=lower<=upper;bounds+=1
   value=float(e['pressure_mach']);assert math.isfinite(value);pressure.append(value)
  elif outcome=='open':
   assert fl&64 and int(e['emission_generation'])==g and pid not in open_ids;open_ids.add(pid)
  else:assert outcome=='partial'
assert outcomes==collections.Counter({'complete':367074,'open':999});assert dict(outcomes)==report['outcomes'];assert report['counters']['partial']==0
hist=collections.Counter(complete);assert sum(hist.values())==n and sum(k*v for k,v in hist.items())==outcomes['complete'];assert {str(k):v for k,v in hist.items()}==report['completed_count_histogram']
assert dict(flags)==report['flags'];quantiles=np.quantile(pressure,[.1,.5,.9]);assert np.array_equal(quantiles,np.array(report['pressure_mach_quantiles']))
assert sum(outcomes.values())/n==report['starts_per_particle'];assert outcomes['complete']/n==report['complete_per_particle'];assert report['counters']['sampled']==n*(g+1)==49430528
baseline=Path('/private/tmp/supersonic-showcase-paper-20261013/episodes.csv');old=0
with baseline.open() as f:
 for e in csv.DictReader(f):old+=e['case']=='baseline' and float(e['first_inside'])<=.3
assert old==58225 and starts03==71419;assert old==report['common_window_comparison']['baseline64_starts'] and starts03==report['common_window_comparison']['highres128_starts']
p=run/'bin/shock_highres128_long.cons.00001.bin';d=read_binary(p);assert d['n_mbs']==1 and d['cycle']==g and d['time']==horizon
f={k:v[0] for k,v in d['mb_data'].items()};rho=f['dens'];assert rho.shape==(128,128,128) and np.all(rho>0)
momentum2=f['mom1']**2+f['mom2']**2+f['mom3']**2;gas=f['ener']-momentum2/(2*rho);pressure=.3999999999999999*gas;assert np.all(pressure>0)
local_mach2=momentum2/(1.4*pressure*rho);derived={'mass':float(np.mean(rho)),'mach_volume':float(np.sqrt(np.mean(local_mach2))),'supersonic_volume_fraction':float(np.mean(local_mach2>1)),'kinetic':float(np.mean(momentum2/(2*rho))),'gas_energy':float(np.mean(gas))}
reference=report['field_evolution'][-1];errors={}
for key,value in derived.items():
 errors[key]=abs(value-reference[key])
 if key=='supersonic_volume_fraction':assert value==reference[key]
 else:assert np.isclose(value,reference[key],rtol=1e-12,atol=1e-12),(key,value,reference[key])
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
result={'status':'passed_independent_csv_and_conserved_field_checks','scope':'Final exported episode stream and final native conserved volume only; no full trace scan or scientific/HPC sign-off','generation':g,'time':horizon,'particles_in_denominator':n,'outcomes':dict(outcomes),'unique_contiguous_sequences':True,'complete_uncensored_brackets_checked':bounds,'open_right_censored_identities':len(open_ids),'flags':dict(flags),'pressure_mach_quantiles':quantiles.tolist(),'common_window_starts':{'t':.3,'baseline64':old,'highres128':starts03},'derived_from_conserved':derived,'primitive_route_absolute_differences':errors,'float_rtol':1e-12,'float_atol':1e-12,'supersonic_fraction_comparison':'exact','source_sha256':{'report':digest(packet/'report.json'),'episodes':digest(packet/'episodes.csv'),'native_final_conserved':digest(p),'baseline_episodes':digest(baseline)}}
# The original second pass checked actual ID presence, beyond the histogram denominator.
ids=set()
with (packet/'episodes.csv').open() as f:
 for e in csv.DictReader(f):ids.add(int(e['particle_id']))
assert ids==set(range(8192))
result['all8192_particle_ids_present_in_episode_export']=True
# Preserve the original report; require byte-identical recomputation instead of rewriting.
original=Path('/private/tmp/highres128-final-independent-check.json')
assert digest(original)=='c9959485d131a3d251966ed77a069575cd31308753bb38580e8483a104e9c0cb'
assert json.dumps(result,indent=2)+'\n'==original.read_text()
print(json.dumps(result,indent=2))
print('PASS frozen computation reproduces original report bytes and SHA256:',digest(original))
