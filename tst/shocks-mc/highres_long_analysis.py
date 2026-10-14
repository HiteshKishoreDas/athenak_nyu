"""Bounded analysis of the native 128^3 long run; no inferred crossing clocks."""
import argparse,csv,gc,html,json,re,sys,subprocess
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--trace',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--generation',type=int);p.add_argument('--pdf',action='store_true');a=p.parse_args()
root=a.root.resolve();run=a.run.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
sys.path[:0]=[str(root/'tst/shocks-mc'),str(root/'vis/python')]
from common import sha,read_thermo
from shock_history import manifests,read_chunk,read_suspended,parse_observation
from bin_convert import read_binary
from turbulence_pilot import FLAGS
N=8192;GAMMA=1.4;IDS=[4575,0,2487,4551]
execution=json.loads((run/'execution.json').read_text());progress=json.loads((a.trace/'progress.json').read_text())
g=a.generation if a.generation is not None else progress['last_generation']
assert g<=progress['last_generation']
final=a.generation is None
if final:assert execution['status']=='completed' and execution['exit']==0 and progress['status']=='completed' and progress['last_time']==1
metadata=manifests(run/'shock_events');assert list(metadata)==[0];lineage,contract=metadata[0]
# Stream events, keeping bounded summaries rather than nested observations for every event.
last=[-1]*N;counts=[0]*N;flags=Counter();status=Counter();starts=[];mp=[];mt=[];bounds=[];entry=[];exitw=[];selected={pid:[] for pid in IDS}
keys=['outcome','particle_id','sequence','flags','emission_generation','pre_time','first_inside','last_inside','post_time','pressure_mach','temperature_mach']
with (out/'episodes.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(keys)
 def consume(e,outcome):
  pid=e['particle_id'];assert 0<=pid<N and e['sequence']==last[pid]+1
  last[pid]=e['sequence'];status[outcome]+=1;starts.append(e['first']['time'])
  for name,bit in FLAGS.items():flags[name]+=bool(e['flags']&bit)
  first,end=e['first'],e['last'];pre=e['pre'];post=e['post']
  if outcome=='complete':
   assert pre['valid'] and post['valid'] and not pre['inside'] and not post['inside']
   counts[pid]+=1;lo=end['time']-first['time'];hi=post['time']-pre['time'];assert 0<=lo<=hi
   bounds.append((lo,hi));entry.append(first['time']-pre['time']);exitw.append(post['time']-end['time']);mp.append(first['mach_pressure']);mt.append(first['mach_temperature'])
  row=[outcome,pid,e['sequence'],e['flags'],e['emission_generation'],pre['time'] if pre['valid'] else '',first['time'],end['time'],post['time'] if post['valid'] else '',first['mach_pressure'],first['mach_temperature']]
  w.writerow(row)
  if pid in selected:selected[pid].append(e)
 for cycle in range(g+1):
  for e in read_chunk(run/'shock_events'/f'cycle{cycle}.rank0.csv',lineage,contract,cycle):consume(e,'complete' if e['status']==0 else 'partial')
 view=run if final else run/'checkpoint-views'/f'generation-{g}'
 for e in read_suspended(view/'shock_events',g):consume(e,'open')
trace_rows={pid:[] for pid in IDS}
with (a.trace/'selected-traces.csv').open() as f:
 for row in csv.DictReader(f):
  if int(row['generation'])<=g:trace_rows[int(row['particle_id'])].append(row)
assert all(len(rows)==g+1 for rows in trace_rows.values())
end_time=float(trace_rows[IDS[0]][-1]['time'])
frames=[]
with (a.trace/'frames.jsonl').open() as f:
 for line in f:
  row=json.loads(line)
  if row['generation']<=g:frames.append(row)
assert len(frames)==g+1 and frames[-1]['time']==end_time
logs=(run if final else view)/'run.log';log=logs.read_text();matches=re.findall(r'^Shock rank=0 (.+)$',log,re.M);assert len(matches)==1
counters={k:float(v) if 'seconds' in k else int(v) for k,v in re.findall(r'(\w+)=([^ ]+)',matches[0])}
assert counters['sampled']==N*(g+1) and counters['invalid']==0
assert counters['complete']==status['complete'] and counters['partial']==status['partial'] and counters['suspended']==status['open']
paths=list((run/'bin').glob('*hydro_w*.bin'))
for folder in sorted((run/'checkpoint-views').glob('generation-*')):paths+=list((folder/'bin').glob('*hydro_w*.bin'))
snapshots={};mask_records=[]
for path in paths:
 d=read_binary(path);cycle=d['cycle']
 if cycle>g:continue
 f={k:v[0] for k,v in d['mb_data'].items()};assert d['n_mbs']==1
 if cycle in snapshots:
  assert all(np.array_equal(f[k],snapshots[cycle]['fields'][k]) for k in f);continue
 rho=f['dens'];pressure=(GAMMA-1)*f['eint'];assert rho.shape==(128,128,128) and all(np.isfinite(v).all() for v in f.values()) and np.all(rho>0) and np.all(pressure>0)
 v2=sum(f[k]**2 for k in ('velx','vely','velz'));mach=np.sqrt(v2/(GAMMA*pressure/rho));t=frames[cycle]['time'];assert float(format(t,'.6g'))==d['time']
 stats={'generation':cycle,'time':t,'mass':float(rho.mean()),'mach_volume':float(np.sqrt((mach**2).mean())),'mach_mass':float(np.sqrt((rho*mach**2).sum()/rho.sum())),'supersonic_volume_fraction':float((mach>1).mean()),'supersonic_mass_fraction':float(rho[mach>1].sum()/rho.sum()),'kinetic':float((rho*v2/2).mean()),'gas_energy':float(f['eint'].mean())}
 assert abs(stats['mass']-1)<=1e-12
 actual={}
 with (run/'shock_events'/f'trace{cycle}.rank0.csv').open() as stream:
  for row in csv.reader(stream):
   pid=int(row[0]);o=parse_observation(row[4:]);i,j,k=np.floor(np.array(o['x'])*128).astype(int)
   assert o['rho']==rho[k,j,i] and abs(o['pressure']-pressure[k,j,i])<=1e-12*max(1,o['pressure'])
   actual[pid]=o
 assert set(actual)==set(range(N))
 if path.parent.parent!=run:
  thp=list((path.parent.parent/'prtcl_thermo_history').glob('*.thp'));assert len(thp)==1
  obs,clocks=read_thermo(thp[0],allow_identical_replay=True);assert clocks=={cycle:t} and len(obs)==N
  for pid,o in actual.items():
   native=obs[cycle,pid];assert native['time']==t and list(native['x'])==o['x'] and native['fields']['density']==o['rho']
   assert np.isclose(native['fields']['pressure'],o['pressure'],rtol=1e-12,atol=1e-12)
 # Native detector on a crop with enough real support for this plane only.
 # Central z=64 uses gradients +/-1 and endpoint samples +/-8: z55..73 suffices.
 fixture=out/f'slice-fixture-{cycle}.txt';native_csv=out/f'slice-detector-{cycle}.csv';cache=out/f'slice-{cycle}.npz'
 if cache.exists() and 'rho_up' in np.load(cache):
  mask=dict(np.load(cache));record=json.loads((out/f'slice-{cycle}.json').read_text());assert record['fluid_sha256']==sha(path) and record['detector_executable_sha256']==sha(root/'build-release/shock_detector_test')
 else:
  arrays=[f[k][55:74].transpose(2,1,0).ravel() for k in ('dens','velx','vely','velz','eint')]
  header='128 128 19\n1.4 1.1 1000 1e-30 0 8 0 0 '+ ' '.join([format(1/128,'.17g')]*3)
  np.savetxt(fixture,np.column_stack(arrays),fmt='%.17g',header=header,comments='')
  with native_csv.open('w') as stream:subprocess.run([str(root/'build-release/shock_detector_test'),str(fixture)],stdout=stream,check=True,timeout=180)
  with native_csv.open() as stream:names=next(line for line in stream if line.startswith('# i,j,k')).removeprefix('# ').strip().split(',')
  needed=['valid','inside','mach_pressure','mach_temperature','rho_up','rho_down','pressure_up','pressure_down','temperature_up','temperature_down']+[s for s in names if s.startswith('normal')]
  data=np.loadtxt(native_csv,delimiter=',',comments='#',usecols=[names.index(k) for k in needed]);mask={k:data[:,n].reshape(128,128,19)[:,:,9].T for n,k in enumerate(needed)}
  assert np.all(mask['valid']==1);np.savez_compressed(cache,**mask)
  record={'generation':cycle,'time':t,'fluid_sha256':sha(path),'fixture_sha256':sha(fixture),'native_csv_sha256':sha(native_csv),'detector_executable_sha256':sha(root/'build-release/shock_detector_test'),'scope':'Only z=64 plane: crop z55..73 includes actual +/-8 endpoints and +/-1 gradients. Full x/y periodic neighbors; artificial z wrap cannot affect central plane. No whole-volume mask assertion.'}
  (out/f'slice-{cycle}.json').write_text(json.dumps(record,indent=2)+'\n')
 checks=0
 for o in actual.values():
  i,j,k=np.floor(np.array(o['x'])*128).astype(int)
  if k==64:
   assert o['inside']==bool(mask['inside'][j,i])
   for key in ('mach_pressure','mach_temperature'):
    assert np.isclose(o[key],mask[key][j,i],rtol=1e-12,atol=1e-12,equal_nan=True)
   for dim,key in enumerate(('normal_x','normal_y','normal_z')):assert np.isclose(o['normal'][dim],mask[key][j,i],rtol=1e-12,atol=1e-12)
   for side,suffix in (('upstream','up'),('downstream','down')):
    for dim,key in enumerate(('rho','pressure','temperature')):assert np.isclose(o[side][dim],mask[key+'_'+suffix][j,i],rtol=1e-12,atol=1e-12)
   checks+=1
 record['actual_trace_particles_checked']=checks;(out/f'slice-{cycle}.json').write_text(json.dumps(record,indent=2)+'\n');mask_records.append(record)
 snapshots[cycle]={'fields':f,'mach':mach,'mask':mask,'stats':stats};del actual;gc.collect()
assert 0 in snapshots and abs(snapshots[0]['stats']['mach_volume']-5)<=1e-12
report={'status':'completed_scientific_review_pending' if final else 'checkpoint_preview_run_still_running','through_generation':g,'end_time':end_time,'particles':N,'configuration':execution['configuration'],'execution':execution,'trace_validation':progress,'counters':counters,'outcomes':dict(status),'starts_per_particle':len(starts)/N,'complete_per_particle':status['complete']/N,'completed_count_histogram':dict(sorted(Counter(counts).items())),'flags':dict(flags),'mask_slice_validation':mask_records,'field_evolution':[snapshots[k]['stats'] for k in sorted(snapshots)],'scientific_signoff':'pending','limitations':['Adiabatic evolving forced realization; no stationarity or physical shock calibration.','Inherited native RK2 forcing thermal oracle fails; forcing source unchanged.','Endpoint no_plateau/cadence flags retained. MC cell hops are not continuous material trajectories.','Width8 at128 matches width4 at64 physical half-span1/16; resolution/CFL/MC dynamics still coupled.','Current run dcycle=-1 presence suppresses intermediate volume outputs; zero-step native checkpoint recovery supplies quarter-time maps. Corrected future input removes dcycle keys.','Local Serial double CPU run, no new MPI/GPU/Rusty/Frontier qualification.'],'figures':[],'source_sha256':{str(Path(__file__).name):sha(__file__),'input.athinput':sha(run/'input.athinput')}}
with (out/'evolution.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(report['field_evolution'][0]));w.writeheader();w.writerows(report['field_evolution'])
(out/'frames.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in frames))
with (out/'selected-traces.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(trace_rows[IDS[0]][0]));w.writeheader()
 for cycle in range(g+1):w.writerows(trace_rows[pid][cycle] for pid in IDS)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
def save(fig,name,title,caption):
 fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180);fig.savefig(out/(name+'.svg'));
 if a.pdf:fig.savefig(out/(name+'.pdf'))
 plt.close(fig);report['figures'].append({'name':name,'title':title,'caption':caption})
def cdf(ax,x,label):
 x=np.sort(np.asarray(x));assert len(x);ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=label)
selected_snaps=[snapshots[k] for k in sorted(snapshots) if k>0]
fig,axes=plt.subplots(2,len(selected_snaps),figsize=(max(6,3*len(selected_snaps)),5.4),squeeze=False)
for col,s in enumerate(selected_snaps):
 for row,values in enumerate((np.log10(s['fields']['dens'][64]),s['mach'][64])):
  im=axes[row,col].imshow(values,origin='lower',extent=(0,1,0,1),cmap='magma' if row==0 else 'viridis',vmin=-1 if row==0 else 0,vmax=1 if row==0 else max(5,float(values.max())))
  mask=s['mask']['inside']
  if mask.min()!=mask.max():axes[row,col].contour((np.arange(128)+.5)/128,(np.arange(128)+.5)/128,mask,levels=[.5],colors='cyan',linewidths=.45)
  axes[row,col].set_title(f't={s["stats"]["time"]:.4f}');axes[row,col].set_xlabel('x / box length');fig.colorbar(im,ax=axes[row,col],label='log10 density' if row==0 else 'Local flow Mach',shrink=.8)
axes[0,0].set_ylabel('y / box length');axes[1,0].set_ylabel('y / box length')
save(fig,'01_flow','128-cubed flow and detected regions','Actual z=64 cell-center slices (z=0.50390625). Top: log10 density; bottom: local speed/sound speed. Cyan outlines use the native unreduced width8 detector. The central crop has full real support; masks and finite/NaN strength values agree with actual live particle samples on this plane at rtol=atol=1e-12. Checkpoint maps were recovered with the same binary and zero evolved steps. These detector masks are numerical evidence, not independently certified physical shock surfaces.')
fig,ax=plt.subplots(1,3,figsize=(10,3.4));evo=report['field_evolution'];ts=[d['time'] for d in evo]
ax[0].plot(ts,[d['mach_volume'] for d in evo],'o-',label='Volume rms local Mach');ax[0].axhline(1,color='k',ls=':');ax[0].legend(fontsize=7);ax[1].plot(ts,[d['supersonic_volume_fraction'] for d in evo],'o-',label='Volume fraction');ax[1].plot(ts,[d['supersonic_mass_fraction'] for d in evo],'s-',label='Mass fraction');ax[1].legend(fontsize=7);ax[2].plot(ts,[d['kinetic'] for d in evo],'o-',label='Kinetic');ax[2].plot(ts,[d['gas_energy'] for d in evo],'s-',label='Internal');ax[2].legend(fontsize=7)
for x in ax:x.set_xlabel('Code time')
ax[0].set_ylabel('Mach');ax[1].set_ylabel('Fraction above local Mach 1');ax[2].set_ylabel('Energy / box volume')
save(fig,'02_evolution','Flow evolution under continuous forcing','Volume rms Mach = sqrt(mean(v^2/(gamma P/rho))). The gas is adiabatic: initial Mach5 and continuous edot50 do not hold its sound speed or Mach number fixed. Lines only connect the saved checkpoint states. Actual mass stays within1e-12 of1. Forcing energy is not independently calibrated because the inherited RK2 source thermal oracle fails; no stationary regime is asserted.')
fig,ax=plt.subplots(1,3,figsize=(10,3.4));hist=Counter(counts);ax[0].bar(sorted(hist),[hist[x]/N for x in sorted(hist)]);ss=np.sort(starts);ax[1].step(ss,np.arange(1,len(ss)+1)/N,where='post');bins=np.linspace(0,end_time,21);h,_=np.histogram(ss,bins);ax[2].stairs(h/(N*np.diff(bins)),bins)
ax[0].set_xlabel('Complete episodes per particle');ax[0].set_ylabel('Fraction of all8192 particles');ax[1].set_xlabel('First-inside recording time');ax[1].set_ylabel('Cumulative starts per particle');ax[2].set_xlabel('First-inside recording time');ax[2].set_ylabel('Starts / particle / code time')
save(fig,'03_encounters','Accumulated particle encounter statistics','All8192 IDs are included, including zero-complete particles. Starts count every committed or terminal-open identity once at its first-inside sample. Partial/open outcomes are retained separately. Rate bins cover the stated time window; these are numerical recording rates, not exact physical shock crossing rates or independent Poisson trials.')
fig,ax=plt.subplots(2,2,figsize=(10,6));times=np.array([d['time'] for d in frames]);ax[0,0].plot(times,[d['inside']/N for d in frames],label='Raw inside mask');ax[0,0].plot(times,[(d['phase_active']+d['phase_pending'])/N for d in frames],label='Active + pending episode');ax[0,0].legend(fontsize=8);ax[0,0].set_ylabel('Fraction of particles');ax[0,0].set_xlabel('Code time');ax[0,1].barh(list(FLAGS),[flags[k]/len(starts) for k in FLAGS]);ax[0,1].set_xlabel('Fraction with flag; overlaps allowed')
finite=np.isfinite(mp);cdf(ax[1,0],np.array(mp)[finite],'Pressure Mach');valid=np.isfinite(mp)&np.isfinite(mt);hh=ax[1,1].hexbin(np.array(mp)[valid],np.array(mt)[valid],gridsize=35,bins='log',mincnt=1);fig.colorbar(hh,ax=ax[1,1],label='Complete episodes per bin');ax[1,0].set_ylabel('Completed-event CDF');ax[1,0].set_xlabel('First-inside pressure Mach');ax[1,1].set_xlabel('Pressure Mach');ax[1,1].set_ylabel('Temperature Mach')
save(fig,'04_membership_strength','Shock flags, episode state and strength evidence','The full per-step trace explicitly saves each raw inside flag. Raw membership and encounter phase differ: a pending exit can be outside while the episode remains open, and a suppressed downstream reentry can be inside while its phase is outside. Both are plotted from recorded states. All validity bits are retained. Finite first-inside strength estimates compare two native estimators; their agreement is not a physical calibration. Either estimator can accept membership, so a pressure estimate belowmin_mach1.1 can remain present.')
fig,ax=plt.subplots(3,1,figsize=(10,7),sharex=True)
for row,pid in enumerate(IDS):
 for e in selected[pid]:
  first,end=e['first']['time'],e['last']['time'];ax[0].plot([first,end],[row,row],lw=5,color='C0')
  if e['pre']['valid']:ax[0].plot([e['pre']['time'],first],[row,row],lw=2,color='C1')
  if e['post']['valid']:ax[0].plot([end,e['post']['time']],[row,row],lw=2,color='C1')
  if e['flags']&64:ax[0].plot(end,row,'>',color='C3')
 rr=trace_rows[pid];t=np.array([float(d['time']) for d in rr]);inside=np.array([int(d['inside']) for d in rr]);ax[1].step(t,row+.6*inside,where='post',label=f'ID{pid}');ax[2].plot(t,[float(d['rho']) for d in rr],label=f'ID{pid}')
ax[0].set_yticks(range(4),[f'ID {pid}: {counts[pid]} complete' for pid in IDS]);ax[1].set_yticks(range(4),[f'ID {pid}: 0/1' for pid in IDS]);ax[1].set_ylabel('Saved raw inside flag');ax[2].set_yscale('log');ax[2].set_ylabel('Sampled particle density');ax[2].set_xlabel('Code time');ax[2].legend(ncol=4,fontsize=8)
save(fig,'05_particle_histories','Actual per-step particle flags and encounter histories','Same IDs as the earlier64-cubed history illustration, with fresh128-cubed realizations. Top: blue first-to-last inside evidence, orange entry/exit brackets, red right-censored terminal marker. Middle: every saved raw0/1 inside flag, vertically offset by particle. Bottom: actual every-step density samples. Connecting samples guides the eye; it does not invent continuous trajectories or exact crossing times. Selected IDs are illustrations, not unbiased ensemble statistics.')
fig,ax=plt.subplots(1,2,figsize=(10,3.8));bd=np.array(bounds);cdf(ax[0],bd[:,0],'Lower span bound');cdf(ax[0],bd[:,1],'Upper span bound');cdf(ax[1],entry,'Entry bracket');cdf(ax[1],exitw,'Exit bracket')
for x in ax:x.set_xlabel('Code time');x.set_ylabel('Completed-event CDF');x.legend()
save(fig,'06_timing','Timing bounds from observed samples','For complete episodes the entry bracket is[pre,first], exit bracket[last,post], and elapsed-span interval[last-first,post-pre]. No midpoint is assigned. Merged flicker may include brief outside intervals, so this span is not continuous inside residence time. Terminal-open records receive no artificial completed duration.')
old=[]
with (a.baseline/'episodes.csv').open() as f:
 for row in csv.DictReader(f):
  if row['case']=='baseline':old.append(float(row['first_inside']))
common_end=min(.3,end_time);new=np.array(starts);new=new[new<=common_end];old=np.array(old);old=old[old<=common_end]
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
for data,label in ((old,'64³, width4'),(new,'128³, width8')):
 data=np.sort(data);ax[0].step(data,np.arange(1,len(data)+1)/N,where='post',label=label);bins=np.linspace(0,common_end,13);h,_=np.histogram(data,bins);ax[1].stairs(h/(N*np.diff(bins)),bins,label=label)
for x in ax:x.set_xlabel('First-inside recording time');x.legend()
ax[0].set_ylabel('Cumulative recorded starts per particle');ax[1].set_ylabel('Starts / particle / code time')
report['common_window_comparison']={'time':common_end,'baseline64_starts':len(old),'highres128_starts':len(new),'physical_half_span':1/16,'claim':'Sensitivity only; not a convergence gate. Closure counts are not compared across unequal endpoints.'}
save(fig,'07_resolution','Matched detector span in the common time window','Both cases use8192 particles, CFL.30 and the same seeds/forcing inputs. Width8 at128³ matches width4 at64³ physical half-span1/16box. Only first-inside starts through the common observed endpoint are compared, including unique terminal-open identities; unequal later closure times do not bias this measure. Resolution still changes fluid evolution, MC random draws and observation cadence together. This is measured sensitivity, with no claimed converged tolerance or independent ensemble error bar.')
report['pressure_mach_quantiles']=np.quantile(np.array(mp)[np.isfinite(mp)],[.1,.5,.9]).tolist()
(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
text='<h1>128-cubed forced box: longer particle histories</h1><p>'+html.escape(report['status'])+f'; through t={end_time:.17g}, generation{g}. Scientific approval pending.</p>'
for item in report['figures']:text+='<article><h2>'+html.escape(item['title'])+'</h2><img src="'+item['name']+'.png"><p>'+html.escape(item['caption'])+'</p><a href="'+item['name']+'.svg">Vector SVG</a></article>'
(out/'index.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:17px/1.5 system-ui;max-width:1200px;margin:2em auto}img{width:100%}article{margin:3em 0}</style>'+text)
print(json.dumps({'status':report['status'],'generation':g,'time':end_time,'outcomes':dict(status),'output':str(out)}))
