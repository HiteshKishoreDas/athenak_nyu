"""Reproducible local forced-box showcase; scientific calibration remains pending."""
import argparse
import csv
from collections import Counter
import gc
import html
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from common import ROOT, edit_input, provenance, sha, thermo
from turbulence_pilot import FLAGS
sys.path.insert(0,str(ROOT/'vis/python'))
from shock_history import read_history, read_suspended
from bin_convert import read_binary
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES={'baseline':{},'tracking_off':{'shock_tracking/enabled':'false'},
       'half_cfl':{'time/cfl_number':.15,'shock_tracking/lineage':'showcase-m5-edot50-64-cfl015'},
       'coarse32':{'shock_tracking/lineage':'showcase-m5-edot50-32'},
       'coarse32_matched':{'shock_tracking/width':2,'shock_tracking/lineage':'showcase-m5-edot50-32-w2'}}
LABELS={'baseline':'64³, CFL 0.30, width 4','half_cfl':'64³, CFL 0.15, width 4',
        'coarse32':'32³, width 4','coarse32_matched':'32³, width 2 (matched span)'}
N=8192;GAMMA=1.4;END=.3


def launch(binary,work):
    work.mkdir(parents=True,exist_ok=False)
    execution={**provenance(binary),'checkout_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'runs':[],'status':'in_progress'}
    base=(ROOT/'inputs/tests/shock_mc_showcase.athinput').read_text()
    for name,configuration in CASES.items():
        d=work/name;d.mkdir();changes=dict(configuration)
        if name.startswith('coarse32'):
            changes.update({section+'/'+axis:32 for section in ('mesh','meshblock') for axis in ('nx1','nx2','nx3')})
        changes['shock_tracking/directory']=str(d/'shock_events')
        inp=d/'input.athinput';inp.write_text(edit_input(base,changes));command=[str(binary),'-i',str(inp)]
        start=time.monotonic()
        with (d/'run.log').open('w') as log:
            try:code=subprocess.run(command,cwd=d,stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
            except subprocess.TimeoutExpired:code='timeout'
        row={'name':name,'command':command,'exit':code,'wall_seconds':time.monotonic()-start,'input_sha256':sha(inp),'log_sha256':sha(d/'run.log')}
        execution['runs'].append(row);(work/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
        assert code==0,row
    execution['status']='simulations_completed';(work/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')


def fields(directory):
    # Existing native binary reader; all arrays are double and finite.
    snapshots={}
    for p in sorted((directory/'bin').glob('*hydro_w*.bin')):
        d=read_binary(p);assert d['n_mbs']==1
        f={k:v[0] for k,v in d['mb_data'].items()};key=d['cycle']
        if key in snapshots:
            assert all(np.array_equal(f[k],snapshots[key]['fields'][k]) for k in f);continue
        rho=f['dens'];pressure=(GAMMA-1)*f['eint'];n=rho.shape[0]
        assert rho.shape==(n,n,n) and rho.dtype==np.float64 and np.all(rho>0) and np.all(pressure>0)
        assert all(np.isfinite(v).all() for v in f.values())
        velocity=[f[k] for k in ('velx','vely','velz')];v2=sum(v*v for v in velocity);m2=v2/(GAMMA*pressure/rho)
        div=sum((np.roll(v,-1,axis=axis)-np.roll(v,1,axis=axis))*(n/2) for v,axis in zip(velocity,(2,1,0)))
        snapshots[key]={'cycle':key,'metadata_time':d['time'],'path':p,'fields':f,'n':n,'rho':rho,'pressure':pressure,'mach_local':np.sqrt(m2),'div':div,
                        'mass':float(rho.mean()),'mach_volume':float(np.sqrt(m2.mean())),
                        'mach_mass':float(np.sqrt((rho*m2).sum()/rho.sum())),
                        'supersonic_volume_fraction':float((m2>1).mean()),
                        'supersonic_mass_fraction':float(rho[m2>1].sum()/rho.sum()),
                        'kinetic':float((rho*v2/2).mean()),'gas_energy':float(f['eint'].mean())}
    assert snapshots and min(snapshots)==0
    return snapshots


def analyze(directory,snapshots):
    events=read_history(directory/'shock_events');open_events=read_suspended(directory/'shock_events');all_events=events+open_events
    assert len({e['identity'] for e in all_events})==len(all_events)
    by_particle={pid:[] for pid in range(N)}
    for e in all_events:
        assert e['particle_id'] in by_particle;by_particle[e['particle_id']].append(e)
    duration=[];waiting=[]
    for group in by_particle.values():
        group.sort(key=lambda e:e['sequence']);assert [e['sequence'] for e in group]==list(range(len(group)))
        for e in group:
            if e['status']==0:
                assert e['pre']['valid'] and e['post']['valid'] and not e['flags']&(1|64|128)
                low=e['last']['time']-e['first']['time'];high=e['post']['time']-e['pre']['time'];assert 0<=low<=high
                duration.append((low,high))
        for a,b in zip(group,group[1:]):
            if a['status']==0 and b['status']==0:
                assert b['sequence']==a['sequence']+1
                low=max(0,b['pre']['time']-a['post']['time']);high=max(0,b['first']['time']-a['last']['time']);assert low<=high
                waiting.append((low,high))
    counts=Counter(e['particle_id'] for e in events)
    ranked=sorted(range(N),key=lambda p:(counts[p],p));representative_ids=[ranked[int(q*(N-1))] for q in (0,.25,.5,.75,.9,1)]
    observations,times=thermo(directory);populations={}
    for g,pid in observations:populations.setdefault(g,set()).add(pid)
    assert all(pids==set(range(N)) for pids in populations.values())
    for g,d in snapshots.items():
        d['time']=times[g];assert float(format(times[g],'.6g'))==d['metadata_time']
    assert times[max(times)]==END
    selected={pid:[o for (g,p),o in observations.items() if p==pid] for pid in representative_ids}
    final_positions=np.array([observations[max(times),pid]['x'] for pid in range(N)])
    for pid,o in enumerate(final_positions):
        f=snapshots[max(snapshots)]['fields'];n=f['dens'].shape[0];i,j,k=np.floor(o*n).astype(int)
        native=observations[max(times),pid]['fields']
        assert native['density']==f['dens'][k,j,i]
        assert abs(native['pressure']-(GAMMA-1)*f['eint'][k,j,i])<=1e-12*max(1,native['pressure'])
    del observations;gc.collect()
    log=(directory/'run.log').read_text();line=re.findall(r'^Shock rank=0 (.+)$',log,re.M);assert len(line)==1
    counters={k:float(v) if 'seconds' in k else int(v) for k,v in re.findall(r'(\w+)=([^ ]+)',line[0])}
    assert counters['sampled']==N*(max(times)+1) and counters['invalid']==0
    assert counters['complete']==sum(e['status']==0 for e in events) and counters['partial']==sum(e['status']==1 for e in events) and counters['suspended']==len(open_events)
    mass_error=max(abs(d['mass']-1) for d in snapshots.values());assert mass_error<=1e-12
    strength=np.array([e['first']['mach_pressure'] for e in events]);temp=np.array([e['first']['mach_temperature'] for e in events])
    summary={'particles':N,'completed':counters['complete'],'partial':counters['partial'],'open':len(open_events),'unique_event_identities':True,'contiguous_sequences':True,'all_particle_ids_preserved':True,
             'emitted_per_particle':len(events)/N,'recorded_starts_per_particle':len(all_events)/N,'particles_with_emitted_events':len(counts),'zero_emitted_particles':N-len(counts),
             'completed_count_histogram':dict(sorted(Counter(counts[p] for p in range(N)).items())),
             'flags':{name:sum(bool(e['flags']&bit) for e in all_events) for name,bit in FLAGS.items()},
             'mach_pressure_quantiles':np.quantile(strength[np.isfinite(strength)],[.1,.5,.9]).tolist(),
             'finite_pressure_estimates':int(np.isfinite(strength).sum()),'finite_temperature_estimates':int(np.isfinite(temp).sum()),
             'completed_duration_bounds':len(duration),'adjacent_complete_waiting_bounds':len(waiting),'max_mass_error':mass_error,'counters':counters,
             'field_evolution':[{k:d[k] for k in ('time','cycle','mass','mach_volume','mach_mass','supersonic_volume_fraction','supersonic_mass_fraction','kinetic','gas_energy')} for d in snapshots.values()]}
    return summary,{'events':events,'open':open_events,'all':all_events,'counts':counts,'by_particle':by_particle,'duration':np.array(duration),'waiting':np.array(waiting),'selected_ids':representative_ids,'selected_observations':selected,'final_positions':final_positions}


def detector_snapshot(binary,snapshot,work):
    n=snapshot['n'];f=snapshot['fields'];fixture=work/('detector-'+str(snapshot['cycle'])+'.txt');output=fixture.with_suffix('.csv')
    xyz=lambda a:a.transpose(2,1,0).ravel()
    values=np.column_stack([xyz(f[k]) for k in ('dens','velx','vely','velz','eint')])
    header=f'{n} {n} {n}\n1.4 1.1 1000 1e-30 0 4 0 0 {1/n:.17g} {1/n:.17g} {1/n:.17g}'
    np.savetxt(fixture,values,fmt='%.17g',header=header,comments='')
    with output.open('w') as stream:subprocess.run([str(binary),str(fixture)],stdout=stream,check=True,timeout=180)
    with output.open() as stream:header=next(line for line in stream if line.startswith('# i,j,k')).removeprefix('# ').strip().split(',')
    keys=('valid','candidate','inside','mach_pressure','mach_temperature')
    data=np.loadtxt(output,delimiter=',',comments='#',usecols=[header.index(k) for k in keys])
    result={k:data[:,i].reshape((n,n,n)).transpose(2,1,0) for i,k in enumerate(keys)}
    assert np.all(result['valid']==1)
    return result,{'generation':snapshot['cycle'],'fixture_sha256':sha(fixture),'native_csv_sha256':sha(output),'detector_executable_sha256':sha(binary),'valid_cells':n**3,'inside_cells':int(result['inside'].sum()),'width_cells':4,'boundary':'native fixture uses true periodic neighbors','snapshot_recomputation':True}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True);p.add_argument('--detector',type=Path,required=True);p.add_argument('--workdir',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--reuse',action='store_true');a=p.parse_args()
    binary=a.executable.resolve();work=a.workdir.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    if not a.reuse:launch(binary,work)
    execution=json.loads((work/'execution.json').read_text());assert execution['status']=='simulations_completed' and execution['executable_sha256']==sha(binary)
    assert {r['name'] for r in execution['runs']}==set(CASES)
    assert all(r['exit']==0 for r in execution['runs'])
    for row in execution['runs']:
        d=work/row['name'];assert row['input_sha256']==sha(d/'input.athinput') and row['log_sha256']==sha(d/'run.log')
    # Whole native thermo bytes and every saved fluid/force array must agree.
    thp=lambda name:next((work/name/'prtcl_thermo_history').glob('*.thp'))
    assert sha(thp('baseline'))==sha(thp('tracking_off'))
    assert sha(thp('coarse32'))==sha(thp('coarse32_matched'))
    for on,off in [('baseline','tracking_off'),('coarse32','coarse32_matched')]:
        paths={p.name:p for p in (work/on/'bin').glob('*.bin')}
        other={p.name:p for p in (work/off/'bin').glob('*.bin')};assert paths.keys()==other.keys()
        for name,path in paths.items():
            av,bv=read_binary(path),read_binary(other[name])
            assert av['time']==bv['time'] and av['cycle']==bv['cycle'] and av['mb_data'].keys()==bv['mb_data'].keys()
            assert all(np.array_equal(av['mb_data'][field],bv['mb_data'][field]) for field in av['mb_data'])
        del av,bv;gc.collect()
    assert not (work/'tracking_off'/'shock_events').exists()
    report={**provenance(binary),'scope':'Evolving, continuously forced adiabatic supersonic demonstration; not stationary turbulence or calibrated physical shock census.','scientific_signoff':'pending','feature_source_commit':'ad4225edea78d255b67df97fe8f593c0c268fa1e','execution':execution,
            'configuration':{'baseline_cells':64,'particles':N,'initial_rms_mach':5,'gamma':GAMMA,'dedt':50,'tlim':END,'forcing_seed':2718,'initial_velocity_seed':24680,'particle_seed':11,'width_cells':4,'min_mach':1.1,'tcorr':.1,'dt_update':.01,'npeak':2,'sol_fraction':.5},
            'non_interference':'Entire native thermo bytes and all saved double primitive/conserved/force fields exact, baseline tracking-on/off. Same coarse32 flow/particle bytes with widths4/2.',
            'limitations':['Native RK2 forcing thermodynamic oracle fails; driver unchanged.','No stationary ensemble, physical shock truth or individual material crossing times.','All episodes carry no_plateau and cadence_limited; merged spans can contain brief outside intervals.','CFL changes fluid evolution, MC random draws and observation cadence jointly.','32³ width 4 changes physical detector half-span; width 2 matches64³ width 4 half-span but still couples mesh/MC/fluid effects.','No independent random-seed ensemble, HPC/MPI/GPU or accepted physical convergence tolerance.'],
            'cases':{},'snapshot_detectors':[],'figures':[]}
    fs={};details={}
    for name in LABELS:
        fs[name]=fields(work/name);report['cases'][name],details[name]=analyze(work/name,fs[name]);assert abs(fs[name][0]['mach_volume']-5)<=1e-12
    baseline=fs['baseline'];selected=[min(baseline.values(),key=lambda d:abs(d['time']-t)) for t in (.05,.15,.3)]
    masks={}
    for snapshot in selected:
        masks[snapshot['cycle']],record=detector_snapshot(a.detector.resolve(),snapshot,work);report['snapshot_detectors'].append(record)
    np.savez_compressed(out/'flow-and-mask-data.npz',**{str(d['cycle'])+'_'+key:value for d in selected for key,value in {**d['fields'],**masks[d['cycle']]}.items()})
    with (out/'episodes.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['case','outcome','particle_id','sequence','flags','pre_time','first_inside','last_inside','post_time','pressure_mach','temperature_mach'])
        for name,dd in details.items():
            for outcome,group in [('emitted',dd['events']),('open',dd['open'])]:
                for e in group:writer.writerow([name,outcome,e['particle_id'],e['sequence'],e['flags'],e['pre']['time'] if e['pre']['valid'] else '',e['first']['time'],e['last']['time'],e['post']['time'] if e['post']['valid'] else '',e['first']['mach_pressure'],e['first']['mach_temperature']])
    with (out/'evolution.csv').open('w',newline='') as stream:
        keys=list(report['cases']['baseline']['field_evolution'][0]);writer=csv.writer(stream);writer.writerow(['case']+keys)
        for name,s in report['cases'].items():
            for row in s['field_evolution']:writer.writerow([name]+[row[k] for k in keys])
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    def save(fig,name,title,caption):
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=250);fig.savefig(out/(name+'.pdf'));fig.savefig(out/(name+'.svg'));plt.close(fig)
        report['figures'].append({'name':name,'title':title,'caption':caption})
    def cdf(ax,x,label,**kw):
        x=np.sort(np.asarray(x));assert len(x)>0;ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=label,**kw)
    fig,axes=plt.subplots(2,3,figsize=(8,5.3),sharex=True,sharey=True)
    density=[np.log10(d['rho'][32]) for d in selected];mach=[d['mach_local'][32] for d in selected]
    for col,d in enumerate(selected):
        for row,values in enumerate((density,mach)):
            im=axes[row,col].imshow(values[col],origin='lower',extent=(0,1,0,1),cmap='magma' if row==0 else 'viridis',vmin=min(x.min() for x in values),vmax=max(x.max() for x in values))
            mask=masks[d['cycle']]['inside'][32];axes[row,col].contour((np.arange(64)+.5)/64,(np.arange(64)+.5)/64,mask,levels=[.5],colors=['#00e5ff'],linewidths=.45)
            axes[row,col].set_title(f't={d["time"]:.4f}, generation {d["cycle"]}');fig.colorbar(im,ax=axes[row,col],shrink=.75,label='log10 density' if row==0 else 'Local |v| / sound speed')
        axes[1,col].set_xlabel('x / box length')
    axes[0,0].set_ylabel('y / box length');axes[1,0].set_ylabel('y / box length')
    save(fig,'01_flow','Forced supersonic flow and shock regions','Actual z=0.5078125 slices of the 64³ baseline. Colors show log10 density (top) and local speed/sound speed (bottom); each row uses shared color limits. Cyan outlines are the unreduced native detector region recomputed from the exact saved primitive state, with true periodic neighbors, width 4 and min_mach=1.1. The three generations are the available saved states nearest t=.05,.15,.30. Local flow Mach is distinct from detector shock-Mach estimates. Masks demonstrate implemented numerical membership, not a certified physical shock catalogue.')
    fig,ax=plt.subplots(1,2,figsize=(8,4));final=selected[-1];dd=details['baseline'];pos=dd['final_positions'];which=np.floor(pos[:,2]*64).astype(int)==32
    ax[0].imshow(np.log10(final['rho'][32]),origin='lower',extent=(0,1,0,1),cmap='Greys');im=ax[0].scatter(pos[which,0],pos[which,1],c=[dd['counts'][pid] for pid in np.flatnonzero(which)],s=12,cmap='plasma',vmin=0,vmax=max(dd['counts'].values()));fig.colorbar(im,ax=ax[0],label='Completed episodes per particle')
    ax[1].imshow(masks[final['cycle']]['inside'][32],origin='lower',extent=(0,1,0,1),cmap='Blues',vmin=0,vmax=1);ax[1].scatter(pos[which,0],pos[which,1],s=5,c='tomato')
    for x in ax:x.set_xlabel('x / box length');x.set_ylabel('y / box length')
    save(fig,'02_particles','Particle histories over the final flow',f'Only the {which.sum()} actual particles whose final parent lies in the plotted z-cell are shown, avoiding a projection of unrelated depths. Left: final density and per-particle accumulated completed-episode count. Right: recomputed native region and the same particles. Counts include their whole recorded history and are not instantaneous shock membership. All 8,192 particle IDs survive; cell-center MC samples are not smooth material trajectories.')
    fig,ax=plt.subplots(1,3,figsize=(8,3.1))
    for color,name in enumerate(('baseline','half_cfl','coarse32_matched')):
        evo=report['cases'][name]['field_evolution'];t=[x['time'] for x in evo];ax[0].plot(t,[x['mach_volume'] for x in evo],color=f'C{color}',label=LABELS[name]);ax[1].plot(t,[x['supersonic_volume_fraction'] for x in evo],color=f'C{color}');ax[2].plot(t,[x['kinetic'] for x in evo],color=f'C{color}');ax[2].plot(t,[x['gas_energy'] for x in evo],color=f'C{color}',ls='--')
    ax[0].axhline(1,color='k',ls=':');ax[0].set_ylabel('Volume rms local Mach');ax[1].set_ylabel('Supersonic volume fraction');ax[2].set_ylabel('Energy / box volume');ax[0].legend(fontsize=6)
    for x in ax:x.set_xlabel('Code time')
    save(fig,'03_evolution','Supersonic evolution and heating','Volume rms local Mach is sqrt(mean(|v|²/(gamma P/rho))); supersonic volume fraction counts cells with local Mach>1. Energy panel: solid kinetic, dashed gas internal energy, using the same case colors. Continuous edot=50 forcing is present, but adiabatic shock/source heating increases sound speed. All saved baseline rms Mach values exceed 1; this is an evolving forced realization, not a stationary turbulence ensemble. Total energy includes actual forcing work and the independently measured inherited RK2 source bias; edot is an input normalization target, not a validated physical injection rate.')
    fig,ax=plt.subplots(1,3,figsize=(8,3.2));bins=np.linspace(0,END,13)
    for name in ('baseline','half_cfl','coarse32_matched'):
        dd=details[name];hist=Counter(dd['counts'][pid] for pid in range(N));ax[0].plot(sorted(hist),[hist[k]/N for k in sorted(hist)],'-o',ms=3,label=LABELS[name]);start=np.array([e['first']['time'] for e in dd['all']]);cdf_times=np.sort(start);ax[1].step(cdf_times,np.arange(1,len(start)+1)/N,where='post');rate,_=np.histogram(start,bins);ax[2].stairs(rate/(N*np.diff(bins)),bins)
    ax[0].set_xlabel('Completed episodes / particle');ax[0].set_ylabel('Fraction of all 8,192 particles');ax[0].legend(fontsize=6);ax[1].set_xlabel('First-inside recording time');ax[1].set_ylabel('Cumulative recorded starts / particle');ax[2].set_xlabel('First-inside recording time');ax[2].set_ylabel('Recorded starts / particle / code time')
    save(fig,'04_counts','How often particles meet detected regions','The count distribution includes all 8,192 particles, including zero-encounter particles. Cumulative/rate panels count each canonical committed or terminal-open identity once at its first-inside observation; left-censored records, if present, are listed separately in the flags table. Rate bins have width 0.025 code time. These are sampled numerical episode recording rates, not exact physical shock-crossing rates; correlated particles and merged flicker prevent a naive independent-event error model. The coarse matched-span case uses width 2 at 32³.')
    fig,ax=plt.subplots(1,2,figsize=(8,3.6))
    for name in ('baseline','half_cfl','coarse32_matched'):
        e=details[name]['events'];v=np.array([x['first']['mach_pressure'] for x in e]);v=v[np.isfinite(v)];cdf(ax[0],v,LABELS[name])
    e=details['baseline']['events'];mp=np.array([x['first']['mach_pressure'] for x in e]);mt=np.array([x['first']['mach_temperature'] for x in e]);valid=np.isfinite(mp)&np.isfinite(mt);h=ax[1].hexbin(mp[valid],mt[valid],gridsize=35,bins='log',mincnt=1,cmap='viridis');fig.colorbar(h,ax=ax[1],label='Episodes per bin');lim=max(mp[valid].max(),mt[valid].max());ax[1].plot([1,lim],[1,lim],'k:',lw=1);ax[0].legend(fontsize=6);ax[0].set_xlabel('First-inside pressure Mach estimate');ax[0].set_ylabel('Completed-event cumulative fraction');ax[1].set_xlabel('Pressure Mach estimate');ax[1].set_ylabel('Temperature Mach estimate')
    save(fig,'05_strength','Recorded strength evidence',f'Left: finite first-inside pressure-Mach estimates for completed episodes, normalized per case. Right: baseline pressure/temperature estimators for {valid.sum()}jointly finite events; {len(e)-valid.sum()}events without both finite estimators remain in counts/CSV and validity flags. The diagonal is estimator agreement, not an independent physical truth test. Fixed endpoints at±width along the dominant normal axis have no certified plateaus. Either estimator can accept membership; a pressure estimate below 1.1 is therefore possible even with min_mach=1.1.')
    fig,ax=plt.subplots(1,2,figsize=(8,3.6));dur=details['baseline']['duration'];cdf(ax[0],dur[:,0],'Lower span bound');cdf(ax[0],dur[:,1],'Upper span bound');e=details['baseline']['events'];entry=[x['first']['time']-x['pre']['time'] for x in e if x['pre']['valid']];exitw=[x['post']['time']-x['last']['time'] for x in e if x['post']['valid']];cdf(ax[1],entry,'Entry bracket width');cdf(ax[1],exitw,'Exit bracket width')
    for x in ax:x.legend(fontsize=7);x.set_xlabel('Code time');x.set_ylabel('Cumulative fraction')
    save(fig,'06_brackets','Episode spans retain their timing brackets','For uncensored completed episodes, entry lies in[pre,first]and exit in[last,post]. The allowed elapsed-span bounds are[last-first,post-pre]; both endpoint CDFs are shown without inventing midpoints or exact crossing times. Right: widths of the underlying entry/exit observation brackets. Merged flicker can span brief outside intervals, so elapsed span is not continuous residence time. Cadence limitation and stochastic cell hopping remain explicit; these curves are not physical timing-error distributions.')
    fig,ax=plt.subplots(1,2,figsize=(8,3.6));wait=details['baseline']['waiting'];cdf(ax[0],wait[:,0],'Lower gap bound');cdf(ax[0],wait[:,1],'Upper gap bound');op=details['baseline']['open'];lower=[x['last']['time']-x['first']['time'] for x in op];cdf(ax[1],lower,'Open-episode observed inside-span lower bound');ax[0].set_xlabel('Between-episode mask gap (code time)');ax[1].set_xlabel('Observed open span (code time)')
    for x in ax:x.set_ylabel('Cumulative fraction');x.legend(fontsize=7)
    save(fig,'07_recurrence','Recurrence and right censoring',f'Left: {len(wait)}pairs of consecutive uncensored completed sequences for the same particle. Gap bounds are[max(0,next.pre-previous.post),max(0,next.first-previous.last)]. A zero lower bound means overlapping evidence brackets, not a known zero waiting time. Right: {len(op)}terminal-open identities retained separately, with observed inside-span lower bounds; no finite completed duration is assigned. No stationary waiting-time law, survival fit or independent-episode assumption is claimed.')
    fig,ax=plt.subplots(1,2,figsize=(8,4));names=list(FLAGS)
    for i,name in enumerate(('baseline','half_cfl','coarse32_matched')):
        flags=report['cases'][name]['flags'];den=len(details[name]['all']);ax[0].barh(np.arange(len(names))+i*.22,[flags[x]/den for x in names],height=.22,label=LABELS[name]);s=report['cases'][name];ax[1].bar(np.arange(3)+i*.22,[s['completed'],s['partial'],s['open']],width=.22)
    ax[0].set_yticks(np.arange(len(names))+.22,[x.replace('_',' ') for x in names]);ax[0].set_xlabel('Fraction with flag (overlap allowed)');ax[0].legend(fontsize=6);ax[1].set_xticks(np.arange(3)+.22,['Complete','Partial','Open']);ax[1].set_ylabel('Recorded identities')
    save(fig,'08_quality','Validity flags and censored outcomes','Every event flag is shown, including universal no_plateau/cadence_limited; overlapping flags do not form a probability partition. Direction ambiguity, temperature invalidity and merged flicker remain evidence, not filtered failures. The zero downstream_reentry event-flag bar does not imply zero suppression: the separate surviving-state counter records 75,455 suppressed-entry attempts in the baseline (84,032 at half CFL). These counters are not extra encounter identities. Open terminal views retain canonical identity and can resume after restart; they are not counted as a second completed event. Complete means the numerical tracker obtained its closure evidence, not human certification of a physical passage.')
    fig,ax=plt.subplots(2,1,figsize=(8,4.8),sharex=True);dd=details['baseline'];chosen=dd['selected_ids']
    for row,pid in enumerate(chosen):
        for e in dd['by_particle'][pid]:
            first,last=e['first']['time'],e['last']['time'];ax[0].plot([first,last],[row,row],lw=5,color='C0');ax[0].plot([e['pre']['time'],first],[row,row],lw=2,color='C1');
            if e['post']['valid']:ax[0].plot([last,e['post']['time']],[row,row],lw=2,color='C1')
            if e['flags']&64:ax[0].plot(last,row,'>',color='C3',ms=7)
        oo=sorted(dd['selected_observations'][pid],key=lambda o:o['time']);ax[1].plot([o['time'] for o in oo],[o['fields']['density'] for o in oo],label=f'ID{pid}')
    ax[0].set_yticks(range(len(chosen)),[f'ID{pid}: {dd["counts"][pid]}complete' for pid in chosen]);ax[1].set_yscale('log');ax[1].set_ylabel('Native particle density');ax[1].set_xlabel('Code time');ax[1].legend(ncol=3,fontsize=7)
    save(fig,'09_histories','Individual particle records without invented trajectories','Particles are selected deterministically at count ranks 0, 25, 50, 75, 90 and 100 percent of the full population, breaking ties by ID; their IDs/counts are printed. Blue segments are first-to-last inside evidence; orange marks entry/exit timing brackets. Red triangles identify open terminal episodes. The lower panel shows actual sampled native density histories; connecting samples is a visual guide, not a continuous material path. Count selection is illustrative rather than an unbiased ensemble sample.')
    fig,ax=plt.subplots(1,3,figsize=(8,3.6));names=list(LABELS);counts=[report['cases'][n]['completed']/N for n in names];opens=[report['cases'][n]['open']/N for n in names];ax[0].bar(range(4),counts,color=['C0','C1','C2','C3']);ax[1].bar(range(4),opens,color=['C0','C1','C2','C3']);ax[2].bar(range(4),[report['cases'][n]['field_evolution'][-1]['mach_volume'] for n in names],color=['C0','C1','C2','C3'])
    short=['64³','64³ CFL/2','32³ w4','32³ w2']
    for x in ax:x.set_xticks(range(4),short,rotation=30,ha='right')
    ax[0].set_ylabel('Complete / particle');ax[1].set_ylabel('Open / particle');ax[2].set_ylabel('Final volume rms Mach')
    save(fig,'10_sensitivity','Timestep, grid and detector-support sensitivity','64³ CFL/2 changes fluid integration, MC draws and observation cadence together.32³ width 4 also doubles the detector physical half-span;32³ width 2 matches the 64³ width 4 half-span(1/16 of box length). The two 32³ runs have byte-identical native particle histories and exact saved fluid/forcing arrays, isolating their detector-width difference. Grid comparisons still mix fluid resolution, stochastic transport and detection. These are measured sensitivities, not convergence acceptance or independent uncertainty estimates. No tolerances were tuned to reduce differences.')
    assert len(report['figures'])==10
    report['status']='showcase_completed_scientific_review_pending';report['source_sha256']={name:sha(ROOT/name) for name in ('tst/shocks-mc/supersonic_showcase.py','inputs/tests/shock_mc_showcase.athinput','src/shocks/detector.hpp','src/shocks/encounter.hpp','src/particles/shock_particles.cpp','src/srcterms/turb_driver.cpp','tst/shocks/detector_test.cpp')}
    for name in CASES:
        dest=out/'runs'/name;dest.mkdir(parents=True)
        for file in ('input.athinput','run.log'):shutil.copy2(work/name/file,dest/file)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    text='<h1>Forced supersonic turbulence: particle shock encounters</h1><p>64³ periodic cells,8,192 MC particles, initial Mach 5, edot50, t=.3. Evolving adiabatic flow; scientific calibration pending.</p>'
    for f in report['figures']:text+='<article><h2>'+html.escape(f['title'])+'</h2><img src="'+f['name']+'.png"><p>'+html.escape(f['caption'])+'</p><p><a href="'+f['name']+'.pdf">VectorPDF</a> · <a href="'+f['name']+'.svg">SVG</a></p></article>'
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Supersonic shock showcase</title><style>body{font:17px/1.5 system-ui;max-width:1200px;margin:2em auto;padding:0 1em}img{width:100%}article{margin:3em 0}</style>'+text)
    print(json.dumps({'status':report['status'],'cases':{name:{k:s[k] for k in ('completed','partial','open','emitted_per_particle')} for name,s in report['cases'].items()},'output':str(out)}))


if __name__=='__main__':main()
