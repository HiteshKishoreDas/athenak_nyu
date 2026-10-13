"""Native-forcing/MC shock compatibility; independent shear heat oracle."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile

from common import ROOT, edit_input, provenance, sha, thermo
from integration import traces, canonical
sys.path.insert(0,str(ROOT/'vis/python'))
from bin_convert import read_binary
from shock_history import read_history, read_suspended
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def volumes(directory):
    result={}
    for p in sorted((directory/'bin').glob('*.bin')):
        d=read_binary(p)
        assert d['n_mbs']==1
        fields={k:v[0] for k,v in d['mb_data'].items()}
        assert all(v.dtype==np.float64 and np.isfinite(v).all() for v in fields.values())
        kind='force' if 'force1' in fields else 'fluid'
        key=kind,d['cycle']
        if key in result:
            previous=result[key]
            assert previous['time']==d['time']
            assert previous['fields'].keys()==fields.keys()
            assert all(np.array_equal(previous['fields'][k],fields[k]) for k in fields)
            continue  # Native final output can repeat an identical generation.
        result[key]={'path':p,'time':d['time'],'cycle':d['cycle'],'fields':fields}
    return result


def stats(directory):
    events=read_history(directory/'shock_events');suspended=read_suspended(directory/'shock_events')
    all_events=events+suspended
    assert len({e['identity'] for e in all_events})==len(all_events)
    return {'completed':sum(e['status']==0 for e in events),
            'partial':sum(e['status']==1 for e in events),'suspended':len(suspended),
            'unique_emitted_particles':len({e['particle_id'] for e in events}),
            'merged_flicker':sum(bool(e['flags']&4) for e in all_events),
            'direction_ambiguous':sum(bool(e['flags']&256) for e in all_events),
            'strengths': [e['first']['mach_pressure'] for e in events if np.isfinite(e['first']['mach_pressure'])]}


def evolution(v):
    result=[]
    for (kind,g),d in sorted(v.items()):
        if kind!='fluid':continue
        f=d['fields'];rho=f['dens']
        ke=sum(f[k]**2 for k in ('mom1','mom2','mom3'))/(2*rho)
        p=(1.4-1)*(f['ener']-ke)
        assert np.all(rho>0) and np.all(p>0)
        result.append({'time':d['time'],'cycle':g,'mass':float(rho.mean()),
                       'total_energy':float(f['ener'].mean()),'gas_energy':float((f['ener']-ke).mean()),
                       'mach':float(np.sqrt(np.mean(2*ke/(1.4*p))))})
    assert result
    assert max(abs(d['mass']-result[0]['mass']) for d in result)<=1e-12
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();binary=a.executable.resolve();work=a.workdir.resolve();out=a.output.resolve()
    work.mkdir(parents=True,exist_ok=False);out.mkdir(parents=True,exist_ok=False)
    report={**provenance(binary),'status':'in_progress','runs':[],'checks':[],
            'scientific_signoff':'pending','scope':'Serial/double RK2 uniform Newtonian ideal hydro; native accel_rms/edot forcing without opt-in.'}
    base=(ROOT/'inputs/tests/shock_mc_forced.athinput').read_text()
    base=edit_input(base,{'output4/file_type':'rst','output4/dt':-1,
                         'output5/file_type':'bin','output5/variable':'turb_force',
                         'output5/id':'force','output5/dt':.025,'output5/data_precision':'real'})
    def save(): (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    def run(name,changes=None,restart=None,overrides=(),expected=None,remove=()):
        directory=work/name;directory.mkdir()
        inp=directory/'input.athinput'
        if restart is None:
            resolved=dict(changes or {})
            resolved['shock_tracking/directory']=str(directory/'shock_events')
            s=edit_input(base,resolved)
            s='\n'.join(l for l in s.splitlines() if not any(l.startswith(k+' =') for k in remove))+'\n'
            inp.write_text(s);command=[str(binary),'-i',str(inp)]
        else:command=[str(binary),'-r',str(restart)]+list(overrides)
        start=time.monotonic()
        with (directory/'run.log').open('w') as log:
            try:code=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
            except subprocess.TimeoutExpired:code='timeout'
        log=(directory/'run.log').read_text()
        row={'name':name,'command':command,'exit':code,'seconds':time.monotonic()-start,
             'expected_error':expected,'log_sha256':sha(directory/'run.log')}
        if inp.exists():row['input_sha256']=sha(inp)
        report['runs'].append(row);save()
        assert code==0 if expected is None else isinstance(code,int) and code!=0 and expected in log,row
        return directory
    shear={'mesh/nx1':16,'mesh/nx2':16,'mesh/nx3':16,'meshblock/nx1':16,'meshblock/nx2':16,'meshblock/nx3':16,
           'tracer_seed_initial/count_per_event':128,'initial_perturbations/velocity_rms':0,
           'turb_driving/accel_rms':1,'turb_driving/spectrum':'power_law','turb_driving/nlow':1,'turb_driving/nhigh':1,
           'turb_driving/min_kx':1,'turb_driving/max_kx':1,'turb_driving/min_ky':0,'turb_driving/max_ky':0,
           'turb_driving/min_kz':0,'turb_driving/max_kz':0,'turb_driving/sol_fraction':1,
           'shock_tracking/lineage':'forcing-shear-control','time/tlim':.05,
           'output3/dcycle':1,'output3/dt':-1,'output5/dcycle':1,'output5/dt':-1}
    step=run('shear_one_step',dict(shear,**{'time/nlim':2}))
    v=volumes(step);initial=v['fluid',1];final=v['fluid',2];f=final['fields'];dt=final['time']-initial['time']
    accel=v['force',2]['fields'];a2=sum(accel['force'+str(i)]**2 for i in (1,2,3))
    assert np.mean(a2)>0  # Native OU starts with zero force until its first update.
    ke=sum(f['mom'+str(i)]**2 for i in (1,2,3))/(2*f['dens'])
    gas=f['ener']-ke;gas0=initial['fields']['ener']
    assert np.all(initial['fields']['dens']==1) and np.array_equal(f['dens'],initial['fields']['dens'])
    assert all(np.all(initial['fields']['mom'+str(i)]==0) for i in (1,2,3))
    assert np.max(abs(accel['force1']))<=1e-12  # Transverse shear, projection roundoff only.
    expected_bias=3/8*a2*dt**2
    error=float(np.max(abs(gas-gas0-expected_bias)))
    mom_error=max(float(np.max(abs(f['mom'+str(i)]-accel['force'+str(i)]*dt))) for i in (1,2,3))
    assert error<=1e-12 and mom_error<=1e-12
    measured_heat=float(np.mean(gas-gas0));assert measured_heat>1e-12
    report['inherited_forcing_heat']={'physical_oracle':'Divergence-free single-mode transverse acceleration preserves rho and gas internal energy in the continuum.',
        'physical_gas_energy_oracle_passed':False,'measured_mean_excess_gas_energy':measured_heat,
        'predicted_mean_excess_gas_energy':float(expected_bias.mean()),'max_discrete_prediction_error':error,
        'max_momentum_error':mom_error,'dt':dt,'formula':'3 rho |a|^2 dt^2 / 8 for native RK2 at zero flux divergence',
        'driver_modified':False}
    control=run('smooth_shear',shear);assert stats(control)['completed']==0 and stats(control)['partial']==0 and stats(control)['suspended']==0
    report['checks'].append({'ids':['EXT06','PHY09'],'smooth_forcing_events':0,'note':'Negative detection control, not thermodynamic neutrality.'})
    timing=run('timing',dict(shear,**{'initial_perturbations/velocity_rms':3,
        'turb_driving/accel_rms':10,'time/nlim':6,'shock_tracking/trace':'true',
        'output2/dcycle':1,'output3/variable':'hydro_w','shock_tracking/lineage':'forced-timing'}))
    obs=traces(timing/'shock_events');native,times=thermo(timing);vv=volumes(timing)
    assert obs.keys()==native.keys()
    max_p_error=0
    for (g,pid),o in obs.items():
        n=native[g,pid];f=vv['fluid',g]['fields']
        ijk=[min(15,int(x*16)) for x in o['x']];i,j,k=ijk
        assert o['time']==times[g] and tuple(o['x'])==n['x']
        # Binary metadata uses six significant digits; thermo retains full doubles.
        assert float(format(o['time'],'.6g'))==vv['fluid',g]['time']
        assert o['rho']==n['fields']['density']==f['dens'][k,j,i]
        expected=(1.4-1)*f['eint'][k,j,i]
        max_p_error=max(max_p_error,abs(o['pressure']-expected))
        assert abs(o['pressure']-expected)<=1e-12*max(1,abs(expected))
    assert max(abs(vv['force',6]['fields'][key]).max() for key in vv['force',6]['fields'])>0
    report['checks'].append({'ids':['INT06','INT08','EXT06'],'forced_final_step_samples':len(obs),'one_per_generation_particle':True,
                             'time_position_density_exact':True,'binary_time_metadata':'exact six-significant-digit serialization match','pressure_relative_absolute_tolerance':1e-12,'max_pressure_absolute_error':float(max_p_error)})
    def assert_on_off(on,off,physical_outflow=False):
        onv,offv=volumes(on),volumes(off);assert onv.keys()==offv.keys()
        for key in onv:
            x,y=onv[key],offv[key];assert x['time']==y['time']
            assert x['fields'].keys()==y['fields'].keys()
            assert all(np.array_equal(x['fields'][k],y['fields'][k]) for k in x['fields'])
        tracked,times=thermo(on);native,native_times=thermo(off)
        assert times==native_times
        lost=set()
        if physical_outflow:
            # Existing enabled outflow deletes; the legacy disabled path wraps.
            # Derive the expected exits from source position and realized cell hop.
            transfers=sorted((on/'shock_events').glob('transfer*.rank0.csv'),
                             key=lambda p:int(p.name.split('.')[0].removeprefix('transfer')))
            assert {int(p.name.split('.')[0].removeprefix('transfer')) for p in transfers}==set(times)-{0}
            for path in transfers:
                generation=int(path.name.split('.')[0].removeprefix('transfer'))
                lines=path.read_text().splitlines()
                rows=[list(map(float,row)) for row in csv.reader(lines[1:-1])]
                assert lines[-1]=='# END '+str(len(rows))
                assert all(len(row)==17 for row in rows)
                assert len({int(row[0]) for row in rows})==len(rows)
                assert {int(row[0]) for row in rows}=={pid for g,pid in tracked if g==generation-1}
                for row in rows:
                    assert row[16].is_integer() and 0<=row[16]<=6
                    pid,direction=int(row[0]),int(row[16])
                    assert tuple(row[12:15])==tracked[generation-1,pid]['x']
                    x=row[12]+(-1/16 if direction==1 else 1/16 if direction==2 else 0)
                    if x<0 or x>=1:lost.add(pid)
                assert {pid for g,pid in tracked if g==generation}=={pid for g,pid in native if g==generation}-lost
            assert lost  # The fixture must exercise a real exit.
            assert native.keys()>=tracked.keys()
            assert all(row==native[key] for key,row in tracked.items())
        else:assert tracked==native
        assert not (off/'shock_events').exists()
        return sorted(lost)

    def assert_replay(on,checkpoint,resumed,generation):
        assert canonical(read_history(checkpoint/'shock_events'))==canonical(read_history(on/'shock_events'))
        restored,_=thermo(resumed);reference,_=thermo(on)
        assert restored=={key:value for key,value in reference.items() if key[0]>generation}
        assert canonical(read_suspended(checkpoint/'shock_events'))==canonical(read_suspended(on/'shock_events'))
        for kind in ('fluid','force'):
            end=max((d for (k,g),d in volumes(on).items() if k==kind),key=lambda d:d['cycle'])
            final=max((d for (k,g),d in volumes(resumed).items() if k==kind),key=lambda d:d['cycle'])
            assert end['cycle']==final['cycle'] and end['time']==final['time']
            assert end['fields'].keys()==final['fields'].keys()
            assert all(np.array_equal(end['fields'][key],final['fields'][key]) for key in end['fields'])
        return len(restored)

    # The ordinary input has no experimental_forcing parameter at all.
    assert 'experimental_forcing' not in base
    on=run('driven_on');off=run('driven_off',{'shock_tracking/enabled':'false'})
    assert_on_off(on,off)
    report['checks'].append({'ids':['CFG01','SYS02','EXT06'],'normalization':'accel_rms',
        'no_opt_in_required':True,'on_off_all_saved_fluid_force_particle_records_exact':True})
    half=run('half_cfl',{'time/cfl_number':.15,'shock_tracking/lineage':'forced-cfl015'})
    edot={'turb_driving/normalization':'edot','turb_driving/dedt':10,'shock_tracking/lineage':'forced-edot10'}
    edon=run('edot_on',edot,remove=('accel_rms',))
    edoff=run('edot_off',dict(edot,**{'shock_tracking/enabled':'false'}),remove=('accel_rms',))
    assert_on_off(edon,edoff)
    report['checks'].append({'ids':['CFG01','SYS02','EXT06'],'normalization':'edot',
        'no_opt_in_required':True,'on_off_all_saved_fluid_force_particle_records_exact':True})
    report['cases']={name:{**stats(path),'evolution':evolution(volumes(path))}
        for name,path in [('driven_on',on),('half_cfl',half),('edot_on',edon)]}
    for name,path in [('driven_on',on),('half_cfl',half),('edot_on',edon)]:
        data,times=thermo(path)
        assert all({pid for g,pid in data if g==cycle}==set(range(4096)) for cycle in times)
        report['cases'][name]['all_particle_ids_preserved']=True
    for name,reference,config,remove in [('accel_rms',on,{},()),('edot',edon,edot,('accel_rms',))]:
        checkpoint=run(name+'_checkpoint',dict(config,**{'time/nlim':30,'output4/dt':1}),remove=remove)
        cp=sorted((checkpoint/'rst').glob('*.rst'))[-1]
        prefix=read_history(checkpoint/'shock_events')
        resumed=run(name+'_resumed',restart=cp,overrides=('time/nlim=-1','output4/dt=-1','job/basename='+name+'_resumed'))
        rows=assert_replay(reference,checkpoint,resumed,30)
        report['checks'].append({'ids':['RST01','RST07','EXT06'],'normalization':name,
            'checkpoint_generation':30,'prefix_events':len(prefix),'replayed_particle_rows_exact':rows,
            'all_event_content_exact':True,'all_final_suspended_content_exact':True,
            'all_expected_resumed_rows_present':True,'final_fluid_and_ou_force_exact':True})
    checkpoint=work/'accel_rms_checkpoint';cp=sorted((checkpoint/'rst').glob('*.rst'))[-1]
    for name,override,error in [('changed_seed','turb_driving/rseed=2719','contract changed'),
                                ('changed_acceleration','turb_driving/accel_rms=11','configuration differs')]:
        run(name,restart=cp,overrides=(override,'output4/dt=-1'),expected=error)
    assert canonical(read_history(checkpoint/'shock_events'))==canonical(read_history(on/'shock_events'))
    small={'mesh/nx1':16,'mesh/nx2':16,'mesh/nx3':16,
           'meshblock/nx1':16,'meshblock/nx2':16,'meshblock/nx3':16,
           'tracer_seed_initial/count_per_event':128,'time/tlim':.04,
           'output2/dcycle':1,'output3/dcycle':1,'output3/dt':-1,
           'output5/dcycle':1,'output5/dt':-1}
    variants={'tiled':{'turb_driving/tile_nx':2},
              'localized':{'turb_driving/localization':'include','turb_driving/sigma_x1':.2},
              'delayed':{'turb_driving/tdriv_start':.01},
              'finite_duration':{'turb_driving/turb_flag':1,'turb_driving/tdriv_duration':.02},
              'outflow':{'mesh/ix1_bc':'outflow','mesh/ox1_bc':'outflow','shock_tracking/trace':'true'},
              'planar_2d':{'mesh/nx3':1,'meshblock/nx3':1,'turb_driving/driving_type':0,
                           'turb_driving/min_kz':0,'turb_driving/max_kz':0},
              'native_planar':{'turb_driving/driving_type':1,
                               'turb_driving/min_kz':1,'turb_driving/max_kz':1}}
    for name,config in variants.items():
        changes=dict(small,**config)
        variant_on=run(name+'_on',changes)
        variant_off=run(name+'_off',dict(changes,**{'shock_tracking/enabled':'false'}))
        lost=assert_on_off(variant_on,variant_off,physical_outflow=name=='outflow')
        force=volumes(variant_on)
        assert any(np.any(d['fields'][key]!=0) for (kind,g),d in force.items() if kind=='force' for key in d['fields'])
        observations,times=thermo(variant_on)
        preserved=all({pid for g,pid in observations if g==cycle}==set(range(128)) for cycle in times)
        if name!='outflow':assert preserved
        inactive=0
        if name in ('delayed','finite_duration'):
            for (kind,g),d in force.items():
                if kind!='force' or g==0:continue
                start=times[g-1]  # Saved force at t_(n+1) was used during [t_n,t_(n+1)].
                if (name=='delayed' and start<.01) or (name=='finite_duration' and start>=.02):
                    assert all(np.all(value==0) for value in d['fields'].values())
                    inactive+=1
            assert inactive>0
        report['checks'].append({'ids':['CFG01','SYS02','EXT06'],'variant':name,
            'on_off_all_saved_fluid_force_exact':True,'on_off_surviving_particle_records_exact':True,
            'nonzero_force_observed':True,'inactive_step_force_zero_checks':inactive,'physically_exited_ids':lost,
            'all_particle_ids_preserved':preserved,'outflow_deletion_allowed':name=='outflow',**stats(variant_on)})
    run('reject_nonfinite_tcorr',{'turb_driving/tcorr':'nan'},expected='native forcing parameters must be finite')
    run('reject_nonfinite_dedt',dict(edot,**{'turb_driving/dedt':'nan'}),remove=('accel_rms',),
        expected='native forcing parameters must be finite')
    report['count_change_fraction']=report['cases']['half_cfl']['completed']/report['cases']['driven_on']['completed']-1
    report['status']='passed_forcing_engineering_compatibility_with_inherited_forcing_heat_bias'
    report['limitations']=['Source thermodynamics have the independently measured native RK2 heating bias; no scientific strength approval.',
                          'One short32^3 realization; no stationary ensemble or convergence claim.',
                          'CFL jointly changes fluid, RNG draws and sampling; not a pure cadence test.',
                          'Both native normalizations allowed; short variant checks are non-interference evidence, not full physical qualification.',
                          'Other sources/MPI/GPU/AMR and physical shock truth remain unqualified.']
    fig,ax=plt.subplots(1,3,figsize=(14,4))
    for name in report['cases']:
        e=report['cases'][name]['evolution'];t=[x['time'] for x in e]
        ax[0].plot(t,[x['mach'] for x in e],'-o',label=name)
        ax[1].plot(t,[x['total_energy'] for x in e],label=name)
        strength=sorted(report['cases'][name]['strengths']);ax[2].step(strength,np.arange(1,len(strength)+1)/len(strength),label=name)
    ax[0].set_ylabel('Volume rms local Mach');ax[1].set_ylabel('Total energy / box volume');ax[2].set_ylabel('Event fraction');ax[2].set_xlabel('First-inside pressure Mach estimate')
    for x in ax[:2]:x.set_xlabel('Code time')
    for x in ax:x.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'driven_checks.png',dpi=150);plt.close(fig)
    for row in report['runs']:
        directory=work/row['name'];target=out/'runs'/row['name'];target.mkdir(parents=True)
        for name in ('run.log','input.athinput'):
            if (directory/name).exists():shutil.copy2(directory/name,target/name)
    report['source_sha256']={name:sha(ROOT/name) for name in ['CMakeLists.txt','src/particles/shock_particles.cpp','src/srcterms/turb_driver.cpp','src/srcterms/turb_driver.hpp','src/utils/random.hpp','src/hydro/hydro_tasks.cpp','src/hydro/hydro_update.cpp','tst/shocks-mc/forcing_checks.py','tst/shocks-mc/common.py','inputs/tests/shock_mc_forced.athinput']}
    save()
    with zipfile.ZipFile(out/'events-and-checkpoint.zip','x',zipfile.ZIP_DEFLATED) as z:
        for name in ('driven_on','half_cfl','edot_on','accel_rms_checkpoint','edot_checkpoint'):
            for p in (work/name/'shock_events').glob('*'):z.write(p,str(p.relative_to(work)))
        for name in ('accel_rms_checkpoint','edot_checkpoint'):
            for p in (work/name/'rst').glob('*'):z.write(p,str(p.relative_to(work)))
    (out/'sha256.json').write_text(json.dumps({str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},indent=2)+'\n')
    print(json.dumps({'status':report['status'],'counts':{k:{n:v[n] for n in ('completed','partial','suspended')} for k,v in report['cases'].items()},'native_heat_bias':measured_heat,'count_change_fraction':report['count_change_fraction']}))


if __name__=='__main__':
    main()
