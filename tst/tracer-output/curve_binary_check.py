"""Unified .crv native parity/oracles/restart/reader checks; retain every run/failure.

Use vis/.venv Python. Optional reference-runs is the previous SHOUT gate directory.
No initial-condition parameters or physical tolerances are fitted here.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import edit_input,provenance,sha
from bin_convert import read_binary
import curve_binary
import material_curves
from tracer_history import read_run
from tracer_thermo_history import read_curves,read_particles
from check import normalized

BASE='rho,p,t,s,v1,v2,v3'
BUDGET=2e-12  # unchanged CT independent physical-oracle budget


def require(condition,message):
    if not condition:raise ValueError(message)


def core(curves):
    return [dict(c,nodes=[{k:v for k,v in n.items() if k not in ('fields','gid','level')}
                         for n in c['nodes']]) for c in curves]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--prescribed-executable',type=Path,required=True)
    parser.add_argument('--reference-runs',type=Path)
    parser.add_argument('--workdir',type=Path,required=True)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    reference=args.reference_runs.resolve() if args.reference_runs else None
    report=dict(checks=[],runs=[],budget=BUDGET,runner_sha256=sha(__file__),
                scope='local Serial/double; not MPI/GPU/scientific approval')
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save();print(row,flush=True)
    def run(name,fixture='planar',changes=None,restart=None,overrides=(),reject=None):
        directory=work/name;directory.mkdir()
        source=(reference/fixture/'input.athinput' if reference and fixture in ('planar','amr_shock','forced') else
                ROOT/'inputs/tests/shocks-curves'/('amr.athinput' if fixture=='amr_shock' else fixture+'.athinput'))
        binary=(args.executable if fixture in ('planar','forced') else args.prescribed_executable).resolve()
        options={'output2/file_type':'prtcl_thermo_history','output2/source':'both',
                 'output2/variables':BASE+',shock','output2/dcycle':10 if fixture=='forced' else 1,
                 'output2/id':'samples','output4/dcycle':0,'output4/dt':-1}
        if fixture=='smooth':options.update({'output2/file_type':'crv','shock_tracking/enabled':'false',
            'output2/variables':'default,vmag,cs,mach','curve1/id':255})
        options.update(changes or {})
        options={k:v for k,v in options.items() if v is not None}
        path=directory/'input.athinput';path.write_text(edit_input(source.read_text(),options))
        command=[str(binary),'-r' if restart else '-i',str(restart or path),*overrides]
        with (directory/'run.log').open('w') as log:
            code=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
        report['runs'].append(dict(name=name,command=command,exit=code,input_used=not bool(restart),
            input_sha256=sha(path),restart_sha256=sha(restart) if restart else None,**provenance(binary)));save()
        require((code!=0 and reject in (directory/'run.log').read_text()) if reject else code==0,'native outcome: '+name)
        return directory
    cases={}
    def parity(fixture):
        directory=run(fixture,fixture);cases[fixture]=directory;current=read_run(directory)
        require(current['geometry'] and current['curve_nodes'],'unified streams missing')
        require(not list(directory.rglob('*.curves.*.txt')) and not (directory/'curve_thermo_history').exists(),
                'unified mode emitted automatic split geometry/thermo')
        states={s['cycle']:s for s in current['geometry']}
        rows,_=read_curves(directory/'curves')
        require(normalized(rows)==normalized(current['curve_nodes']),'direct CRV thermo adapter differs')
        count=0
        for state in states.values():
            require(all(key not in state for key in ('curve_events','encounters')),'complete events unexpectedly embedded')
            for curve in state['curves']:
                for node in curve['nodes']:
                    row=rows[state['cycle'],node['id']]
                    require(row['x']==tuple(node['x']) and normalized(row['fields'])==normalized(node['fields']),'geometry/thermo alignment differs')
                    require((row['gid'],row['level'])==(node['gid'],node['level']),'mesh location differs')
                    count+=1
        if reference:
            previous=read_run(reference/fixture)
            for key in ('curve_nodes','curve_times','particles','particle_times','curve_events','curve_suspended','particle_events','particle_suspended'):
                require(normalized(current[key])==normalized(previous[key]),'pre-change '+key+' differs')
            oldframes={s['cycle']:s for s in previous['geometry']}
            for cycle,state in states.items():
                require(normalized(core(state['curves']))==normalized(oldframes[cycle]['curves']),
                        'pre-change geometry/diagnostics differ')
            old=next((reference/fixture/'prtcl_thermo_history').glob('*.thp'))
            new=next((directory/'prtcl_thermo_history').glob('*.thp'))
            require(old.read_bytes()==new.read_bytes(),'particle THP bytes changed')
            for folder,pattern in [('bin','*.bin'),('tab','*.tab')]:
                for p in (directory/folder).glob(pattern):
                    old=reference/fixture/folder/p.name
                    if folder=='bin':
                        a,b=read_binary(p),read_binary(old)
                        require((a['cycle'],a['time'])==(b['cycle'],b['time']),'native BIN clocks changed')
                        require(a['mb_data'].keys()==b['mb_data'].keys(),'native BIN field schema changed')
                        for key,values in a['mb_data'].items():
                            require(np.array_equal(values,b['mb_data'][key],equal_nan=True),
                                    'native fluid/mesh arrays changed')
                    else:
                        require(p.read_bytes()==old.read_bytes(),'native fluid text bytes changed')
        return dict(frames=len(states),node_samples=count,reference_parity=bool(reference))
    check('CRV01_planar_exact_previous_geometry_thermo_shocks_particle_bytes',lambda:parity('planar'))
    check('CRV02_live_AMR_negative_ID_insertion_exact_parity',lambda:parity('amr_shock'))
    def explicit():
        directory=run('explicit_crv',changes={'output2/source':'particles','output4/file_type':'crv',
            'output4/id':'samples','output4/dcycle':1,'output4/variables':BASE+',shock'})
        a,b=read_run(cases['planar']),read_run(directory)
        for key in ('curve_nodes','geometry','curve_events','curve_suspended'):
            left=a[key];right=b[key]
            if key=='geometry':left=[s['curves'] for s in left];right=[s['curves'] for s in right]
            require(normalized(left)==normalized(right),'explicit CRV differs: '+key)
        return dict(explicit_and_source_both_equal=True)
    check('CRV03_explicit_native_output_without_mesh_variable',explicit)
    def oracle():
        directory=run('smooth','smooth',{'problem/flow':'acceleration','problem/a1':.3,'problem/a2':.2})
        cases['smooth']=directory;data=read_run(directory);maximum=0
        require(all(node<0 for _,node in data['curve_nodes']),'negative IDs not exercised')
        for row in data['curve_nodes'].values():
            t=row['time'];v=[.2+.3*t,-.1+.2*t,0];p=(4-max(map(abs,v)))**2/(5/3)
            expected=dict(density=1,pressure=p,temperature=p,entropy=math.log(p),
                internal_energy=p/(2/3),v1=v[0],v2=v[1],v3=0,
                vmag=math.hypot(*v),sound_speed=math.sqrt(5*p/3),mach=math.hypot(*v)/math.sqrt(5*p/3))
            maximum=max(maximum,max(abs(row['fields'][k]-v) for k,v in expected.items()))
        require(maximum<=BUDGET,'independent prescribed thermo oracle failed')
        return dict(samples=len(data['curve_nodes']),max_error=maximum,budget=BUDGET)
    check('CRV04_independent_accelerating_thermo_oracle',oracle)
    def diagnostics():
        maximum=0;births=0;levels=set()
        for s in read_run(cases['amr_shock'])['geometry']:
            for c in s['curves']:
                g=material_curves.geometry(c,s['domain'])
                maximum=max(maximum,abs(g['length']-c['length']))
                for i,n in enumerate(c['nodes']):
                    for field in ('edge_length','stretch','curvature'):
                        a,b=n[field],g[field][i]
                        require((math.isnan(a) and math.isnan(b)) or math.isfinite(a) and math.isfinite(b),'diagnostic NaN mask differs')
                        if math.isfinite(a):maximum=max(maximum,abs(a-b))
                    require(n['curvature_mask']==g['curvature_mask'][i],'curvature mask differs')
                    levels.add(n['level']);births+=n['birth']>0
        require(maximum<=BUDGET and births>0 and len(levels)>1,'AMR/geometry fixture ineffective or oracle failed')
        return dict(max_error=maximum,post_birth_samples=births,levels=sorted(levels),budget=BUDGET)
    check('CRV05_independent_geometry_and_actual_AMR_births',diagnostics)
    def periodic():
        # Reuse accepted production.py translation_wrap physical parameters.
        directory=run('periodic','smooth',{'problem/flow':'translation','problem/v1':2,
            'problem/v2':-1,'time/tlim':2.2})
        cases['periodic']=directory
        frames=read_run(directory)['geometry']
        require(any(n['image'][0] for s in frames for c in s['curves'] for n in c['nodes']),'periodic images not exercised')
        return dict(frames=len(frames),images_preserved=True,links_audited=True)
    check('CRV06_accepted_periodic_translation_lifts_and_links',periodic)
    def circle():
        # Accepted production.py circle_32 fixture and its existing tolerances.
        directory=run('circle','smooth',{'curve1/type':'circle','curve1/count':32,
            'curve1/x1':0,'curve1/x2':0,'curve1/radius':.7,'problem/v1':0,'problem/v2':0})
        cases['circle']=directory;frames=read_run(directory)['geometry']
        length=32*1.4*math.sin(math.pi/32);curvature=1/.7
        le=max(abs(c['length']-length) for s in frames for c in s['curves'])
        ce=max(abs(n['curvature']-curvature) for s in frames for c in s['curves'] for n in c['nodes'])
        require(le<=2e-11 and ce<=2e-10,'independent accepted circle oracle failed')
        require(all(c['closed'] and all(n['curvature_mask']==0 for n in c['nodes'])
                    for s in frames for c in s['curves']),'closed-loop fixture ineffective')
        return dict(frames=len(frames),length_error=le,curvature_error=ce,
                    accepted_length_budget=2e-11,accepted_curvature_budget=2e-10)
    check('CRV13_accepted_closed_circle_length_curvature_oracle',circle)
    def cadence():
        directory=run('cadence',changes={'output2/dcycle':5})
        data=read_run(directory);last=max(data['curve_times'])
        require(set(data['curve_times'])=={0,last}|set(range(5,last+1,5)),'snapshot cadence mismatch')
        require(len(list((directory/'curves').glob('*.crv')))==len(data['curve_times']),'duplicate terminal snapshot')
        return dict(saved_cycles=sorted(data['curve_times']))
    check('CRV07_sparse_cadence_and_final_exact_replay',cadence)
    def restarted():
        full=read_run(cases['planar'])
        stop=min(g for (g,n),o in full['curve_nodes'].items() if o['fields']['shock_region']==1)
        seed=run('restart_seed',changes={'time/nlim':stop,'output1/dcycle':stop,'output1/id':'checkpoint'})
        journal=work/'restart-journal';shutil.copytree(seed/'shock_events',journal)
        checkpoints=list((seed/'rst').glob('*.rst'));require(checkpoints,'missing paired checkpoint')
        checkpoint=max(checkpoints,key=lambda p:int(p.name.split('.')[-2]))
        # Select by native checkpoint cycle using retained seed final output clock.
        source=read_run(seed);require(max(source['curve_times'])==stop,'restart seed clock mismatch')
        from curve_shock_history import read_trace
        trace=read_trace(seed/'shock_events/curves'/('trace'+str(stop)+'.rank0.csv'))
        require(any(o['phase']==1 for o in trace.values()),'checkpoint is not in an active encounter')
        suffix=run('restart_resume',restart=checkpoint,overrides=['time/nlim=-1','shock_tracking/directory='+str(journal)])
        zero_journal=work/'zero-journal';shutil.copytree(seed/'shock_events',zero_journal)
        zero=run('restart_zero',restart=checkpoint,overrides=['time/nlim='+str(stop),'shock_tracking/directory='+str(zero_journal)])
        full=read_run(cases['planar']);continued=read_run(suffix);z=read_run(zero)
        for key in ('curve_nodes','particles'):
            require(normalized(source[key]|continued[key])==normalized(full[key]),'restart continuation differs '+key)
            require(normalized(z[key])==normalized({k:v for k,v in source[key].items() if k[0]==stop}),'zero-step differs '+key)
        for s in [*continued['geometry'],*z['geometry']]:
            expected=next(f for f in full['geometry'] if f['cycle']==s['cycle'])
            require(normalized(s['curves'])==normalized(expected['curves']),'restart topology differs')
        cases['restart']=suffix;cases['zero']=zero
        return dict(active_restart_cycle=stop,continuation_exact=True,zero_step_exact=True)
    check('CRV08_active_and_zero_step_restart_exact',restarted)
    def empty():
        directory=run('empty',changes={'material_curves/ncurves':0,'tracer_seed_initial/count_per_event':0,'time/nlim':2})
        data=read_run(directory)
        require(not data['curve_nodes'] and not data['particles'] and len(data['geometry'])==3,'empty metadata missing')
        require(all(s['field_names'] and s['curve_count']==s['node_count']==0 for s in data['geometry']),'empty rank schema missing')
        return dict(empty_clocks=3)
    check('CRV09_empty_populations_explicit_schema',empty)
    def bits_zero():
        directory=run('bits_zero',changes={'material_curves/curve_id_bits':0,'curve1/id':1,'time/nlim':4})
        data=read_run(directory)
        require(data['curve_id_layout']==(32,0,True),'zero-bit layout lost')
        return dict(full_32bit_layout=True,samples=len(data['curve_nodes']))
    check('CRV10_zero_curve_bits',bits_zero)
    def guards():
        run('guard_disabled',changes={'shock_tracking/enabled':'false'},reject='shock fields require')
        run('guard_missing',changes={'material_curves/enabled':'false','output2/source':'particles',
            'output4/file_type':'crv','output4/dcycle':1},reject='curve output requires')
        return dict(disabled_tracking_rejected=True,missing_curves_rejected=True)
    check('CRV11_configuration_guards',guards)
    check('CRV12_short_accepted_forced_box_exact_parity',lambda:parity('forced'))
    def dt_restart():
        options={'output2/dcycle':None,'output2/dt':.05}
        full=read_run(run('dt_full','smooth',options))
        times=list(full['curve_times'].values())
        require(len(times)==9 and times[0]==0 and times[-1]==.4,'time cadence changed')
        for stop in (18,20):
            seed=run('dt_seed_'+str(stop),'smooth',dict(options,**{'time/nlim':stop,'output1/dt':1}))
            checkpoint=max((seed/'rst').glob('*.rst'),key=lambda p:int(p.name.split('.')[-2]))
            resumed=run('dt_resume_'+str(stop),'smooth',restart=checkpoint,overrides=['time/nlim=-1'])
            before,after=read_run(seed),read_run(resumed)
            joined=before['curve_nodes']|after['curve_nodes']
            require(all(normalized(joined.get(k))==normalized(v) for k,v in full['curve_nodes'].items()),
                    'dt restart lost/changed a scheduled sample')
            require(not before['curve_nodes'].keys()&after['curve_nodes'].keys(),'dt restart duplicated identity')
        return dict(saved_times=times,on_and_off_cadence_restart_exact=True)
    check('CRV14_time_cadence_restart_and_no_duplicate_identity',dt_restart)
    def immutable():
        directory=run('immutable','smooth',{'time/nlim':2})
        command=report['runs'][-1]['command']
        paths=list((directory/'curves').glob('*.crv'));before={p.name:sha(p) for p in paths}
        def replay(name):
            with (directory/(name+'.log')).open('w') as log:
                code=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=90).returncode
            report['runs'].append(dict(name=name,command=command,exit=code,
                log_sha256=sha(directory/(name+'.log')),**provenance(args.prescribed_executable)))
            save();return code
        require(replay('identical_replay')==0,'identical replay rejected')
        require(before=={p.name:sha(p) for p in paths},'identical replay modified chunks')
        path=paths[0];path.write_bytes(path.read_bytes()+b'conflicting output');damaged=sha(path)
        require(replay('conflicting_replay')!=0 and 'existing file differs' in (directory/'conflicting_replay.log').read_text(),
                'conflicting replay accepted')
        require(sha(path)==damaged,'conflicting output silently overwritten')
        return dict(identical_replay_exact=True,conflicting_output_rejected=True)
    check('CRV15_immutable_replay_no_silent_overwrite',immutable)
    def ascii_coexistence():
        directory=run('explicit_ascii',changes={'output4/file_type':'curves','output4/id':'shape','output4/dcycle':4})
        current=read_run(directory);full=read_run(cases['planar'])
        require(normalized(current['curve_nodes'])==normalized(full['curve_nodes']),'explicit ASCII changed node samples')
        require([normalized(s['curves']) for s in current['geometry']]==
                [normalized(s['curves']) for s in full['geometry']],'explicit ASCII/CRV integration differs')
        require(list((directory/'curves').glob('*.txt')),'explicit ASCII path not exercised')
        return dict(explicit_ASCII_preserved=True,matching_geometry_validated=True)
    check('CRV16_explicit_ASCII_and_unified_reader_coexistence',ascii_coexistence)


    save();return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
