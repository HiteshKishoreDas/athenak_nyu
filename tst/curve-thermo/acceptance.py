"""Laptop checks for shared particle/curve thermo output. Retain every native run.

Use the existing built-in and material_curve_test builds; --workdir must be new.
Oracles and tolerance (2e-12) fixed before executing tests, inherited from LT/CSH.
"""
import argparse
import json
import math
import re
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'), str(ROOT/'tst/shocks-mc')]
from common import edit_input, provenance, sha, thermo
from tracer_thermo_history import read_curves
from curve_shock_history import read_trace, read_history
import material_curves

BUDGET = 2e-12
FIELDS = 'rho,p,t,s,eint,vx,vy,vz,vmag,cs,mach,rho'
NAMES = ['density','pressure','temperature','entropy','internal_energy',
         'v1','v2','v3','vmag','sound_speed','mach']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--prescribed-executable', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve(); work.mkdir(parents=True,exist_ok=False)
    report = dict(schema='curve-thermo-1', runs=[], checks=[], budget=BUDGET,
                  hardware='local Serial/double', scientific_approval='pending',
                  runner_sha256=sha(__file__))

    def save():
        (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')

    def check(name, action):
        row = dict(id=name)
        try:
            row.update(action() or {}); row['status']='pass'
        except Exception as exc:
            row.update(status='fail',error=repr(exc))
        report['checks'].append(row); save()

    def run(name, changes=None, fixture='smooth', restart=None, overrides=(), reject=None):
        directory=work/name; directory.mkdir()
        binary=(args.executable if fixture in ('planar','forced') else args.prescribed_executable).resolve()
        path=directory/'input.athinput'
        options={'shock_tracking/enabled':'false','output3/file_type':'prtcl_thermo_history',
                 'output3/source':'curves','output3/variables':FIELDS,'output3/dcycle':1}
        options.update(changes or {})
        options={key:value for key,value in options.items() if value is not None}
        path.write_text(edit_input((ROOT/('inputs/tests/shocks-curves/'+fixture+'.athinput')).read_text(),options))
        command=([str(binary),'-r',str(restart)] if restart else [str(binary),'-i',str(path)])+list(overrides)
        started=time.monotonic()
        with (directory/'run.log').open('w') as log:
            status=subprocess.run(command,cwd=directory,stdout=log,
                                  stderr=subprocess.STDOUT,timeout=90).returncode
        report['runs'].append(dict(name=name,command=command,exit=status,
            seconds=time.monotonic()-started,input_used=not bool(restart),input_sha256=sha(path),
            restart_sha256=sha(restart) if restart else None,log_sha256=sha(directory/'run.log'),
            **provenance(binary)))
        save()
        log=(directory/'run.log').read_text()
        require((status!=0 and reject.lower() in log.lower()) if reject else status==0,
                name+': unexpected native exit/diagnostic; see run.log')
        return directory

    particles={'particles/particle_type':'lagrangian_tracer','particles/pusher':'lagrangian_tracer',
        'particles/random_seed':11,'tracer_seed_initial/id':1,'tracer_seed_initial/count_per_event':128,
        'tracer_seed_initial/start_time':0,'tracer_seed_initial/seed':903,
        'tracer_seed_initial/region':'box','tracer_seed_initial/weight':'volume',
        'tracer_seed_initial/x1min':-.4,'tracer_seed_initial/x1max':.4,
        'tracer_seed_initial/x2min':-.4,'tracer_seed_initial/x2max':.4}

    def oracle(rows, acceleration=False):
        maximum=0
        for sample in rows.values():
            t=sample['time']; v=[.2+.3*t,-.1+.2*t,0] if acceleration else [.2,-.1,0]
            gamma=5/3; rho=1; pressure=(4-max(map(abs,v)))**2/gamma
            speed=math.sqrt(sum(x*x for x in v)); cs=math.sqrt(gamma*pressure/rho)
            expected=[rho,pressure,pressure/rho,math.log(pressure/rho**gamma),
                      pressure/(gamma-1),*v,speed,cs,speed/cs]
            require(list(sample['fields'])==NAMES,'field alias/order/dedup mismatch')
            maximum=max(maximum,max(abs(sample['fields'][k]-value) for k,value in zip(NAMES,expected)))
        require(maximum<=BUDGET,'independent prescribed thermo oracle failed: '+str(maximum))
        return maximum

    def curve_only():
        directory=run('curve_only',{'curve1/id':255,'problem/flow':'acceleration',
                                   'problem/a1':.3,'problem/a2':.2})
        rows,times=read_curves(directory/'curve_thermo_history')
        require(rows and -(1<<31) <= min(p for g,p in rows) and max(p for g,p in rows)<0,'signed high-bit node IDs not exercised')
        require(not (directory/'prtcl_thermo_history').exists(),'curve-only output requires particles')
        require(not (directory/'shock_events').exists(),'thermo unexpectedly depends on shocks')
        error=oracle(rows,True)
        return dict(observations=len(rows),cycles=len(times),max_thermo_error=error,
                    minimum_node_id=min(p for g,p in rows))
    check('CT01_curve_only_shocks_off_exact_u32',curve_only)

    def both():
        directory=run('both',dict(particles,**{'output3/source':'both','problem/flow':'acceleration',
                                              'problem/a1':.3,'problem/a2':.2}))
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'curve/particle clocks differ')
        return dict(curve_observations=len(curves),particle_observations=len(regular),
                    curve_error=oracle(curves,True),particle_error=oracle(regular,True))
    check('CT02_both_shared_fields_and_clock',both)

    def passivity():
        options=dict(particles,**{'problem/flow':'acceleration','problem/a1':.3,'problem/a2':.2})
        legacy=run('legacy',dict(options,**{'output3/source':None}))
        both=work/'both'
        require(not (legacy/'curve_thermo_history').exists(),'default particle output changed')
        a=next((both/'prtcl_thermo_history').glob('*.thp'))
        b=next((legacy/'prtcl_thermo_history').glob('*.thp'))
        require(a.read_bytes()==b.read_bytes(),'curve output changed legacy particle bytes')
        require((both/'material.curves.final.txt').read_bytes()==
                (legacy/'material.curves.final.txt').read_bytes(),'curve output changed final curve state')
        return dict(particle_binary_and_curve_state_byte_equal=True)
    check('CT03_passive_legacy_output',passivity)

    def mc():
        options=dict(particles,**{'particles/particle_type':'lagrangian_mc',
                                 'particles/pusher':'lagrangian_mc','output3/source':'both'})
        directory=run('mc_and_curves',options)
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'MC/curve output clocks differ')
        return dict(curve_error=oracle(curves),particle_error=oracle(regular))
    check('CT04_MC_and_curves_thermo_without_shocks',mc)

    def amr():
        directory=run('amr',dict(particles,**{'output3/source':'both'}),fixture='amr')
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'AMR clocks differ')
        states=[material_curves.read(p) for p in directory.glob('material.curves.*.txt')]
        counts=[len(s['amr']['leaves']) for s in states]
        require(max(counts)>min(counts) and len(material_curves.read(directory/'material.curves.final.txt')['amr']['leaves'])<max(counts),
                'live refinement/coarsening not exercised')
        earliest={}
        for (g,p),o in sorted(curves.items()): earliest.setdefault(p,o)
        newborn=[o for o in earliest.values() if o['birth']>0]
        require(newborn and all(o['time']==o['birth'] for o in newborn),'insertion birth samples missing/invented')
        return dict(curve_error=oracle(curves),particle_error=oracle(regular),
                    inserted_nodes=len(newborn),leaf_count_range=[min(counts),max(counts)],
                    curve_observations=len(curves),particle_observations=len(regular))
    check('CT05_live_AMR_insertion_births',amr)

    def restart():
        seed=run('restart_seed',dict(particles,**{'output3/source':'both','time/nlim':20}),fixture='amr')
        checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
        resumed=run('restart_resume',fixture='amr',restart=checkpoint,
                    overrides=('time/nlim=-1',))
        reference,rt=read_curves(work/'amr/curve_thermo_history')
        prefix,st=read_curves(seed/'curve_thermo_history')
        suffix,et=read_curves(resumed/'curve_thermo_history')
        require(not (prefix.keys() & suffix.keys()),'restart duplicated curve observation identity')
        require(dict(prefix,**{})|suffix==reference,'restart curve thermo differs from uninterrupted run')
        regular,pt=thermo(work/'amr'); pr,_=thermo(seed); sr,_=thermo(resumed)
        require(not (pr.keys() & sr.keys()) and pr|sr==regular,'restart particle thermo differs')
        require(st|et==rt,'restart clocks changed')
        return dict(curve_observations=len(reference),particle_observations=len(regular),
                    exact_continuation=True,restart_cycle=max(st))
    check('CT06_live_AMR_restart_exact',restart)

    def shocks():
        directory=run('planar',{'output2/dt':-1,'output2/dcycle':0,
            'output3/source':'both','shock_tracking/enabled':'true',
            'shock_tracking/directory':str(work/'planar/shock_events')},fixture='planar')
        rows,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'shock/thermo populations have differing clocks')
        signals={}
        for path in (directory/'shock_events/curves').glob('trace*.rank*.csv'):
            for p,s in read_trace(path).items(): signals[s['generation'],p]=s
        require(rows.keys()==signals.keys(),'curve thermo/shock sample identities differ')
        require(any(s['inside'] for s in signals.values()) and
                read_history(directory/'shock_events/curves'),'planar fixture lacks shock encounters')
        error=0; position_error=0
        for key,s in rows.items():
            o=signals[key]
            require(s['time']==o['time'] and s['birth']==o['birth'],
                    'shock/thermo time/birth mismatch')
            # Thermo uses authoritative continuous lifts; shock evidence is wrapped.
            for x,y,length in zip(s['x'],o['x'],(1,.1,.1)):
                position_error=max(position_error,abs((x-y+length/2)%length-length/2))
            for field,raw in [('density','rho'),('pressure','pressure'),('temperature','temperature')]:
                error=max(error,abs(s['fields'][field]-o[raw]))
        require(error<=BUDGET and position_error<=BUDGET,'shock/thermo sampling differs')
        return dict(observations=len(rows),max_shared_sample_error=error,
                    max_periodic_position_error=position_error)
    check('CT07_shock_thermo_synchronization',shocks)

    def cadence():
        directory=run('cadence',{'output3/dcycle':5})
        _,clocks=read_curves(directory/'curve_thermo_history')
        final=material_curves.read(directory/'material.curves.final.txt')['cycle']
        require(set(clocks)=={0,final}|set(range(5,final+1,5)),'thermo cadence ignored')
        # Exact final-cycle replay is native (dcycle plus Finalize); one immutable chunk remains.
        return dict(saved_cycles=sorted(clocks),final_cycle=final)
    check('CT08_independent_cadence_final_replay',cadence)

    def guards():
        for name,change,reason in [
            ('bad_source',{'output3/source':'wrong'},'source must be'),
            ('missing_curves',{'material_curves/enabled':'false'},'requires material_curves'),
            ('missing_particles',{'output3/source':'both'},'particle thermo source requires'),
            ('bad_field',{'output3/variables':'magical'},'unknown field'),
            ('bad_mhd',{'output3/variables':'bmag'},'requires MHD'),
            ('bad_scalar',{'output3/variables':'scalar0'},'out of range')]:
            run(name,change,reject=reason)
        return dict(explicit_rejections=6)
    check('CT09_input_guards',guards)

    def reader_guards():
        source=next((work/'curve_only/curve_thermo_history').glob('*.thp'))
        directory=work/'reader_corruption';directory.mkdir(); path=directory/source.name
        text=source.read_bytes()
        variants=[text[:-1],text[:32]+(9).to_bytes(4,sys.byteorder,signed=True)+text[36:],
                  text.replace(b'birth\na\nlevel',b'birth\na\nbirth',1)]
        for changed in variants:
            require(changed != text,'ineffective reader corruption fixture')
            path.write_bytes(changed)
            try: read_curves(directory)
            except ValueError: pass
            else: raise ValueError('reader accepted corrupt chunk')
        return dict(rejected_corruptions=len(variants))
    check('CT10_reader_integrity',reader_guards)

    def nonuniform():
        # Same accepted nonuniform-density fixture as material-curves/production.py.
        directory=run('nonuniform',dict(particles,**{'output3/source':'both',
                                                    'problem/density_amplitude':.5}))
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'nonuniform population clocks differ')
        maximum=0; pairs={}; comparisons=0
        for source,rows in [('curves',curves),('particles',regular)]:
            for (cycle,p),sample in rows.items():
                x=-2+(sample['x'][0]+2)%4
                cell=min(31,math.floor((x+2)/.125)); center=-2+(cell+.5)*.125
                rho=1+.5*math.sin(center); pressure=rho*3.8**2/(5/3)
                maximum=max(maximum,abs(sample['fields']['density']-rho),
                             abs(sample['fields']['pressure']-pressure))
                key=cycle,cell
                if key in pairs:
                    require(sample['fields']==pairs[key],'same-cell thermo fields disagree')
                    comparisons+=1
                pairs[key]=sample['fields']
        require(maximum<=BUDGET,'containing-cell spatial oracle failed: '+str(maximum))
        require(comparisons>0,'same-cell comparison not exercised')
        return dict(max_spatial_thermo_error=maximum,same_cell_comparisons=comparisons)
    check('CT11_nonuniform_containing_cell_oracle',nonuniform)

    def immutable():
        directory=run('immutable')
        command=report['runs'][-1]['command']
        def replay(name):
            with (directory/(name+'.log')).open('w') as log:
                status=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=90).returncode
            report['runs'].append(dict(name=name,command=command,exit=status,
                directory=str(directory),log_sha256=sha(directory/(name+'.log')),
                **provenance(args.prescribed_executable)))
            save(); return status
        before={p.name:sha(p) for p in (directory/'curve_thermo_history').glob('*.thp')}
        require(replay('identical_replay')==0,'identical curve replay rejected')
        require(before=={p.name:sha(p) for p in (directory/'curve_thermo_history').glob('*.thp')},
                'identical replay modified/duplicated chunks')
        path=next((directory/'curve_thermo_history').glob('*.cycle0.rank0of1.bits*.thp'))
        path.with_suffix('.original').write_bytes(path.read_bytes())
        path.write_bytes(path.read_bytes()+b'deliberate conflicting output')
        damaged=sha(path)
        require(replay('conflicting_replay')!=0 and 'existing file differs' in
                (directory/'conflicting_replay.log').read_text(),'conflicting output not rejected')
        require(sha(path)==damaged,'conflicting output silently overwritten')
        return dict(identical_replay_byte_equal=True,conflicting_replay_rejected=True)
    check('CT12_immutable_output_replay',immutable)

    def time_cadence():
        directory=run('time_cadence',dict(particles,**{'output3/source':'both',
            'output3/dcycle':None,'output3/dt':.05}))
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(ct==pt,'dt-cadenced clocks differ')
        times=[t for c,t in sorted(ct.items())]
        require(times[0]==0 and times[-1]==.4 and len(times)==9,'dt cadence/final output mismatch')
        # Native outputs occur at the first completed timestep reaching each threshold.
        return dict(saved_times=times,curve_error=oracle(curves),particle_error=oracle(regular))
    check('CT13_shared_time_cadence',time_cadence)

    def forced():
        # Output-only changes to the accepted 32^3 forced-box showcase.
        changes={'shock_tracking/enabled':'true','output2/source':'both',
            'output2/variables':FIELDS,'output3/source':None,'output3/variables':None,
            'output3/file_type':'bin','output3/dcycle':None}
        both=run('forced_both',changes,fixture='forced')
        legacy=run('forced_legacy',dict(changes,**{'output2/source':None}),fixture='forced')
        curves,ct=read_curves(both/'curve_thermo_history'); regular,pt=thermo(both)
        require(ct==pt,'forced population clocks differ')
        a=next((both/'prtcl_thermo_history').glob('*.thp'))
        b=next((legacy/'prtcl_thermo_history').glob('*.thp'))
        require(a.read_bytes()==b.read_bytes(),'curve thermo changed forced particle output')
        require((both/'shock_turb_forced.curves.final.txt').read_bytes()==
                (legacy/'shock_turb_forced.curves.final.txt').read_bytes(),
                'curve thermo changed forced curve geometry')
        require(json.dumps(read_history(both/'shock_events/curves'),sort_keys=True)==
                json.dumps(read_history(legacy/'shock_events/curves'),sort_keys=True),
                'thermo output changed forced curve encounters')
        return dict(curve_observations=len(curves),particle_observations=len(regular),
                    sampled_cycles=len(ct),particle_bytes_geometry_and_events_equal=True)
    check('CT14_native_forced_box_coexistence',forced)

    def scalar():
        directory=run('passive_scalar',{'hydro/nscalars':1,'time/nlim':0,
            'output2/dt':-1,'output2/dcycle':0,'output3/source':'both',
            'output3/variables':'rho,p,scalar0','curve1/x1':.2},fixture='planar')
        curves,_=read_curves(directory/'curve_thermo_history'); regular,_=thermo(directory)
        values=set()
        for rows in (curves,regular):
            for sample in rows.values():
                center=(math.floor(sample['x'][0]*512)+.5)/512
                expected=1 if center<.3 else 0
                require(sample['fields']['scalar0']==expected,'initial passive scalar oracle mismatch')
                values.add(expected)
        require(values=={0,1},'scalar discontinuity not sampled')
        return dict(curve_observations=len(curves),particle_observations=len(regular),
                    scalar_left_right_exact=True)
    check('CT15_shared_passive_scalar',scalar)

    def dt_restart():
        reference,rt=read_curves(work/'time_cadence/curve_thermo_history')
        regular,_=thermo(work/'time_cadence'); cases=[]
        for stop in (18,20):  # Off cadence, then on cadence with repeated final output.
            seed=run('dt_seed_'+str(stop),dict(particles,**{'output3/source':'both',
                'output3/dcycle':None,'output3/dt':.05,'output1/dt':1,'time/nlim':stop}))
            checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
            resumed=run('dt_resume_'+str(stop),restart=checkpoint,overrides=('time/nlim=-1',))
            prefix,_=read_curves(seed/'curve_thermo_history')
            suffix,_=read_curves(resumed/'curve_thermo_history'); joined=prefix|suffix
            pr,_=thermo(seed); sr,_=thermo(resumed); joined_particles=pr|sr
            require(not (prefix.keys() & suffix.keys()),'dt restart duplicated curve identity')
            require(all(joined.get(k)==v for k,v in reference.items()),'dt restart lost/changed scheduled curve sample')
            require(all(joined_particles.get(k)==v for k,v in regular.items()),'dt restart lost/changed scheduled particle sample')
            require(all(g==stop for g,p in joined.keys()-reference.keys()),'unexpected extra dt restart sample')
            cases.append(dict(stop_cycle=stop,scheduled_samples_exact=True))
        return dict(cases=cases)
    check('CT16_restart_preserves_time_cadence',dt_restart)

    def overdue():
        directory=run('overdue',dict(particles,**{'output3/source':'both',
            'output3/dcycle':None,'output3/dt':.005,'output1/dcycle':18,'time/nlim':18}))
        def clock(path):
            header=path.read_bytes().split(b'<par_end>',1)[0].decode()
            block=header.split('<output3>',1)[1].split('<',1)[0]
            match=re.search(r'^last_time\s*=\s*(\S+)',block,re.M)
            require(match is not None,'checkpoint lacks output clock')
            return float(match[1])
        checkpoints=sorted((directory/'rst').glob('*.rst'))
        require(len(checkpoints)==3,'fixture lacks initial, scheduled and terminal checkpoints')
        ordinary,terminal=checkpoints[-2:]
        scheduled,final=clock(ordinary),clock(terminal)
        require(scheduled==final,'same-cycle Finalize advanced overdue output clock twice')
        require(abs(final-.09)<=BUDGET,'unexpected native overdue-threshold schedule')
        resumed=run('overdue_zero_step',restart=terminal)
        last=sorted((resumed/'rst').glob('*.rst'))[-1]
        require(clock(last)==final,'zero-step restart advanced persisted output clock')
        return dict(scheduled_checkpoint_last_time=scheduled,final_checkpoint_last_time=final,
                    output_dt=.005,fluid_dt=.0125,zero_step_clock_exact=True)
    check('CT17_overdue_threshold_replay',overdue)

    def empty():
        directory=run('empty',dict(particles,**{'output3/source':'both',
            'material_curves/ncurves':0,'tracer_seed_initial/count_per_event':0}))
        curves,ct=read_curves(directory/'curve_thermo_history'); regular,pt=thermo(directory)
        require(not curves and not regular and ct==pt and len(ct)>1,'empty output batches failed')
        return dict(empty_chunks=len(ct),empty_particle_and_node_records=True)
    check('CT18_empty_populations',empty)
    save()
    print(json.dumps(report['checks'],indent=2))
    return 0 if all(c['status']=='pass' for c in report['checks']) else 1


if __name__=='__main__':
    raise SystemExit(main())
