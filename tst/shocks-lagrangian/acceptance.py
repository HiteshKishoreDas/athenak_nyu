"""Native continuous tracer qualification. All raw failures retained; stdlib only.

Run after both builds finish. Existing CSH gates additionally cover curve histories.
Position budget 2e-12 is inherited from affine curve controls, never fit to outcomes.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'), str(ROOT/'tst/shocks-mc'),
                str(ROOT/'tst/shocks-curves')]
from common import edit_input, moment_oracle, provenance, sha, thermo, uniform_changes
from integration import traces, canonical
from acceptance import density_fronts, snapshots
from shock_history import read_history, read_suspended

BUDGET = 2e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--prescribed-executable', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve(); work.mkdir(parents=True, exist_ok=False)
    report = dict(schema='lagrangian-shocks-1', runs=[], checks=[],
                  scientific_approval='pending', hardware='local Serial/double',
                  position_budget=BUDGET, runner_sha256=sha(__file__))

    def save():
        (work/'report.json').write_text(json.dumps(report, indent=2, allow_nan=True)+'\n')

    def check(name, action):
        row = dict(id=name)
        try:
            row.update(action() or {}); row['status']='pass'
        except Exception as exc:
            row.update(status='fail', error=repr(exc))
        report['checks'].append(row); save()

    def run(name, changes=None, prescribed=False, fixture=None, restart=None,
            overrides=(), reject=None):
        directory=work/name; directory.mkdir()
        binary=(args.prescribed_executable if prescribed else args.executable).resolve()
        fixture=fixture or 'inputs/tests/shock_lagrangian.athinput'
        options={'shock_tracking/directory':str(directory/'shock_events'),
                 'output3/dcycle':1, 'output3/data_format':' %24.17e'}
        options.update(changes or {})
        path=directory/'input.athinput'
        path.write_text(edit_input((ROOT/fixture).read_text(),options))
        command=([str(binary),'-r',str(restart)] if restart else
                 [str(binary),'-i',str(path)])+list(overrides)
        begin=time.monotonic()
        with (directory/'run.log').open('w') as log:
            try:
                status=subprocess.run(command,cwd=directory,stdout=log,
                                      stderr=subprocess.STDOUT,timeout=60).returncode
            except subprocess.TimeoutExpired:
                status='timeout'
        report['runs'].append(dict(name=name, command=command, exit=status,
            seconds=time.monotonic()-begin, input_sha256=None if restart else sha(path),
            generated_input_used=not bool(restart),
            restart_sha256=sha(restart) if restart else None,
            log_sha256=sha(directory/'run.log'), **provenance(binary)))
        save()
        if reject:
            require(status != 0 and reject.lower() in (directory/'run.log').read_text().lower(),
                    name+': missing explicit rejection: '+reject)
        else:
            require(status==0,name+': failed native simulation; inspect run.log')
        return directory

    def position_error(directory, velocity, acceleration=(0,0,0), shear=0, length=1, low=0):
        rows,times=thermo(directory); first=min(times); final=max(times)
        initial={p:o for (g,p),o in rows.items() if g==first}
        terminal={p:o for (g,p),o in rows.items() if g==final}
        require(initial and initial.keys()==terminal.keys(),'tracer ID loss or no tracers')
        lengths=[length]*3 if isinstance(length,(int,float)) else length
        errors=[]; expected=[]; observed=[]; changed=0
        for (g,p),o in rows.items():
            t=o['time']-times[first]; x=initial[p]['x']; target=[]
            for d in range(3):
                value=x[d]+velocity[d]*t+.5*acceleration[d]*t*t
                if d==0:value+=shear*(x[1]*t+.5*velocity[1]*t*t)
                target.append(low+(value-low)%lengths[d])
            diff=[(o['x'][d]-target[d]+lengths[d]/2)%lengths[d]-lengths[d]/2 for d in range(3)]
            errors.append(max(map(abs,diff)))
            if g==final:
                changed += o['gid'] != initial[p]['gid']
                observed.append(o['x']);expected.append(target)
        error=max(errors)
        require(error<=BUDGET,'continuous trajectory error '+str(error)+' exceeds '+str(BUDGET))
        return dict(max_error=error, observations=len(rows), particles=len(initial),
                    changed_leaf=changed, times=list(times.values()),
                    initial=[initial[p]['x'] for p in sorted(initial)],
                    final=observed, expected=expected)

    uniform=uniform_changes(count=128,blocks=4)
    uniform.update({'particles/particle_type':'lagrangian_tracer',
                    'particles/pusher':'lagrangian_tracer'})
    def uniform_test(three=False):
        options=dict(uniform)
        if three:options.update({'mesh/nx3':4,'meshblock/nx3':4,
                                'mesh/ix3_bc':'periodic','mesh/ox3_bc':'periodic',
                                'problem/wl':.4,'problem/wr':.4})
        directory=run('uniform3d' if three else 'uniform2d',options)
        result=position_error(directory,(1,0,.4 if three else 0),length=(1,.1,.1))
        if three:
            samples,_=thermo(directory);last={};wraps=0
            for (g,p),o in sorted(samples.items()):
                if p in last and o['x'][2]<last[p]-.05:wraps+=1
                last[p]=o['x'][2]
            require(wraps>0,'3D fixture never crossed the z periodic seam')
            result['z_wraps']=wraps
        require(result['changed_leaf']>0,'fixture never crossed a block boundary')
        return result
    check('LT01_uniform_periodic2d',uniform_test)
    check('LT01_uniform_periodic3d',lambda:uniform_test(True))

    smooth={'shock_tracking/enabled':'false', 'particles/particle_type':'lagrangian_tracer',
        'particles/pusher':'lagrangian_tracer', 'particles/random_seed':11,
        'tracer_seed_initial/id':1, 'tracer_seed_initial/count_per_event':128,
        'tracer_seed_initial/start_time':0, 'tracer_seed_initial/seed':903,
        'tracer_seed_initial/region':'box', 'tracer_seed_initial/weight':'volume',
        'tracer_seed_initial/x1min':-.4, 'tracer_seed_initial/x1max':.4,
        'tracer_seed_initial/x2min':-.4, 'tracer_seed_initial/x2max':.4,
        'output3/file_type':'prtcl_thermo_history', 'output3/dcycle':1,
        'output3/variables':'rho,p,v1,v2,v3'}
    def affine(name, flow, acceleration=(0,0,0), shear=0, amr=False):
        options=dict(smooth, **{'problem/flow':flow})
        for d,a in enumerate(acceleration):options['problem/a'+str(d+1)]=a
        directory=run(name,options,prescribed=True,
                      fixture='inputs/tests/shocks-curves/'+('amr' if amr else 'smooth')+'.athinput')
        result=position_error(directory,(.2,-.1,0),acceleration,shear,length=4,low=-2)
        if amr:
            counts=[len(s['amr']['leaves']) for _,s in sorted(snapshots(directory).items())]
            # Prescribed fixture's basename is material.
            if not counts:
                import material_curves
                counts=[len(material_curves.read(p)['amr']['leaves'])
                        for p in sorted(directory.glob('material.curves.*.txt'))]
            require(counts and max(counts)>min(counts),'AMR did not change the leaf tree')
            final=__import__('material_curves').read(directory/'material.curves.final.txt')
            require(len(final['amr']['leaves'])<max(counts),'AMR did not coarsen')
            result['leaf_counts']=counts
        return result
    check('LT02_acceleration_RK_stage_time',lambda:affine('acceleration','acceleration',(.3,.2,0)))
    check('LT02_affine_shear',lambda:affine('shear','shear',shear=.2))
    check('LT04_refine_coarsen_position',lambda:affine('amr','translation',amr=True))

    def mixed_amr():
        options=dict(smooth, **{'problem/flow':'translation','problem/v1':.97,
            'problem/v2':.97,'problem/refine_time':0,
            'mesh_refinement/refinement_interval':1,'mesh_refinement/ncycle_check':1,'problem/refine_x1max':0,
            'tracer_seed_initial/x1min':-.124,'tracer_seed_initial/x1max':-.001,
            'tracer_seed_initial/x2min':-1.124,'tracer_seed_initial/x2max':-1.001})
        directory=run('mixed_amr',options,prescribed=True,
                      fixture='inputs/tests/shocks-curves/amr.athinput')
        result=position_error(directory,(.97,.97,0),length=4,low=-2)
        import material_curves
        states=[material_curves.read(p) for p in directory.glob('material.curves.*.txt')]
        levels=[{leaf['level'] for leaf in state['amr']['leaves']} for state in states]
        require(any(len(level)>1 for level in levels),'mixed-level AMR not actually exercised')
        require(result['changed_leaf']>0,'tracer never changed leaf')
        require('Lagrangian containing-leaf fallback' in (directory/'run.log').read_text(),
                'mixed diagonal fixture never exercised the containing-leaf fallback')
        result['mixed_levels_exercised']=True
        result['containing_leaf_fallback_exercised']=True
        return result
    check('LT04_mixed_level_diagonal_crossing',mixed_amr)


    def planar():
        # Same accepted RH states/resolution as original fixture. Restrict births
        # to upstream slab that is swept within its accepted duration.
        directory=run('planar',{'output3/dcycle':1})
        rows,times=thermo(directory); observed=traces(directory/'shock_events')
        require(rows.keys()==observed.keys(),'shock/thermo sample identity mismatch')
        require(all(o['time']==rows[key]['time'] and tuple(o['x'])==rows[key]['x']
                    for key,o in observed.items()),'shock/thermo synchronization mismatch')
        events=read_history(directory/'shock_events');require(events,'no planar encounters')
        counts=Counter(e['particle_id'] for e in events)
        require(max(counts.values())==1,'duplicate planar passage')
        require(all(e['status']==0 for e in events),'unexpected partial planar event')
        fronts=density_fronts(directory);front={g:x for g,t,x in fronts};covered=0;brackets=[]
        crossed={}
        all_ids={p for g,p in observed}
        for p in all_ids:
            samples=[(o['time'],o['x'][0]-front[g]) for (g,q),o in sorted(observed.items())
                     if q==p and g in front]
            crossings=[(a[0],b[0]) for a,b in zip(samples,samples[1:]) if a[1]>0>=b[1]]
            require(len(crossings)==1,'seeded particle lacks exactly one independent density crossing')
            crossed[p]=crossings[0]
        require(counts.keys()==crossed.keys(),'particle shock events missing/excess independent crossings')
        for e in events:
            crossing=crossed[e['particle_id']]
            require(e['pre']['time']<=crossing[1] and e['post']['time']>=crossing[0],
                    'event envelope excludes independent density crossing bracket')
            covered+=1;brackets.append(crossing)
        return dict(events=len(events),density_brackets_covered=covered,density_brackets=brackets,
            density_fronts=fronts,mach=[e['representative']['mach_pressure'] for e in events])
    check('LT05_PHY08_INT06_RH_passage',planar)

    def restart():
        base=work/'planar';ref,_=thermo(base);events=read_history(base/'shock_events')
        observed=traces(base/'shock_events');phases={}
        for name,phase in [('outside',0),('active',1),('pending',2)]:
            generations=[g for (g,p),o in observed.items() if g>0 and o['phase']==phase]
            require(generations,'missing standalone phase '+name);g=min(generations)
            seed=run('seed_'+name,{'time/nlim':g,'output1/dt':1})
            cp=sorted((seed/'rst').glob('*.rst'))[-1]
            resumed=run('resume_'+name,restart=cp,
                        overrides=('time/nlim=-1','job/basename=resumed_'+name))
            actual,_=thermo(resumed)
            require(actual=={k:v for k,v in ref.items() if k[0]>g},'native restart values differ')
            require(canonical(read_history(seed/'shock_events'))==canonical(events),
                    'event evidence or identity changed on restart')
            phases[name]=g
        return dict(exact_native_and_event_replay=True,phases=phases)
    check('LT06_IO01_IO06_restart_phases',restart)

    def live_amr_restart():
        options={'problem/flow':'shock_tube','problem/user_srcs':'false',
            'tracer_seed_initial/slab_min':.38,'tracer_seed_initial/slab_max':.42,
            'problem/refine_time':.02,'problem/coarsen_time':.06,
            'meshblock/nx1':64,'meshblock/nx2':4,
            'mesh_refinement/refinement':'adaptive','mesh_refinement/num_levels':2,
            'mesh_refinement/refinement_interval':1,'mesh_refinement/ncycle_check':1,
            'mesh_refinement/max_nmb_per_rank':128,'amr_criterion0/method':'user'}
        fixture='inputs/tests/shocks-curves/planar.athinput'
        base=run('live_amr',options,prescribed=True,fixture=fixture)
        expected=read_history(base/'shock_events');suspended=read_suspended(base/'shock_events')
        ref,clocks=thermo(base);observed=traces(base/'shock_events')
        states=snapshots(base);counts={g:len(s['amr']['leaves']) for g,s in states.items()}
        require(max(counts.values())>counts[min(counts)] and
                counts[max(counts)]<max(counts.values()),'no actual refine/coarsen')
        active=[g for (g,p),o in observed.items() if o['phase']==1 and
                counts[g]==max(counts.values())]
        require(active,'no active particle encounter on refined mesh')
        generation=min(active);options=dict(options, **{'time/nlim':generation,'output1/dt':1})
        seed=run('live_amr_seed',options,prescribed=True,fixture=fixture)
        cp=sorted((seed/'rst').glob('*.rst'))[-1]
        resumed=run('live_amr_resume',prescribed=True,restart=cp,
                    overrides=('time/nlim=-1','job/basename=amr_continued'))
        require(thermo(resumed)[0]=={k:v for k,v in ref.items() if k[0]>generation},
                'AMR native particle continuation changed')
        require(canonical(read_history(seed/'shock_events'))==canonical(expected),
                'AMR particle event evidence/identity changed')
        require(canonical(read_suspended(seed/'shock_events'))==canonical(suspended),
                'AMR active terminal state changed')
        require(canonical([o for _,o in sorted(traces(seed/'shock_events').items())])==
                canonical([o for _,o in sorted(observed.items())]),'AMR complete traces changed')
        # Callback configuration must remain bound even without the curve consumer.
        run('reject_amr_callback',prescribed=True,restart=cp,
            overrides=('problem/refine_time=.03','job/basename=reject_amr_callback'),reject='contract')
        return dict(events=len(expected),suspended=len(suspended),checkpoint_cycle=generation,
                    exact_particle_events_state_and_samples=True,leaf_counts=counts)
    check('LT04_LT06_live_AMR_active_restart',live_amr_restart)


    def transport_rejection():
        comparison=run('mc_comparison',dict(uniform, **{
            'particles/particle_type':'lagrangian_mc','particles/pusher':'lagrangian_mc'}))
        samples,clocks=thermo(comparison)
        moments=moment_oracle(samples,clocks)
        require(moments['mean_pass'] and moments['variance_pass'],
                'native MC mass-flux transport regressed')
        mc=run('mc_native',dict(uniform, **{'particles/particle_type':'lagrangian_mc',
                 'particles/pusher':'lagrangian_mc','time/nlim':2,'output1/dt':1}))
        cp=sorted((mc/'rst').glob('*.rst'))[-1]
        run('reject_transport',restart=cp,overrides=('particles/particle_type=lagrangian_tracer',
            'particles/pusher=lagrangian_tracer'),reject='transport does not match')
        run('reject_mc_history',{'particles/particle_type':'lagrangian_mc',
            'particles/pusher':'lagrangian_mc'},reject='Shock tracking requires')
        run('reject_mixed_pusher',{'particles/pusher':'lagrangian_mc'},reject='must match')
        return dict(mc_native_tracking_off_works=True,mc_histories_rejected=True,
                    restart_transport_change_rejected=True,mc_moments=moments)
    check('LT08_CFG01_MC_excluded',transport_rejection)

    def births():
        directory=run('delayed_inside',{'tracer_seed_initial/start_time':.1,
            'tracer_seed_initial/end_time':.1,'tracer_seed_initial/slab_min':.34,
            'tracer_seed_initial/slab_max':.36})
        observed=traces(directory/'shock_events');first={}
        for (g,p),o in sorted(observed.items()):first.setdefault(p,o)
        require(len(first)==128,'scheduled birth count changed')
        require(all(o['time']>=.1 for o in first.values()),'pre-insertion observation')
        require(all(o['sequence']==0 for o in first.values()),'newborn inherited sequence')
        inside=[o for o in first.values() if o['inside']]
        require(inside and all(o['flags']&1 for o in inside),'inside newborn lacks censor flag')
        return dict(actual_birth_times=sorted({o['time'] for o in first.values()}),
                    newborns=len(first),born_inside=len(inside),fresh_sequence_zero=True)
    check('LT07_EVT02_delayed_birth',births)

    def zero():
        directory=run('zero',{'tracer_seed_initial/count_per_event':0})
        require(not read_history(directory/'shock_events'),'zero particles emitted encounters')
        require(not thermo(directory)[0],'zero particles emitted native observations')
        return dict(events=0)
    check('LT09_CFG06_empty',zero)

    def bounds_and_corruption():
        run('overflow',{'shock_tracking/event_capacity':1,
            'tracer_seed_initial/slab_min':.45,'tracer_seed_initial/slab_max':.46},
            reject='event buffer overflow')
        seed=run('corruption_seed',{'time/nlim':2,'output1/dt':1})
        cp=sorted((seed/'rst').glob('*.rst'))[-1]
        chunk=seed/'shock_events/cycle0.rank0.csv'
        chunk.write_text(chunk.read_text()+'corruption\n')
        run('reject_corrupt_history',restart=cp,
            overrides=('time/nlim=-1','job/basename=corrupt'),reject='corrupted')
        return dict(overflow_rejected=True,history_corruption_rejected=True)
    check('LT10_IO08_bounded_output_integrity',bounds_and_corruption)

    def predictor_exit():
        # The accepted acceleration fixture supplies a spatially uniform reversal.
        # Birth x=0.99609375; predictor exits, corrected endpoint remains inside.
        options=dict(smooth, **{'material_curves/enabled':'false',
            'mesh/x1min':0,'mesh/x1max':1,'mesh/nx1':128,'meshblock/nx1':128,
            'mesh/ix1_bc':'outflow','mesh/ox1_bc':'outflow','problem/v1':1,
            'problem/v2':0,'problem/flow':'acceleration','problem/a1':-2000,
            'problem/signal_speed':1.1,'time/cfl_number':.9,'time/nlim':1,'time/tlim':.01,
            'tracer_seed_initial/x1min':.993,'tracer_seed_initial/x1max':1,
            'tracer_seed_initial/x2min':-.1,'tracer_seed_initial/x2max':.1})
        run('reject_predictor_exit',options,prescribed=True,
            fixture='inputs/tests/shocks-curves/smooth.athinput',
            reject='RK2 predictor left physical domain')
        return dict(predictor_exit_fails_closed=True,no_fabricated_exit_event=True)
    check('LT11_predictor_outflow_guard',predictor_exit)

    report['status']='pass' if all(c['status']=='pass' for c in report['checks']) else 'fail'
    save();print(json.dumps(dict(status=report['status'], report=str(work/'report.json'))))
    return report['status']!='pass'


if __name__=='__main__':
    raise SystemExit(main())
