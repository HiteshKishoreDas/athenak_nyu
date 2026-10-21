"""MC encounter engineering checks; stochastic physical recovery is measured."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re
import shutil
import sys
import time

from common import ROOT, check_uniform_transfers, provenance, run, sha, thermo, uniform_changes

sys.path.insert(0, str(ROOT/'vis/python'))
from shock_history import parse_observation, read_history, read_suspended


def traces(directory):
    result = {}
    for path in Path(directory).glob('trace*.rank*.csv'):
        generation, rank = map(int, re.fullmatch(r'trace(\d+)\.rank(\d+)\.csv', path.name).groups())
        for row in csv.reader(path.read_text().splitlines()):
            if len(row) != 28:
                raise ValueError('trace schema mismatch')
            observation = parse_observation(row[4:])
            key = (generation,int(row[0]))
            if key in result or observation['generation'] != generation:
                raise ValueError('duplicate/misaligned trace generation')
            result[key] = dict(observation, sequence=int(row[1]), phase=int(row[2]),
                               flags=int(row[3]), rank=rank)
    return result


def canonical(events):
    # JSON produces stable NaN spellings; schema already limits permitted NaNs.
    return json.dumps([{k:v for k,v in event.items() if k!='rank'} for event in events],
                      sort_keys=True, allow_nan=True)


def recovery(directory):
    events = read_history(Path(directory)/'shock_events')
    suspended = read_suspended(Path(directory)/'shock_events')
    observed = traces(Path(directory)/'shock_events')
    sampled = {pid for generation,pid in observed}
    inside = {pid for (_,pid),o in observed.items() if o['valid'] and o['inside']}
    emitted = {e['particle_id'] for e in events}
    counts = Counter(e['particle_id'] for e in events)
    return {'sampled_particles':len(sampled), 'raw_region_particles':len(inside),
            'emitted_particles':len(emitted), 'emitted_events':len(events),
            'completed_events':sum(e['status']==0 for e in events),
            'partial_events':sum(e['status']==1 for e in events),
            'terminal_suspended':len(suspended),
            'particles_with_multiple_episodes':dict((str(p),c) for p,c in counts.items() if c>1),
            'raw_membership_recovery_fraction':len(emitted & inside)/len(inside) if inside else None,
            'physical_acceptance':'stochastic recovery, censoring and reentries measured; owner thresholds pending',
            'event_evidence':[{'particle_id':e['particle_id'],'sequence':e['sequence'],
                'status':e['status'],'flags':e['flags'],
                'entry_region_bracket':[e['pre']['time'],e['first']['time']],
                'exit_region_bracket':[e['last']['time'],e['post']['time']],
                'pressure_mach':e['representative']['mach_pressure'],
                'temperature_mach':e['representative']['mach_temperature'],
                'normal':e['representative']['normal']} for e in events]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args = parser.parse_args()
    work = args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    report = args.report.resolve();report.parent.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    results = {'executable_sha256':sha(args.executable),'input_sha256':sha(ROOT/'inputs/tests/shock_mc.athinput'),
               'runs':[],'checks':[],'physical':{},'scope':'RK2 existing MC tracer, serial uniform mesh only',
               'scientific_signoff':'pending','statistical_path_approval':'pending'}
    results.update(provenance(args.executable))
    def save():
        results['elapsed_seconds'] = time.monotonic()-started
        report.write_text(json.dumps(results,indent=2,allow_nan=True)+'\n')
    def execute(label,changes=None,restart=None,overrides=(),expected_error=None):
        if time.monotonic()-started > 300:raise TimeoutError('suite300s budget exhausted')
        changes = dict(changes or {})
        if restart is None:changes['shock_tracking/directory']=str(work/label/'shock_events')
        item = run(args.executable,work/label,changes,restart,overrides)
        item.update(label=label,expected_error=expected_error)
        log = (work/label/'run.log').read_text()
        item['passed']=(item['exit']==0 if expected_error is None else
            isinstance(item['exit'],int) and item['exit']!=0 and expected_error in log)
        results['runs'].append(item);save()
        target=report.parent/'runs'/label;target.mkdir(parents=True,exist_ok=True)
        for name in ('input.athinput','run.log'):
            if (work/label/name).exists():shutil.copy2(work/label/name,target/name)
        (target/'manifest.json').write_text(json.dumps(item,indent=2)+'\n')
        if not item['passed']:raise AssertionError(label+': expected '+str(expected_error)+', got '+str(item['exit']))
        return work/label
    try:
        changes=uniform_changes(count=1024);changes.update({'output3/dt':1,'output3/data_format':' %24.17e'})
        off = execute('uniform_off',changes)
        changes = uniform_changes(count=1024);changes.update({'shock_tracking/enabled':'true',
            'output3/dt':1,'output3/data_format':' %24.17e'})
        on = execute('uniform_on',changes)
        ot,otime = thermo(off);nt,ntime = thermo(on)
        assert ot==nt and otime==ntime
        assert sha(sorted((off/'tab').glob('*.tab'))[-1])==sha(sorted((on/'tab').glob('*.tab'))[-1])
        assert not (off/'shock_events').exists()
        assert read_history(on/'shock_events')==[]
        results['checks'].append(check_uniform_transfers(on,nt,ntime))
        results['checks'].append({'test_ids':['CFG01','PHY08','TRC01'],'status':'passed',
            'comparison':'exact samepartition fixedseed thermo observations; no deterministic physical path'})
        zero=execute('zero',{'tracer_seed_initial/count_per_event':0,'time/tlim':.02})
        assert read_history(zero/'shock_events')==[]
        results['checks'].append({'test_id':'CFG06','status':'passed'})
        contact=execute('contact',{'problem/dl':2,'problem/dr':1,'problem/pl':1,
            'problem/pr':1,'problem/ul':.3,'problem/ur':.3,
            'tracer_seed_initial/slab_min':.29,'tracer_seed_initial/slab_max':.31})
        assert read_history(contact/'shock_events')==[]
        results['checks'].append({'test_id':'PHY09','status':'passed','contact_events':0})
        base=execute('planar',{'output3/dt':1,'output3/data_format':' %24.17e'});events=read_history(base/'shock_events')
        baseoff=execute('planar_off',{'shock_tracking/enabled':'false','output3/dt':1,
            'output3/data_format':' %24.17e'})
        assert thermo(baseoff)==thermo(base)
        assert sha(sorted((baseoff/'tab').glob('*.tab'))[-1])==sha(sorted((base/'tab').glob('*.tab'))[-1])
        results['checks'].append({'test_id':'SYS02','status':'passed',
            'comparison':'exact fixedseed planar particle records and byte-identical17digit conserved slice'})
        results['physical']['planar']=recovery(base)
        assert events, 'no planar encounters recovered'
        observed=traces(base/'shock_events')
        # Match synchronized shock observation to independent existing thermo output.
        bt,times=thermo(base)
        assert observed.keys()==bt.keys()
        assert all(o['time']==bt[key]['time'] and tuple(o['x'])==bt[key]['x']
                   and o['rho']==bt[key]['fields']['density'] for key,o in observed.items())
        results['checks'].append({'test_ids':['INT06','INT08'],'status':'passed',
            'comparison':'one observation per generation/ID, exact MC thermo timestamp/x/rho'})
        multi=execute('multiblock',{'meshblock/nx1':32})
        mt,mtimes=thermo(multi);final=max(mtimes)
        initial={p:o for (g,p),o in mt.items() if g==0}
        surviving={p:o for (g,p),o in mt.items() if g==final}
        assert initial.keys()==surviving.keys()
        moved=sum(initial[p]['gid']!=surviving[p]['gid'] for p in initial)
        assert moved>0
        results['checks'].append({'test_ids':['INT01','IO01'],'status':'passed','local_gid_changes':moved,
            'comparison':'all initial IDs survive; partition-dependent seed mapping not forced identical'})
        results['physical']['multiblock']=recovery(multi)
        inside=execute('initial_inside',{'tracer_seed_initial/slab_min':.292,
                                        'tracer_seed_initial/slab_max':.308})
        results['physical']['initial_inside']=recovery(inside)
        assert any(e['flags']&1 for e in read_history(inside/'shock_events')+read_suspended(inside/'shock_events'))
        delayed=execute('delayed_inside',{'tracer_seed_initial/start_time':.1,
            'tracer_seed_initial/end_time':.1,'tracer_seed_initial/slab_min':.34,
            'tracer_seed_initial/slab_max':.36})
        results['physical']['delayed_inside']=recovery(delayed)
        dt,dtimes=thermo(delayed)
        assert len({p for g,p in dt})==128
        assert min(o['time'] for o in dt.values())>=.1
        assert any(e['flags']&1 for e in read_history(delayed/'shock_events')+read_suspended(delayed/'shock_events'))
        results['checks'].append({'test_ids':['EVT02','EVT05'],'status':'passed','initial_and_delayed_left_censoring':True})
        phases={}
        for phase,number in [('outside',0),('active',1),('pending',2)]:
            generations=[g for (g,p),o in observed.items() if g>0 and o['phase']==number]
            assert generations, 'unpopulated phase '+phase
            generation=min(generations);phases[phase]=generation
            seed=execute('seed_'+phase,{'time/nlim':generation,'output1/dt':1})
            checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
            continued=execute('resume_'+phase,restart=checkpoint,
                overrides=('time/nlim=-1','output1/dt=-1','job/basename=continued_'+phase))
            assert canonical(read_history(seed/'shock_events'))==canonical(events)
            restored,restored_times=thermo(continued)
            assert all(bt[key]==o for key,o in restored.items())
            results['checks'].append({'test_ids':['RST01','RST07'],'status':'passed','phase':phase,
                'generation':generation,'exact_replayed_event_content':True,'exact_particle_observations':True,
                'checkpoint_sha256':sha(checkpoint)})
        results['checkpoint_generations']=phases
        # Corrupt copies preserve successful evidence and original checkpoint.
        seed=work/'seed_active';checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
        for label,action,reason in [('missing_sidecar','missing','sidecar'),
                                  ('corrupt_sidecar','corrupt','checksum'),
                                  ('bad_history','history','history corrupted')]:
            copy=work/(label+'_artifacts');shutil.copytree(seed,copy)
            cp=copy/'rst'/checkpoint.name;sidecar=Path(str(cp)+'.shock.rank0')
            if action=='missing':sidecar.unlink()
            elif action=='corrupt':
                raw=bytearray(sidecar.read_bytes());raw[-1]^=1;sidecar.write_bytes(raw)
            else:(copy/'shock_events/cycle0.rank0.csv').write_text('truncated\n')
            execute(label,restart=cp,overrides=('output1/dt=-1',
                    'shock_tracking/directory='+str(copy/'shock_events')),expected_error=reason)
        results['checks'].append({'test_ids':['RST03','RST04','IO08'],'status':'passed'})
        execute('overflow',{'shock_tracking/event_capacity':1,
            'tracer_seed_initial/count_per_event':4096,'tracer_seed_initial/slab_min':.59,
            'tracer_seed_initial/slab_max':.60},expected_error='event buffer overflow')
        results['checks'].append({'test_ids':['IO03','IO04'],'status':'passed','capacity':1})
        execute('unsupported_rk3',{'time/integrator':'rk3'},expected_error='RK2')
        results['checks'].append({'test_id':'CFG','status':'passed','unsupported_integrator_rejected':'rk3'})
        results['status']='passed_supported_engineering_checks'
    except Exception as exc:
        results['status']='failed';results['error']=repr(exc)
    save();print(json.dumps({'status':results['status'],'report':str(report)}))


if __name__=='__main__':
    main()
