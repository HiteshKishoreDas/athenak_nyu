"""Bounded independent curve-history acceptance; retains every input/log/failure.

Scientific strength/timing calibration and hardware qualification remain pending.
"""
import argparse
from collections import Counter
import importlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import struct
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'), str(ROOT/'tst/shocks-mc')]
from common import edit_input, provenance, sha, thermo
import material_curves


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(events):
    return json.dumps([{k:v for k,v in event.items() if k != 'rank'}
                       for event in events], sort_keys=True, allow_nan=True)


def snapshots(directory):
    result = {}
    for path in Path(directory).glob('curve_shock.curves.*.txt'):
        state = material_curves.read(path)
        # initial/final may repeat one generation exactly.
        result[state['cycle']] = state
    return result


def nodes(state):
    return {n['id']: n for curve in state['curves'] for n in curve['nodes']}


def trace_rows(directory, reader, birth=False):
    result = {}
    prefix = 'birthtrace' if birth else 'trace'
    for path in Path(directory).glob(prefix+'*.rank*.csv'):
        generation = int(re.fullmatch(prefix+r'(\d+)\.rank\d+\.csv', path.name)[1])
        for node_id, observation in reader.read_trace(path).items():
            key = generation, node_id
            require(key not in result and observation['generation'] == generation,
                    'duplicate or misaligned node trace identity')
            result[key] = observation
    return result


def density_fronts(directory):
    """Native central density midpoint roots; independent of detector/events."""
    result = []
    for path in sorted((Path(directory)/'tab').glob('*.tab')):
        text = path.read_text()
        match = re.search(r'time\s*=\s*([^\s]+)', text)
        require(match is not None, 'native profile has no time header')
        t = float(match[1]); rows = []
        for line in text.splitlines():
            if line and not line.startswith('#'):
                fields = line.split()
                rows.append((float(fields[-6]), float(fields[-5])))
        rows.sort(); roots = []
        midpoint = (1 + 8/3)/2
        for (x0,r0),(x1,r1) in zip(rows, rows[1:]):
            if r0 >= midpoint > r1 and .2 < (x0+x1)/2 < .5:
                roots.append(x0+(x1-x0)*(r0-midpoint)/(r0-r1))
        require(len(roots) == 1, 'native central density front not uniquely resolved')
        cycle = re.search(r'cycle\s*=\s*([0-9]+)', text)
        require(cycle is not None, 'native profile lacks cycle')
        result.append((int(cycle[1]), t, roots[0]))
    return result


def audit_native_leaf(directory, state, observation):
    """Derive leaf geometry from native logical tree, independently of sampler."""
    if not state['amr']['active']:
        return
    block = None; values = {}
    for line in (Path(directory)/'input.athinput').read_text().splitlines():
        match = re.fullmatch(r'\s*<([^>]+)>\s*',line)
        if match:block=match[1]
        assignment = re.match(r'\s*([^#=]+?)\s*=\s*([^#]+)',line)
        if block and assignment:values[block+'/'+assignment[1].strip()]=assignment[2].strip()
    dims=state['dimensions']; domain=state['domain']
    perblock=[int(values['meshblock/nx'+str(d+1)]) for d in range(3)]
    roots=[dims[d]//perblock[d] for d in range(3)]
    rootlevel=(max(roots)-1).bit_length()
    x=[domain[2*d]+(observation['x'][d]-domain[2*d])%
       (domain[2*d+1]-domain[2*d]) for d in range(3)]
    matches=[]
    for leaf in state['amr']['leaves']:
        level=leaf['level']; lengths=[]; located=True
        for d in range(3):
            extent=domain[2*d+1]-domain[2*d]
            count=roots[d]*(2**(level-rootlevel) if dims[d]>1 else 1)
            length=extent/count; lower=domain[2*d]+leaf['lx'+str(d+1)]*length
            located = located and lower<=x[d]<lower+length
            lengths.append(length/perblock[d])
        if located:matches.append((level,lengths))
    require(len(matches)==1,'node does not map uniquely into native saved leaf tree')
    level,spacing=matches[0]
    require(observation['level']==level,'sample uses stale or wrong native leaf level')
    require(max(abs(a-b) for a,b in zip(observation['spacing'],spacing))<=2e-12,
            'sample uses stale or wrong native local cell spacing')


def analyze(directory, physical=False):
    reader = importlib.import_module('curve_shock_history')
    directory = Path(directory); history = directory/'shock_events/curves'
    events = reader.read_history(history); suspended = reader.read_suspended(history)
    states = snapshots(directory); require(states, 'missing native curve snapshots')
    observed = trace_rows(history, reader); births = trace_rows(history, reader, True)
    require(observed, 'enabled traces missing')
    comparisons = 0; maximum_coordinate_roundoff=0
    by_generation={}
    for generation,node_id in observed:
        by_generation.setdefault(generation,set()).add(node_id)
    require(set(by_generation)==set(states),'full trace/native snapshot generation coverage differs')
    for generation,state in states.items():
        require(by_generation[generation]==set(nodes(state)),
                'full trace/native snapshot exact node coverage differs')
    for (generation,node_id),o in observed.items():
        require(generation in states, 'node trace lacks native synchronized snapshot')
        state = states[generation]; n = nodes(state).get(node_id)
        require(n is not None, 'traced node absent from synchronized snapshot')
        require(o['time'] == state['time'] and o['birth'] == n['birth'],
                'node sample time/generation/birth mismatch')
        lifted=[o['x'][d]+o['image'][d]*(state['domain'][2*d+1]-state['domain'][2*d])
                for d in range(3)]
        coordinate_roundoff=max(abs(a-b) for a,b in zip(lifted,n['x']))
        maximum_coordinate_roundoff=max(maximum_coordinate_roundoff,coordinate_roundoff)
        require(coordinate_roundoff<=2e-12,
                'node wrapped sample/lifted snapshot position mismatch')
        require(o['time'] >= n['birth'], 'invented prebirth observation')
        audit_native_leaf(directory,state,o)
        comparisons += 1
    final_nodes = nodes(states[max(states)])
    for e in events+suspended:
        n = final_nodes[e['node_id']]
        require(e['node_id'] >> 24 == e['curve_id'], 'curve/node namespace mismatch')
        require(e['birth'] == n['birth'], 'event birth provenance mismatch')
        for name in ('pre','first','last','post','representative'):
            if e[name]['valid']:
                require(e[name]['time'] >= n['birth'], 'event inherited prebirth observation')
    earliest={}
    for (generation,node_id),o in list(observed.items())+list(births.items()):
        if node_id not in earliest or o['time']<earliest[node_id]['time']:
            earliest[node_id]=o
    require(set(earliest)==set(final_nodes),'node lacks initial/birth history observation')
    for node_id,node in final_nodes.items():
        first=earliest[node_id]
        require(first['time']==node['birth'] and first['sequence']==0,
                'earliest node observation is not actual birth/sequencezero')
        if first['inside']:
            require(first['flags']&1,'earliest inside birth not leftcensored')
        else:
            require(first['phase']==0,'earliest outside birth inherited active state')
    # Completed-step adaptation is observed in the full trace; stage-one
    # adaptation is observed in birthtrace. Classify actual earliest births.
    positive_births = {node_id:o for node_id,o in earliest.items() if o['birth'] > 0}
    inside_births = {key:o for key,o in positive_births.items() if o['inside']}
    outside_births = {key:o for key,o in positive_births.items() if not o['inside']}
    for node_id,o in positive_births.items():
        require(o['time'] == o['birth'] and o['sequence'] == 0,
                'newborn history did not start at birth/sequence zero')
        if o['inside']:
            require(o['flags'] & 1, 'inside newborn missing left censor')
        else:
            require(o['phase'] == 0, 'outside newborn inherited active encounter')
    counts = Counter(e['node_id'] for e in events)
    metrics = {'events':len(events), 'completed':sum(e['status']==0 for e in events),
               'suspended':len(suspended), 'nodes':len(final_nodes),
               'synchronized_observations':comparisons,
               'fulltrace_exact_native_node_sets':True,
               'earliest_birth_observations_checked':len(earliest),
               'maximum_lift_roundoff':maximum_coordinate_roundoff,
               'positive_births':len(positive_births),
               'inside_positive_births':len(inside_births),
               'outside_positive_births':len(outside_births),
               'extra_episodes':{str(k):v for k,v in counts.items() if v>1},
               'scientific_approval':'pending'}
    if physical:
        require(events, 'physical shock produced no completed event histories')
        initial = nodes(states[min(states)])
        for node_id in initial:
            require(counts[node_id] == 1, 'initial node lacks exactly one completed passage')
        require(all(e['status']==0 for e in events), 'physical passage unexpectedly partial')
        fronts = density_fronts(directory); by_cycle = {g:x for g,t,x in fronts}
        require(len(by_cycle) >= 2, 'insufficient native density profile evidence')
        density_brackets = {}; covered=0
        for node_id in initial:
            samples = []
            for generation,state in sorted(states.items()):
                if generation in by_cycle and node_id in nodes(state):
                    samples.append((state['time'], nodes(state)[node_id]['x'][0]
                                    - by_cycle[generation]))
            brackets = [[a[0],b[0]] for a,b in zip(samples,samples[1:]) if a[1]>0>=b[1]]
            require(len(brackets)==1, 'native density crossing not uniquely bracketed')
            density_brackets[str(node_id)] = brackets[0]
            event=next(e for e in events if e['node_id']==node_id)
            lower=event['pre']['time'] if event['pre']['valid'] else event['first']['time']
            upper=event['post']['time'] if event['post']['valid'] else event['last']['time']
            require(lower<=brackets[0][1] and upper>=brackets[0][0],
                    'event evidence envelope excludes native density-front crossing bracket')
            covered+=1
        metrics.update(native_density_fronts=fronts, density_crossing_brackets=density_brackets,
            density_brackets_covered_by_event_envelope=covered,
            pressure_mach_errors=[e['representative']['mach_pressure']-2 for e in events],
            region_entry_brackets=[[e['pre']['time'],e['first']['time']] for e in events],
            region_exit_brackets=[[e['last']['time'],e['post']['time']] for e in events])
    return metrics, events, suspended, states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--prescribed-executable', type=Path)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--suite', choices=('core','forced'),default='core')
    args = parser.parse_args()
    require(args.executable or args.baseline, 'executable or baseline required')
    work = args.workdir.resolve(); work.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {'schema':'curve-shocks-acceptance-1', 'runs':[], 'checks':[],
              'scientific_approval':'pending', 'hardware_scope':'laptop serial only',
              'criteria_sha256':sha(ROOT/'docs/curve-shocks-oracles.md'),
              'runner_sha256':sha(Path(__file__)),
              'inspected_worktree_source_sha256_not_binary_binding':{str(p.relative_to(ROOT)):sha(p) for p in
                  sorted((ROOT/'src/material_curves').glob('*.cpp'))+
                  sorted((ROOT/'src/shocks').glob('*.hpp'))}}
    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        (work/'report.json').write_text(json.dumps(report,indent=2,allow_nan=True)+'\n')
        lines=['# Curve shock acceptance run','',
               'Frozen criteria: `'+report['criteria_sha256']+'`. Scientific approval pending.','',
               'Purpose: independent sampled front/path evidence; no exact crossing estimator.','']
        for item in report['runs']:
            lines += ['- Native `'+item['name']+'`: exit `'+str(item['exit'])+'`, '+
                      str(round(item['elapsed_seconds'],3))+' seconds. Input/log hashes in report.',
                      '  Command: `'+ ' '.join(item['command'])+'`.']
        lines += ['', '## Checks','']
        for row in report['checks']:
            lines.append('- '+row['id']+': '+row['status']+'. '+row.get('error',''))
        (work/'steps.md').write_text('\n'.join(lines)+'\n')
    def run(name, binary, fixture='planar', changes=None, restart=None, overrides=(), reject=None):
        require(time.monotonic()-started < 300, '300 second suite budget exhausted')
        directory = work/name; directory.mkdir()
        actual = directory/'input.athinput'
        options = {'job/basename':'curve_shock',
                   'shock_tracking/directory':str(directory/'shock_events')}
        options.update(changes or {})
        actual.write_text(edit_input((ROOT/'inputs/tests/shocks-curves'/
                                     (fixture+'.athinput')).read_text(), options))
        command = ([str(binary.resolve()), '-i', str(actual)] if restart is None else
                   [str(binary.resolve()), '-r',str(restart)])
        command += list(overrides)
        begin = time.monotonic()
        with (directory/'run.log').open('w') as log:
            try:
                code = subprocess.run(command,cwd=directory,stdout=log,
                    stderr=subprocess.STDOUT,timeout=30).returncode
            except subprocess.TimeoutExpired:
                code = 'timeout'
        item = dict(name=name,command=command,exit=code,
                    elapsed_seconds=time.monotonic()-begin,input_sha256=sha(actual),
                    log_sha256=sha(directory/'run.log'),**provenance(binary))
        if restart is not None:
            item['restart_sha256']=sha(restart)
            item['generated_input_used']=False
        report['runs'].append(item); save()
        if reject is not None:
            log = (directory/'run.log').read_text()
            require(isinstance(code,int) and code != 0 and reject.lower() in log.lower(),
                    name+' missing explicit rejection '+reject+'; see '+str(directory/'run.log'))
        else:
            require(code == 0, name+' failed native run; see '+str(directory/'run.log'))
        return directory
    def check(name, action):
        row = {'id':name}
        try:
            row.update(action() or {}); row['status'] = 'pass'
        except Exception as error:
            row.update(status='fail',error=type(error).__name__+': '+str(error))
        report['checks'].append(row);save()
    if args.suite=='forced':
        require(args.executable is not None,'forced suite needs candidate executable')
        def forced():
            # Reuse native binary parser, not input-header byte comparisons.
            import numpy as np
            from bin_convert import read_binary
            reader=importlib.import_module('curve_shock_history')
            from shock_history import read_history as read_mc
            on=run('forced_on',args.executable,'forced',{'shock_tracking/trace':'true'})
            off=run('forced_off',args.executable,'forced',{'shock_tracking/enabled':'false'})
            no=run('forced_no_trace',args.executable,'forced')
            metrics,events,suspended,states=analyze(on)
            def volumes(directory):
                result={}
                for path in sorted((directory/'bin').glob('*.bin')):
                    data=read_binary(path);key=data['cycle']
                    result[key]=(data['time'],data['mb_data'],data['mb_geometry'])
                require(result,'forced fixture lacks native fluid snapshots')
                return result
            reference=volumes(on);reference_mc=thermo(on)
            reference_nodes=material_curves.read(on/'curve_shock.curves.final.txt')
            for other in (off,no):
                actual=volumes(other)
                require(reference.keys()==actual.keys(),'forced output generations changed')
                for generation,(t,fields,geometry) in reference.items():
                    at,af,ag=actual[generation]
                    require(t==at and np.array_equal(geometry,ag) and fields.keys()==af.keys(),
                            'forced fluid metadata/geometry changed')
                    require(all(np.array_equal(fields[k],af[k]) for k in fields),
                            'forcing/tracking changed native fluid arrays')
                require(thermo(other)==reference_mc,'forced tracking changed native continuous-tracer trajectories')
                require(json.dumps(material_curves.read(other/'curve_shock.curves.final.txt'),
                                   sort_keys=True)==json.dumps(reference_nodes,sort_keys=True),
                        'forced tracking changed native material geometry')
            require(canonical(reader.read_history(no/'shock_events/curves'))==canonical(events),
                    'forced trace setting changed curve events')
            require(canonical(read_mc(no/'shock_events'))==canonical(read_mc(on/'shock_events')),
                    'forced trace setting changed MC events')
            require(not list((no/'shock_events/curves').glob('*trace*.csv')),
                    'forced optional trace disabled but files emitted')
            require(not (off/'shock_events').exists(),'forced global disabled writes histories')
            metrics.update(exact_fluid_tracer_curve_passivity=True,trace_optional_events_equal=True,
                physical_settings='inputs/tests/shocks-curves/forced.athinput',
                physical_strength_truth='not asserted',
                inherited_RK2_forcing_heat_bias='unresolved; no scientific approval')
            return metrics
        check('CSH09_forced_tracer_curves_passivity',forced)
        report['status']='pass' if report['checks'][0]['status']=='pass' else 'failed'
        save();print(json.dumps({'status':report['status'],'report':str(work/'report.json')}))
        return 0 if report['status']=='pass' else 1
    if args.baseline:
        def red():
            directory = run('old_coexist_rejected',args.baseline)
            paths = list((directory/'shock_events/curves').glob('manifest.rank*'))
            require(paths, 'red baseline: native run has no per-node shock histories')
        check('RED_native_coexistence',red)
        def curve_only_red():
            directory = run('old_missing_curve_histories',args.baseline,'planar_curveonly')
            require(list((directory/'shock_events/curves').glob('manifest.rank*')),
                    'red baseline: native curve-only run has no per-node histories')
        check('RED_native_missing_feature',curve_only_red)
        def native_control():
            directory=run('old_curve_native_control',args.baseline,'planar_curveonly',
                          {'shock_tracking/enabled':'false'})
            states=snapshots(directory);require(states,'native curve snapshots absent')
            fronts=density_fronts(directory)
            return {'native_fluid_curve_run':True,'steps':max(states),
                    'nodes':len(nodes(states[max(states)])),
                    'central_density_fronts':fronts,
                    'scope':'native curve/fluid control only; history absent'}
        check('ORACLE_native_density_control',native_control)
        report['status']='red_baseline_recorded'
        save(); print(json.dumps({'status':report['status'],'report':str(work/'report.json')}))
        return 0
    def planar():
        directory=run('planar',args.executable)
        metrics,*_=analyze(directory,physical=True)
        return metrics
    check('CSH02_CSH03_planar',planar)
    def curve_only():
        directory=run('curve_only',args.executable,'planar_curveonly')
        metrics,*_=analyze(directory,physical=True)
        return metrics
    check('CSH02_curve_only_global_enable',curve_only)
    def inside():
        directory=run('initial_inside',args.executable,changes={
            'curve1/x1':.3,'curve1/end1':.3,'time/tlim':.04})
        metrics,events,suspended,_=analyze(directory)
        require(events+suspended and all(e['flags']&1 for e in events+suspended),
                'inside initial node lacks left-censored episode')
        return metrics
    check('CSH04_initial_inside',inside)
    def passivity():
        on=work/'planar';require((on/'curve_shock.curves.final.txt').exists(), 'planar missing')
        off=run('tracking_off',args.executable,changes={'shock_tracking/enabled':'false','output1/dt':1})
        notrace=run('trace_off',args.executable,changes={'shock_tracking/trace':'false'})
        for other in (off,notrace):
            require(json.dumps(material_curves.read(on/'curve_shock.curves.final.txt'),sort_keys=True)==
                    json.dumps(material_curves.read(other/'curve_shock.curves.final.txt'),sort_keys=True),
                    'tracking changed native curve snapshot')
            require(thermo(on)==thermo(other), 'tracking changed native continuous-tracer trajectories')
            tables=lambda d:[p.read_bytes() for p in sorted((d/'tab').glob('*.tab'))]
            require(tables(on)==tables(other), 'tracking changed fluid conserved tables')
        reader=importlib.import_module('curve_shock_history')
        require(canonical(reader.read_history(on/'shock_events/curves'))==
                canonical(reader.read_history(notrace/'shock_events/curves')),
                'optional trace changed canonical events')
        require(not list((notrace/'shock_events/curves').glob('*trace*.csv')),
                'disabled optional traces still written')
        require(not (off/'shock_events').exists(), 'tracking disabled still creates events')
        return {'exact_fluid_tracer_geometry_equal':True,'trace_optional_events_equal':True}
    check('CSH05_coexist_passivity',passivity)
    def restart_phases():
        reader=importlib.import_module('curve_shock_history')
        base=work/'planar'; observed=trace_rows(base/'shock_events/curves',reader)
        expected=reader.read_history(base/'shock_events/curves')
        expected_suspended=reader.read_suspended(base/'shock_events/curves')
        encode_trace=lambda rows:json.dumps([[g,p,o] for (g,p),o in sorted(rows.items())],
                                           sort_keys=True,allow_nan=True)
        base_thermo,_=thermo(base); phases={}
        for name,phase in [('outside',0),('active',1),('pending',2)]:
            generations=[g for (g,p),o in observed.items() if g>0 and o['phase']==phase]
            require(generations,'phase not exercised: '+name)
            generation=min(generations);phases[name]=generation
            seed=run('seed_'+name,args.executable,changes={'time/nlim':generation,'output1/dt':1})
            checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
            # Keep prefix directory embedded in checkpoint; new basename avoids
            # overwriting the seed's immutable checkpoint pair.
            resumed=run('resume_'+name,args.executable,restart=checkpoint,
                overrides=('time/nlim=-1','job/basename=curve_resumed_'+name))
            require(canonical(reader.read_history(seed/'shock_events/curves'))==canonical(expected),
                    'continued event identity/evidence differs in '+name)
            actual,_=thermo(resumed)
            expected_mc={k:v for k,v in base_thermo.items() if k[0]>generation}
            require(expected_mc and actual==expected_mc,
                    'native MC continuation missing/excess generations, IDs or changed values in '+name)
            require(canonical(reader.read_suspended(seed/'shock_events/curves'))==
                    canonical(expected_suspended),'continued terminal censored view differs')
            require(encode_trace(trace_rows(seed/'shock_events/curves',reader))==
                    encode_trace(observed),'complete continued curve traces differ')
            final=material_curves.read(resumed/('curve_resumed_'+name+'.curves.final.txt'))
            reference=material_curves.read(base/'curve_shock.curves.final.txt')
            require(json.dumps(nodes(final),sort_keys=True)==json.dumps(nodes(reference),sort_keys=True),
                    'continued material node state differs')
        return {'phase_generations':phases,'samepartition_events_nodes_tracer_exact':True}
    check('CSH06_restart_actual_phases',restart_phases)
    def relative_directory_restart():
        reader=importlib.import_module('curve_shock_history')
        seed=run('relative_seed',args.executable,changes={'time/nlim':1,
            'output1/dt':1,'shock_tracking/directory':'shock_events'})
        cp=sorted((seed/'rst').glob('*.rst'))[-1]
        run('relative_launch',args.executable,restart=cp,
            overrides=('-d',str(seed.resolve()),'time/nlim=-1',
                       'job/basename=curve_relative_resumed'))
        require(canonical(reader.read_history(seed/'shock_events/curves'))==
                canonical(reader.read_history(work/'planar/shock_events/curves')),
                '-d relative historyprefix continuation differs')
        require((seed/'curve_relative_resumed.curves.final.txt').is_file(),
                '-d did not place native resumed output in requested directory')
        return {'relative_historyprefix_after_chdir_replayed':True,
                'checkpoint_argument':'absolute; configured history prefix relative'}
    check('CSH06_d_relative_prefix',relative_directory_restart)

    def sidecar_rejections():
        seed=work/'seed_active'; cp=sorted((seed/'rst').glob('*.rst'))[-1]
        original=Path(str(cp)+'.curves.shock.rank0')
        require(original.is_file(),'curve shock sidecar absent')
        for name in ('missing','corrupt','future','wide_id'):
            copied=work/('artifacts_'+name);shutil.copytree(seed/'rst',copied)
            checkpoint=copied/cp.name;sidecar=Path(str(checkpoint)+'.curves.shock.rank0')
            if name=='missing':sidecar.unlink()
            elif name=='corrupt':
                raw=bytearray(sidecar.read_bytes());raw[-1]^=1;sidecar.write_bytes(raw)
            elif name=='wide_id':
                raw=bytearray(sidecar.read_bytes())
                header=struct.Struct('@8s5I2i8Qd')  # current native-ABI double header
                h=header.unpack_from(raw)
                require(h[0]==b'AKCSHCK1' and h[3]==8 and h[5]==4 and h[9]>0 and
                        header.size+h[10]+h[9]*h[4]+8==len(raw),'unexpected sidecar ABI')
                offset=header.size+h[10]
                node=struct.unpack_from('=Q',raw,offset)[0]
                require(node<1<<32,'source node not compact')
                struct.pack_into('=Q',raw,offset,node+(1<<32))
                digest=14695981039346656037
                for byte in raw[:-8]:digest=((digest^byte)*1099511628211)&((1<<64)-1)
                struct.pack_into('=Q',raw,len(raw)-8,digest)
                sidecar.write_bytes(raw)
            else:
                # Frozen native header: AKCSHCK1 then uint32 version at offset8.
                raw=bytearray(sidecar.read_bytes())
                require(raw[:8]==b'AKCSHCK1','unexpected curve checkpoint magic')
                struct.pack_into('=I',raw,8,2);sidecar.write_bytes(raw)
            run('reject_'+name,args.executable,restart=checkpoint,
                overrides=('time/nlim=-1','job/basename=reject_'+name),
                reject='curve shock node ID exceeds 32 bits' if name=='wide_id' else 'shock')
        copied=work/'artifacts_prefix';shutil.copytree(seed,copied)
        checkpoint=copied/'rst'/cp.name
        chunk=copied/'shock_events/curves/cycle0.rank0.csv'
        chunk.write_bytes(chunk.read_bytes()+b'# corrupt prefix\n')
        run('reject_prefix',args.executable,restart=checkpoint,
            overrides=('time/nlim=-1','job/basename=reject_prefix',
                       'shock_tracking/directory='+str(copied/'shock_events')),
            reject='prefix')
        copied=work/'artifacts_future_chunk';shutil.copytree(seed,copied)
        checkpoint=copied/'rst'/cp.name
        # The native checkpoint cycle is the last committed prefix generation.
        generation=max(snapshots(seed))
        chunk=copied/'shock_events/curves'/('cycle'+str(generation+1)+'.rank0.csv')
        chunk.write_bytes(b'x'*(8*1024*1024))
        run('reject_oversized_future_chunk',args.executable,restart=checkpoint,
            overrides=('time/nlim=-1','job/basename=reject_future_chunk',
                       'shock_tracking/directory='+str(copied/'shock_events')),
            reject='existing file differs')
        return {'missing_and_corrupt_curve_sidecar_rejected':True,
                'future_version_explicitly_rejected':True,'wide_ID_alias_with_repaired_checksum_rejected':True,
                'corrupt_curve_prefix_rejected':True,
                'oversized_uncommitted_future_chunk_rejected':True}
    check('CSH06_sidecar_integrity',sidecar_rejections)
    def overflow():
        run('overflow',args.executable,'planar_curveonly',
            changes={'shock_tracking/event_capacity':1,'curve1/x1':.4,'curve1/end1':.4},
            reject='overflow')
        return {'explicit_overflow_rejection':True}
    check('CSH08_overflow',overflow)
    def old_checkpoint():
        off=run('untracked_curve_seed',args.executable,'planar_curveonly',
                {'shock_tracking/enabled':'false','output1/dt':1,'time/nlim':1})
        cp=sorted((off/'rst').glob('*.rst'))[-1]
        run('reject_old_untracked_checkpoint',args.executable,restart=cp,
            overrides=('shock_tracking/enabled=true','time/nlim=-1',
                       'job/basename=reject_old_untracked_checkpoint'),reject='sidecar')
        return {'old_untracked_checkpoint_activation_rejected':True}
    check('CSH08_old_untracked_checkpoint',old_checkpoint)
    def disabled_consumer():
        cp=sorted((work/'seed_active/rst').glob('*.rst'))[-1]
        run('reject_disabled_curves',args.executable,restart=cp,
            overrides=('material_curves/enabled=false','time/nlim=-1',
                       'job/basename=reject_disabled_curves'),reject='curve shock history cannot be disabled')
        return {'stored_curve_history_cannot_drop_consumer':True}
    check('CSH06_disabled_curve_consumer',disabled_consumer)
    def changed_native_source():
        seed=run('native_source_seed',args.executable,'planar_curveonly',{
            'time/nlim':1,'output1/dt':1,'hydro/const_accel':'true',
            'hydro/const_accel_val':.1,'hydro/const_accel_dir':1})
        cp=sorted((seed/'rst').glob('*.rst'))[-1]
        run('reject_changed_native_source',args.executable,restart=cp,
            overrides=('hydro/const_accel_val=.2','time/nlim=-1',
                       'job/basename=reject_changed_native_source'),reject='contract')
        return {'allowed_native_hydro_settings_bound_to_checkpoint':True,
                'physics_qualification':'configuration integrity only'}
    check('CSH06_native_hydro_contract',changed_native_source)
    def nonfinite_forcing():
        for field in ('dt_update','tcorr'):
            run('reject_curveonly_nan_'+field,args.executable,'planar_curveonly',{
                'turb_driving/normalization':'accel_rms','turb_driving/accel_rms':10,
                'turb_driving/nlow':1,'turb_driving/nhigh':3,
                'turb_driving/'+field:'nan'},reject='finite native forcing parameters')
        return {'curve_only_nonfinite_forcing_rejected_before_evolution':True}
    check('CSH08_curveonly_forcing_input',nonfinite_forcing)
    if args.prescribed_executable:
        for name,flow,v1,v2 in [('zero','zero',0,0),('translation','translation',.2,-.1),
                               ('shear','shear',.2,-.1)]:
            def smooth(name=name,flow=flow,v1=v1,v2=v2):
                directory=run(name,args.prescribed_executable,'smooth',{
                    'problem/flow':flow,'problem/v1':v1,'problem/v2':v2})
                metrics,events,suspended,states=analyze(directory)
                require(not events and not suspended, 'smooth control emitted shock history')
                initial,final=states[min(states)],states[max(states)]
                by_id=nodes(final);error=0
                for node in nodes(initial).values():
                    t=final['time'];expected=[node['x'][0]+v1*t,node['x'][1]+v2*t,node['x'][2]]
                    if flow=='shear':expected[0]+=.2*(node['x'][1]*t+.5*v2*t*t)
                    error=max(error,max(abs(a-b) for a,b in zip(by_id[node['id']]['x'],expected)))
                require(error<=2e-12,'prescribed affine control exceeds existing roundoff budget')
                return dict(metrics,max_position_error=error,roundoff_budget=2e-12)
            check('CSH01_'+name,smooth)
        def amr():
            directory=run('amr_smooth',args.prescribed_executable,'amr')
            metrics,events,suspended,states=analyze(directory)
            require(not events and not suspended,'smooth AMR control emitted shock history')
            leaves=[len(state['amr']['leaves']) for _,state in sorted(states.items())]
            require(max(leaves)>leaves[0] and leaves[-1]<max(leaves),
                    'fixture did not actually refine and coarsen')
            require(metrics['positive_births']>0,'AMR fixture lacks positive-time insertion')
            return dict(metrics,leaf_counts=leaves)
        check('CSH07_actual_refine_coarsen',amr)
        def changed_callback_contract():
            seed=run('callback_seed',args.prescribed_executable,'amr',{'time/nlim':1})
            cp=sorted((seed/'rst').glob('*.rst'))[-1]
            run('reject_changed_refine_time',args.prescribed_executable,'amr',restart=cp,
                overrides=('problem/refine_time=.2','time/nlim=-1',
                           'job/basename=reject_changed_refine_time'),reject='contract')
            return {'authentic_checkpoint_changed_refinement_callback_rejected':True}
        check('CSH06_changed_callback_contract',changed_callback_contract)

        def live_amr(fixture='planar', name='amr_live_shock'):
            directory=run(name,args.prescribed_executable,fixture,{
                'problem/flow':'shock_tube','problem/user_srcs':'false',
                'problem/refine_time':.02,'problem/coarsen_time':.06,
                'meshblock/nx1':64,'meshblock/nx2':4,
                'mesh_refinement/refinement':'adaptive',
                'mesh_refinement/num_levels':2,
                'mesh_refinement/refinement_interval':1,
                'mesh_refinement/ncycle_check':1,
                'mesh_refinement/max_nmb_per_rank':128,
                'amr_criterion0/method':'user',
                'material_curves/adaptive':'true','material_curves/spacing_cells':1,
                'material_curves/max_spacing':0,'material_curves/max_angle':3,
                'material_curves/chord_error':0,
                'curve1/x1':.28,'curve1/end1':.45,'curve1/x2':.03,
                'curve1/end2':.07,'curve1/count':16})
            metrics,events,suspended,states=analyze(directory)
            leaves=[len(state['amr']['leaves']) for _,state in sorted(states.items())]
            require(max(leaves)>leaves[0] and leaves[-1]<max(leaves),
                    'live shock did not actually refine and coarsen')
            require(metrics['inside_positive_births']>0 and metrics['outside_positive_births']>0,
                    'live shock lacks both positive-time inside/outside births')
            return dict(metrics,leaf_counts=leaves)
        check('CSH04_CSH07_curveonly_live_amr_births',
              lambda: live_amr('planar_curveonly','amr_curveonly_live_shock'))
        def amr_birth_restart():
            reader=importlib.import_module('curve_shock_history')
            base=work/'amr_curveonly_live_shock'
            observed=trace_rows(base/'shock_events/curves',reader)
            generation=min(g for (g,p),o in observed.items()
                           if o['birth']>0 and o['inside'] and o['time']==o['birth'])
            seed=run('amr_birth_seed',args.prescribed_executable,str(base/'input'),
                     {'time/nlim':generation,'output1/dt':1})
            cp=sorted((seed/'rst').glob('*.rst'))[-1]
            resumed=run('amr_birth_resumed',args.prescribed_executable,restart=cp,
                        overrides=('time/nlim=-1','job/basename=curve_amr_resumed'))
            def encoded(rows):
                return json.dumps([[g,p,o] for (g,p),o in sorted(rows.items())],
                                  sort_keys=True,allow_nan=True)
            for load in (reader.read_history,reader.read_suspended):
                require(canonical(load(base/'shock_events/curves'))==
                        canonical(load(seed/'shock_events/curves')),
                        'AMR newborn checkpoint continuation differs')
            for birth in (False,True):
                require(encoded(trace_rows(base/'shock_events/curves',reader,birth))==
                        encoded(trace_rows(seed/'shock_events/curves',reader,birth)),
                        'AMR full/birth traces differ after checkpoint')
            reference=material_curves.read(base/'curve_shock.curves.final.txt')
            final=material_curves.read(resumed/'curve_amr_resumed.curves.final.txt')
            require(json.dumps(nodes(reference),sort_keys=True)==
                    json.dumps(nodes(final),sort_keys=True),
                    'AMR newborn restart changed final native node state')
            return {'checkpoint_generation':generation,'active_newborn_exercised':True,
                    'events_censored_views_full_birth_traces_nodes_exact':True,
                    'scope':'curve-only; MC AMR qualification remains blocked'}
        check('CSH06_CSH07_amr_newborn_restart',amr_birth_restart)
        check('CSH04_CSH07_coexisting_live_amr_births',live_amr)
    else:
        report['checks'].append({'id':'CSH01_CSH07','status':'not_run',
                                 'reason':'prescribed executable not supplied'})
    report['status']='pass' if all(c['status']=='pass' for c in report['checks']) else 'incomplete_or_failed'
    save(); print(json.dumps({'status':report['status'],'report':str(work/'report.json')}))
    return 0 if report['status']=='pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
