"""Existing-output scheduling, exact IDs, AMR births, restart and reader integrity.

Reuses accepted CT/CSH physics fixtures. No detector or physical tolerance tuning.
Run directories must be new; failed runs and exact commands are retained.
"""
import argparse
import json
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import edit_input, provenance, sha
from tracer_history import read_run
import curve_shock_history
import shock_history
import material_curves


def require(condition, message):
    if not condition:
        raise ValueError(message)


def normalized(value):
    if isinstance(value,float) and math.isnan(value):
        return 'NaN'
    if isinstance(value,dict):
        return {k:normalized(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [normalized(v) for v in value]
    return value


def geometry_by_cycle(data):
    result = {}
    for state in data['geometry']:
        row = normalized(state['curves'])
        require(state['cycle'] not in result or result[state['cycle']]==row,
                'same-cycle geometry differs across output streams')
        result[state['cycle']] = row
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--prescribed-executable',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--reference-runs',type=Path,help='optional retained pre-removal text/CSV runs')
    args = parser.parse_args()
    work = args.workdir.resolve(); work.mkdir(parents=True,exist_ok=False)
    report = dict(schema='tracer-output-check-1',runs=[],checks=[],
                  hardware='local Serial CPU/double; no MPI/GPU runtime claim',
                  scientific_approval='pending',runner_sha256=sha(__file__))

    def save():
        (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')

    def check(name, action):
        row = dict(id=name)
        try:
            row.update(action() or {}); row['status']='pass'
        except Exception as exc:
            row.update(status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row); save()

    def run(name, fixture='planar', changes=None, restart=None, overrides=(), reject=None, curve_only=False, prescribed=False):
        directory = work/name; directory.mkdir()
        binary = (args.executable if fixture=='planar' and not prescribed else args.prescribed_executable).resolve()
        options = {'output2/file_type':'prtcl_thermo_history',
                   'output2/source':'both','output2/variables':'rho,p,t,s,v1,v2,v3',
                   'output2/dcycle':1,'output4/file_type':'curves','output4/id':'shape',
                   'output4/dcycle':4}
        if fixture!='planar':
            options.update({'shock_tracking/enabled':'false',
                'particles/particle_type':'lagrangian_tracer','particles/pusher':'lagrangian_tracer',
                'particles/random_seed':11,'tracer_seed_initial/id':1,
                'tracer_seed_initial/count_per_event':128,'tracer_seed_initial/start_time':0,
                'tracer_seed_initial/seed':903,'tracer_seed_initial/region':'box',
                'tracer_seed_initial/weight':'volume','tracer_seed_initial/x1min':-.4,
                'tracer_seed_initial/x1max':.4,'tracer_seed_initial/x2min':-.4,
                'tracer_seed_initial/x2max':.4})
        options.update(changes or {})
        input_path = directory/'input.athinput'
        text=edit_input((ROOT/f'inputs/tests/shocks-curves/{fixture}.athinput').read_text(),options)
        if curve_only:
            text=re.sub(r'(?ms)^<(?:particles|tracer_seed[^>]*)>.*?(?=^<|\Z)','',text)
            text=edit_input(text,{'output2/source':'curves'})
        input_path.write_text(text)
        command = [str(binary),'-r' if restart else '-i',str(restart or input_path),*overrides]
        started = time.monotonic()
        with (directory/'run.log').open('w') as log:
            status = subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=120).returncode
        report['runs'].append(dict(name=name,command=command,exit=status,seconds=time.monotonic()-started,
            input_sha256=sha(input_path),restart_sha256=sha(restart) if restart else None,**provenance(binary)))
        save()
        log = (directory/'run.log').read_text()
        require((status!=0 and reject.lower() in log.lower()) if reject else status==0,
                name+': unexpected native exit/diagnostic; see run.log')
        return directory

    data = {}
    def pair(fixture):
        legacy = run(fixture+'_control',fixture,changes={'output4/dcycle':0})
        native = run(fixture,fixture)
        a,b = read_run(legacy),read_run(native)
        data[fixture] = (a,b)
        require(not list(native.rglob('*.tbin')),'removed format emitted')
        if args.reference_runs:
            historical=read_run(args.reference_runs/(fixture+'_legacy'))
            for key in ('particles','curve_nodes','particle_times','curve_times',
                        'particle_events','curve_events','particle_suspended','curve_suspended'):
                require(normalized(historical[key])==normalized(b[key]),'historical '+key+' changed')
            require(geometry_by_cycle(historical)==geometry_by_cycle(b),'historical geometry changed')
            old=next((args.reference_runs/(fixture+'_legacy')/'prtcl_thermo_history').glob('*.thp'))
            current=next((native/'prtcl_thermo_history').glob('*.thp'))
            require(old.read_bytes()==current.read_bytes(),'historical THP bytes changed')
            if fixture=='planar':
                for path in (native/'tab').glob('*.tab'):
                    require(path.read_bytes()==(args.reference_runs/'planar_legacy/tab'/path.name).read_bytes(),
                            'historical fluid bytes changed')
        for key in ('particles','curve_nodes','particle_times','curve_times',
                    'particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(a[key])==normalized(b[key]),fixture+': '+key+' changed')
        require(geometry_by_cycle(a)==geometry_by_cycle(b),fixture+': geometry changed')
        thp_a=next((legacy/'prtcl_thermo_history').glob('*.thp'))
        thp_b=next((native/'prtcl_thermo_history').glob('*.thp'))
        require(thp_a.read_bytes()==thp_b.read_bytes(), 'scheduled geometry changed THP bytes')
        if fixture=='planar':
            tabs=sorted((legacy/'tab').glob('*.tab'))
            require(tabs and all(p.read_bytes()==(native/'tab'/p.name).read_bytes() for p in tabs),
                    'conserved fluid tables changed')
        return dict(particle_samples=len(b['particles']),curve_samples=len(b['curve_nodes']),
                    particle_events=len(b['particle_events']),curve_events=len(b['curve_events']),
                    sampled_cycles=len(b['curve_times']),geometry_cycles=len(geometry_by_cycle(b)),
                    thp_byte_equal=True,all_records_exact=True)

    check('TO01_IO03_IO04_planar_all_records_exact',lambda:pair('planar'))
    check('TO02_AMR04_PROP06_AMR_all_records_exact',lambda:pair('amr'))

    def physical_oracle():
        maximum=0
        for d in data['amr']:
            for sample in list(d['particles'].values())+list(d['curve_nodes'].values()):
                fields=sample['fields']; pressure=(4-.2)**2/(5/3)
                for name,truth in {'density':1.,'pressure':pressure,'temperature':pressure,
                                   'entropy':math.log(pressure),'v1':.2,'v2':-.1,'v3':0.}.items():
                    maximum=max(maximum,abs(fields[name]-truth))
        require(maximum<=2e-12, 'existing manufactured-field error budget exceeded')
        births={}
        for (cycle,node),sample in sorted(data['amr'][1]['curve_nodes'].items()):
            births.setdefault(node,sample)
        inserted=[s for s in births.values() if s['birth']>0]
        require(inserted and all(s['birth']==s['time'] for s in inserted), 'inserted-node histories missing or predate birth')
        tables=[material_curves.read(p) for p in (work/'amr').glob('*.curves.*.txt')]
        leaves=[len(t['amr']['leaves']) for t in tables]
        require(min(leaves)<max(leaves), 'no live AMR change exercised')
        return dict(maximum_thermo_error=maximum,budget=2e-12,inserted_nodes=len(inserted),leaf_count_range=[min(leaves),max(leaves)])
    check('TO03_independent_thermo_oracle_and_births',physical_oracle)

    def large_ids():
        directory=run('large_ids','smooth',changes={'curve1/id':4294967295,'time/nlim':4})
        d=read_run(directory)
        ids={node for cycle,node in d['curve_nodes']}
        require(min(ids)>2**63 and max(ids)<2**64,'large uint64 fixture not exercised')
        geometry_ids={n['id'] for s in d['geometry'] for c in s['curves'] for n in c['nodes']}
        require(ids==geometry_ids,'uint64 thermo/geometry identity mismatch')
        return dict(minimum_id=min(ids),maximum_id=max(ids),nodes=len(ids))
    check('TO04_IO04_exact_uint64_above_2pow63',large_ids)

    def restart():
        # Use an absolute journal location exactly as the existing shock-restart tests do.
        shared=work/'restart-journal'
        seed=run('restart_seed',changes={'time/nlim':20,'shock_tracking/directory':str(shared),
                                       'output1/dcycle':-1})
        checkpoint=next((seed/'rst').glob('*.rst'))
        resumed=run('restart_resume',restart=checkpoint,overrides=('time/nlim=-1',))
        a,b=read_run(seed,shock_directory=str(shared)),read_run(resumed,shock_directory=str(shared))
        reference=data['planar'][1]
        for key in ('particles','curve_nodes','particle_times','curve_times'):
            merged=dict(a[key]); merged.update(b[key])
            require(normalized(merged)==normalized(reference[key]),'restart changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(b[key])==normalized(reference[key]),'restart changed '+key)
        require(geometry_by_cycle(b)[max(b['curve_times'])]==geometry_by_cycle(reference)[max(reference['curve_times'])],
                'restart changed final geometry')
        before={str(p.relative_to(shared)):sha(p) for p in shared.rglob('*') if p.is_file()}
        zero=run('restart_zero',restart=checkpoint,overrides=('time/nlim=20',))
        after={str(p.relative_to(shared)):sha(p) for p in shared.rglob('*') if p.is_file()}
        require(before==after, 'zero-step replay changed durable event files')
        return dict(restart_cycle=20,exact_continuation=True,zero_step_exit=0,
                    particle_events=len(b['particle_events']),curve_events=len(b['curve_events']))
    check('TO05_IO06_IO07_RST02_restart_replay',restart)

    def aliases():
        folder=work/'duplicate'; shutil.copytree(work/'planar/shock_events',folder)
        source=next(folder.glob('cycle*.csv')); shutil.copyfile(source,folder/source.name.replace('cycle','cycle0',1))
        try:shock_history.read_history(folder)
        except ValueError:return dict(duplicate_generation_rejected=True)
        raise ValueError('accepted duplicate generation alias')
    check('TO07_IO08_duplicate_generation_rejected',aliases)

    def cadence():
        d=run('cadence','smooth',changes={'time/nlim':10,'material_curves/output_dcycle':0,
            'output2/dcycle':2,'output4/dcycle':3})
        a=read_run(d)
        require(set(a['curve_times'])=={0,2,4,6,8,10},'thermo cycle cadence changed')
        require(set(geometry_by_cycle(a))=={0,3,6,9,10},'geometry native scheduler cadence changed')
        return dict(thermo_cycles=sorted(a['curve_times']),geometry_cycles=sorted(geometry_by_cycle(a)))
    check('TO10_PROP06_independent_output_clocks',cadence)

    def manifest_checks():
        results={}
        for carrier in ('particle','curve'):
            source=work/'planar/shock_events'
            if carrier=='curve':source=source/'curves'
            chunk=source/'cycle0.rank0.csv'
            orphan=work/('orphan-'+carrier);where=orphan/'shock_events'
            if carrier=='curve':where=where/'curves'
            where.mkdir(parents=True);shutil.copyfile(chunk,where/chunk.name)
            try:read_run(orphan)
            except ValueError as exc:
                require('orphan' in str(exc),'wrong orphan rejection');results[carrier+'_orphan']=True
            else:raise ValueError('silently hid orphan event files')
        return results
    check('TO12_IO09_missing_manifest_regressions',manifest_checks)

    def traces():
        counts={}
        for carrier,module in (('particle',shock_history),('curve',curve_shock_history)):
            relative=Path('shock_events')/('curves' if carrier=='curve' else '')
            native=work/'planar'/relative
            count=0
            for path in native.glob('*trace*.csv'):
                legacy=work/'planar_control'/relative/path.with_suffix('.csv').name
                a,b=module.read_trace(legacy),module.read_trace(path)
                require(normalized(a)==normalized(b),carrier+' trace changed')
                count+=len(b)
            require(count>0,'trace fixture absent');counts[carrier+'_trace_samples']=count
        return counts
    check('TO15_IO03_raw_flags_phase_trace_parity',traces)

    def missing_curves():
        run('missing_curves',changes={'material_curves/enabled':'false','output2/source':'particles'},
            reject='file_type=curves requires')
        return dict(missing_curve_module_rejected=True)
    check('TO16_CFG_geometry_without_curves_rejected',missing_curves)

    def mc():
        options={'particles/particle_type':'lagrangian_mc','particles/pusher':'lagrangian_mc','time/nlim':4}
        b=read_run(run('mc','smooth',changes=options))
        require(b['particles'] and b['curve_nodes'] and not b['particle_events'],'MC thermo fixture missing')
        return dict(particle_samples=len(b['particles']),shock_tracking=False)
    check('TO17_fixed_mesh_MC_thermo_only',mc)

    def curve_only():
        b=read_run(run('curve_only',curve_only=True))
        require(not b['particles'] and not b['particle_events'] and b['curve_events'],
                'curve-only shock fixture not exercised')
        return dict(curve_events=len(b['curve_events']),particle_module=False)
    check('TO18_curve_only_shocks_thermo_geometry',curve_only)

    def empty():
        d=run('empty_particles',changes={'tracer_seed_initial/target':'rho',
                                         'tracer_seed_initial/target_min':100})
        data=read_run(d)
        require(not data['particles'] and data['curve_nodes'] and not data['particle_events'],
                'empty particle population not exercised')
        return dict(particle_samples=0,curve_samples=len(data['curve_nodes']),empty_tables=True,
                    note='single rank with zero particles; not MPI empty-rank qualification')
    check('TO19_IO02_empty_particle_population',empty)

    def shared_directory():
        d=run('journal_in_run_root',changes={'shock_tracking/directory':'.'})
        result=read_run(d,shock_directory='.')
        reference=data['planar'][1]
        for key in ('particles','curve_nodes','particle_events','curve_events'):
            require(normalized(result[key])==normalized(reference[key]),'shared directory changed '+key)
        require(geometry_by_cycle(result)==geometry_by_cycle(reference),
                'confused curve journal tables with scheduled geometry')
        return dict(journal_in_run_root=True,curve_events_and_geometry_share_directory=True)
    check('TO20_native_geometry_discovery_with_shared_journal_directory',shared_directory)

    def live_shock_birth_restart():
        # Exact existing CSH04/CSH07 setup; test serialization, not a new AMR model.
        options={'problem/flow':'shock_tube','problem/user_srcs':'false',
            'problem/refine_time':.02,'problem/coarsen_time':.06,
            'meshblock/nx1':64,'meshblock/nx2':4,'mesh_refinement/refinement':'adaptive',
            'mesh_refinement/num_levels':2,'mesh_refinement/refinement_interval':1,
            'mesh_refinement/ncycle_check':1,'mesh_refinement/max_nmb_per_rank':128,
            'amr_criterion0/method':'user','material_curves/adaptive':'true',
            'material_curves/spacing_cells':1,'material_curves/max_spacing':0,
            'material_curves/max_angle':3,'material_curves/chord_error':0,
            'curve1/x1':.28,'curve1/end1':.45,'curve1/x2':.03,'curve1/end2':.07,'curve1/count':16}
        a=read_run(run('amr_shock_control',changes=dict(options,**{'output4/dcycle':0}),prescribed=True))
        directory=run('amr_shock',changes=options,prescribed=True)
        b=read_run(directory)
        for key in ('particles','curve_nodes','particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(a[key])==normalized(b[key]),'live shock AMR changed '+key)
        require(geometry_by_cycle(a)==geometry_by_cycle(b),'live shock AMR geometry differs')
        earliest={}
        for pattern in ('trace*.csv','birthtrace*.csv'):
            for path in (directory/'shock_events/curves').glob(pattern):
                generation=int(re.fullmatch(r'(?:birthtrace|trace)([0-9]+)\.rank[0-9]+\.csv',path.name)[1])
                for node,o in curve_shock_history.read_trace(path).items():
                    if node not in earliest or o['time']<earliest[node][1]['time']:
                        earliest[node]=(generation,o)
        births=[(g,o) for g,o in earliest.values() if o['birth']>0]
        for g,o in births:
            require(o['time']==o['birth'] and o['sequence']==0, 'newborn history not initialized at actual birth')
            require((o['flags']&1) if o['inside'] else o['phase']==0,
                    'newborn inherited encounter or lacks incomplete-entry flag')
        inside=[(g,o) for g,o in births if o['inside']]
        outside=[(g,o) for g,o in births if not o['inside']]
        require(inside and outside,'AMR lacks inside/outside positive-time birth observations')
        generation=min(g for g,o in inside)
        shared=work/'amr-shock-restart-journal'
        seed=run('amr_shock_seed',changes=dict(options,**{'time/nlim':generation,
            'output1/dcycle':-1,'shock_tracking/directory':str(shared)}),prescribed=True)
        cp=next((seed/'rst').glob('*.rst'))
        resumed=run('amr_shock_resume',restart=cp,overrides=('time/nlim=-1',),prescribed=True)
        prefix=read_run(seed,shock_directory=str(shared));suffix=read_run(resumed,shock_directory=str(shared))
        for key in ('particles','curve_nodes'):
            require(normalized(prefix[key]|suffix[key])==normalized(b[key]),'AMR newborn restart changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(suffix[key])==normalized(b[key]),'AMR newborn restart changed '+key)
        require(geometry_by_cycle(suffix)[max(suffix['curve_times'])]==geometry_by_cycle(b)[max(b['curve_times'])],
                'AMR newborn restart changed final geometry')
        def birth_rows(folder):
            return normalized({(int(re.fullmatch(r'birthtrace([0-9]+)\.rank[0-9]+\.csv',p.name)[1]),node):row
                for p in folder.glob('birthtrace*.csv')
                for node,row in curve_shock_history.read_trace(p).items()})
        require(birth_rows(shared/'curves')==birth_rows(directory/'shock_events/curves'),
                'AMR newborn restart changed complete birth traces')
        return dict(inside_births=len(inside),outside_births=len(outside),checkpoint_cycle=generation,
                    separate_birthtrace_files=len(list((directory/'shock_events/curves').glob('birthtrace*.csv'))),
                    all_records_exact=True,newborn_restart_exact=True)
    check('TO21_AMR04_RST02_live_shocks_insertion_restart',live_shock_birth_restart)
    save()
    print(json.dumps(report['checks'],indent=2))
    return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':
    raise SystemExit(main())
