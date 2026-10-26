"""Sample-at-output versus opt-in binary encounters: native local qualification.

Accepted planar, live-AMR and forced-box inputs are reused unchanged in physics.
All commands, inputs, binary fingerprints, failed cases and pending packets remain.
"""
import argparse
from collections import Counter
import copy
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import edit_input,provenance,sha
from tracer_history import read_run
import curve_shock_history
import shock_history
import shock_binary
import material_curves
from curve_binary_check import normalized,core

BASE='rho,p,t,s,v1,v2,v3'
SHOCK=['shock_valid','shock_region','shock_mach_pressure','shock_mach_temperature',
       'shock_mach_pressure_valid','shock_mach_temperature_valid']


def require(condition,message):
    if not condition:raise ValueError(message)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--prescribed-executable',type=Path,required=True)
    p.add_argument('--reference-runs',type=Path,required=True,help='previous CRV16 accepted runs')
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args();work=a.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    reference=a.reference_runs.resolve()
    report=dict(checks=[],runs=[],runner_sha256=sha(__file__),scope='local Serial/double',
                scientific_approval='pending',limits=['MPI runtime/GPU/production budgets unqualified'])
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save();print(row,flush=True)
    def run(name,fixture='planar',full=(False,False),changes=None,restart=None,
            overrides=(),reject=None,source=None):
        d=work/name;d.mkdir()
        input_source=Path(source) if source else reference/fixture/'input.athinput'
        binary=(a.prescribed_executable if fixture=='amr_shock' else a.executable).resolve()
        options={'output1/dcycle':0,'output2/dcycle':5,'output2/variables':BASE+',shock',
          'output3/dcycle':0,'output4/dcycle':0,'output5/dcycle':0,'output6/dcycle':0,
          'shock_tracking/trace':'false','output7/file_type':'prtcl_thermo_history',
          'output7/mode':'encounters','output7/source':'particles','output7/id':'particle_events',
          'output7/dcycle':11 if full[0] else 0,'output8/file_type':'crv',
          'output8/mode':'encounters','output8/id':'curve_events','output8/dcycle':17 if full[1] else 0}
        options.update(changes or {})
        if options.get('output2/dcycle') is None:
            text=input_source.read_text();block=None;lines=[]
            for line in text.splitlines():
                m=re.fullmatch(r'\s*<([^>]+)>\s*',line)
                if m:block=m[1]
                if block=='output2' and re.match(r'\s*dcycle\s*=',line):continue
                lines.append(line)
            text='\n'.join(lines)+'\n'
        else:text=input_source.read_text()
        path=d/'input.athinput';path.write_text(edit_input(text,{k:v for k,v in options.items() if v is not None}))
        command=[str(binary),'-r' if restart else '-i',str(restart or path),*overrides]
        with (d/'run.log').open('w') as f:
            code=subprocess.run(command,cwd=d,stdout=f,stderr=subprocess.STDOUT,timeout=180).returncode
        report['runs'].append(dict(name=name,command=command,exit=code,input_used=not bool(restart),
             input_sha256=sha(path),restart_sha256=sha(restart) if restart else None,**provenance(binary)));save()
        require((code!=0 and reject in (d/'run.log').read_text()) if reject else code==0,'native result '+name)
        return d
    def calls(directory):
        log=(directory/'run.log').read_text()
        def number(pattern):
            m=re.search(pattern+r'=(\d+)',log);return int(m[1]) if m else 0
        return number('particle detector cache calls this process'),number('curve stencil batches this process')
    def exact_samples(left,right):
        for key in ('particles','curve_nodes'):
            require(all(normalized(o)==normalized(right[key].get(k)) for k,o in left[key].items()),
                    'sample properties/time/ID changed: '+key)
    def histories(current,old):
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(current[key])==normalized(old[key]),'complete history changed: '+key)
    def checkpoints(directory):
        return {material_curves.read(str(p)+'.curves')['cycle']:p for p in sorted((directory/'rst').glob('*.rst'))}
    def copy_prefix(source,target,through):
        shutil.copytree(source/'shock_events',target)
        for path in list(target.rglob('*')):
            if path.suffix in ('.thp','.crv') and shock_binary.read(path)['generation']>through:path.unlink()
    def pending_packet(checkpoint,curve,folder):
        path=Path(str(checkpoint)+('.curves.shock.rank0' if curve else '.shock.rank0'))
        data=path.read_bytes();offset=data.find(shock_binary.MAGIC)
        require(offset>=0,'checkpoint lacks pending packet')
        result=folder/('curve.pending.crv' if curve else 'particle.pending.thp')
        result.write_bytes(data[offset:-8]);packet=shock_binary.read(result)
        require(packet['kind']==2,'checkpoint packet kind mismatch')
        return packet
    cases={}
    def sampled():
        d=run('sampled');cases['sampled']=d;x=read_run(d);clocks=len(x['curve_times'])
        require(calls(d)==(clocks,clocks),'snapshot mode evaluated between saved clocks')
        require(not (d/'shock_events').exists(),'snapshot mode created encounter journals')
        exact_samples(x,read_run(reference/'planar'))
        return dict(completed_steps=max(x['curve_times']),saved_clocks=clocks,
                    particle_cache_calls=calls(d)[0],curve_stencil_batches=calls(d)[1])
    check('SM01_TRC_sampling_only_at_snapshot_clocks',sampled)
    def no_fields():
        d=run('no_shock_fields',changes={'output2/variables':BASE});x=read_run(d)
        require(calls(d)==(0,0) and not (d/'shock_events').exists(),'unused detector still evaluated')
        for key in ('particles','curve_nodes'):
            baseline=read_run(cases['sampled'])[key]
            for k,o in x[key].items():
                old=copy.deepcopy(baseline[k]);old['fields']={f:v for f,v in old['fields'].items() if f not in SHOCK}
                require(normalized(o)==normalized(old),'turning off unused detector changed fluid samples')
        return dict(particle_calls=0,curve_calls=0)
    check('SM02_CFG_unused_detector_zero_work',no_fields)
    def complete():
        d=run('complete',full=(True,True));cases['complete']=d;x=read_run(d)
        old=read_run(reference/'planar');histories(x,old)
        require(calls(d)[0]==108 and calls(d)[1]>=108,'complete mode skipped timestep observations')
        require(not list((d/'shock_events').rglob('*.csv')),'completed CSV journal remains')
        exact_samples(x,old);require(len(x['curve_times'])==23,'event blocks duplicated snapshots')
        counts={}
        for carrier,folder in [(0,d/'shock_events'),(1,d/'shock_events/curves')]:
            packets=[shock_binary.read(p) for p in shock_binary.binary_files(folder,carrier)]
            require(all('fields' not in p for p in packets),'event output duplicated thermo snapshots')
            counts[str(carrier)]=sum(p['kind']==0 for p in packets)
        return dict(particle_events=len(x['particle_events']),curve_events=len(x['curve_events']),
                    particle_cache_calls=calls(d)[0],curve_stencil_batches=calls(d)[1],batch_counts=counts)
    check('SM03_EVT_binary_complete_history_exact_CSV_parity',complete)
    def carrier_scope():
        for name,mode,expected in [('curve_only_history',(False,True),(23,131)),
                                    ('particle_only_history',(True,False),(108,23))]:
            d=run(name,full=mode);x=read_run(d)
            require(calls(d)==expected,'history activation leaked to other carrier')
            require(bool(x['particle_events'])==mode[0] and bool(x['curve_events'])==mode[1],
                    'wrong history carrier enabled')
        return dict(independent_carriers=True)
    check('SM04_CFG_independent_population_activation',carrier_scope)
    def event_only():
        d=run('events_only',full=(True,True),changes={'output2/dcycle':0});x=read_run(d)
        require(not x['particles'] and not x['curve_nodes'] and not x['geometry'],'event-only outputs wrote snapshots')
        histories(x,read_run(reference/'planar'))
        require(calls(d)==(108,108),'event-only mode evaluated snapshots')
        return dict(particle_events=len(x['particle_events']),curve_events=len(x['curve_events']),snapshots=0)
    check('SM05_IO_event_only_outputs_without_snapshots',event_only)
    def bytes_and_mesh():
        d=run('every_step_sampled',changes={'output2/dcycle':1,'output5/dcycle':1});x=read_run(d)
        old=read_run(reference/'planar');exact_samples(x,old)
        require(next((d/'prtcl_thermo_history').glob('*.thp')).read_bytes()==
                next((reference/'planar/prtcl_thermo_history').glob('*.thp')).read_bytes(),
                'ordinary particle THP bytes changed')
        from shock_output_check import mesh_frames,mesh_sample
        frames=mesh_frames(d,'shock_mesh')
        for key in ('particles','curve_nodes'):
            for (g,n),o in x[key].items():
                for field in SHOCK:
                    a,b=o['fields'][field],mesh_sample(frames[g],o,field)
                    require(a==b or (a!=a and b!=b),'sample/cache/native mesh mismatch')
        return dict(particle_bytes_exact=True,mesh_samples_exact=True)
    check('SM06_IO_existing_THP_bytes_and_mesh_sampling_parity',bytes_and_mesh)
    def active_restart():
        old=read_run(reference/'planar');stop=min(e['first']['generation'] for e in old['curve_events'])
        seed=run('active_seed',full=(True,True),changes={'time/nlim':stop,'output1/dcycle':stop,
                 'shock_tracking/trace':'true'})
        checkpoint=checkpoints(seed)[stop];x=read_run(seed)
        trace=curve_shock_history.read_trace(seed/'shock_events/curves'/f'trace{stop}.rank0.csv')
        require(any(o['phase']==1 for o in trace.values()),'checkpoint not active')
        shared=work/'active-journal';copy_prefix(seed,shared,stop)
        resume=run('active_resume',full=(True,True),restart=checkpoint,
                   overrides=['time/nlim=-1','shock_tracking/directory='+str(shared)])
        continued=read_run(resume,shock_directory=str(shared));histories(continued,old)
        merged=copy.deepcopy(continued)
        for key in ('particles','curve_nodes'):merged[key]=x[key]|continued[key]
        exact_samples(read_run(cases['sampled']),merged)
        # The terminal prefix must replay exactly from a zero-step restart.
        zero_dir=work/'active-zero-journal';copy_prefix(seed,zero_dir,stop)
        before={str(p.relative_to(zero_dir)):sha(p) for p in zero_dir.rglob('*') if p.is_file()}
        zero=run('active_zero',full=(True,True),restart=checkpoint,
                 overrides=['time/nlim='+str(stop),'shock_tracking/directory='+str(zero_dir)])
        after={str(p.relative_to(zero_dir)):sha(p) for p in zero_dir.rglob('*') if p.is_file()}
        require(before==after,'zero-step changed/duplicated event files')
        z=read_run(zero,shock_directory=str(zero_dir))
        for key in ('particles','curve_nodes'):
            require(normalized(z[key])==normalized({k:v for k,v in x[key].items() if k[0]==stop}),
                    'zero-step snapshot changed')
        expect=next(s for s in x['geometry'] if s['cycle']==stop)
        require(normalized(z['geometry'][0]['curves'])==normalized(expect['curves']),'zero-step topology changed')
        cases['active']=seed;cases['resume']=resume
        return dict(active_cycle=stop,zero_step_file_replay_exact=True,continuation_exact=True)
    check('SM07_RST_active_and_terminal_zero_step_replay',active_restart)
    def pending_restart():
        d=run('pending_seed',full=(True,True),changes={'output7/dcycle':100,'output8/dcycle':100,
                                                    'output1/dcycle':80})
        checkpoint=checkpoints(d)[80];folder=work/'pending-packets';folder.mkdir()
        packets=[pending_packet(checkpoint,c,folder) for c in (False,True)]
        require(packets[1]['count']>0,'checkpoint did not retain unflushed completed events')
        require(all(p['previous_flush_generation']==0 for p in packets),'pending checkpoint unexpectedly flushed')
        shared=work/'pending-journal';copy_prefix(d,shared,80)
        resume=run('pending_resume',restart=checkpoint,
                 overrides=['shock_tracking/directory='+str(shared)])
        x=read_run(resume,shock_directory=str(shared));histories(x,read_run(reference/'planar'))
        cases['pending']=d;cases['pending_checkpoint']=checkpoint
        return dict(checkpoint_cycle=80,pending_counts=[p['count'] for p in packets],
                    original_emission_clocks_preserved=True)
    check('SM08_RST_unflushed_completed_records_survive_checkpoint',pending_restart)
    def amr():
        sampled=run('amr_sampled','amr_shock');complete=run('amr_complete','amr_shock',full=(True,True))
        a0,b0=read_run(sampled),read_run(complete);old=read_run(reference/'amr_shock')
        exact_samples(a0,old);exact_samples(b0,old);histories(b0,old)
        require(all(n<0 for g,n in a0['curve_nodes']),'negative IDs not exercised')
        require(any(o['birth']>0 for o in a0['curve_nodes'].values()),'live births missing')
        require(calls(sampled)==(len(a0['curve_times']),len(a0['curve_times'])),'AMR sample mode tracked births between outputs')
        cases['amr']=complete
        return dict(sample_clocks=len(a0['curve_times']),negative_IDs=True,fresh_births=True,
                    completed_events_exact=True)
    check('SM09_AMR_births_migration_and_negative_ID_parity',amr)
    def forced():
        d=run('forced','forced',full=(True,True),changes={'output2/dcycle':10})
        x=read_run(d);old=read_run(reference/'forced');exact_samples(x,old);histories(x,old)
        cases['forced']=d
        return dict(curve_samples=len(x['curve_nodes']),particle_events=len(x['particle_events']),
                    curve_events=len(x['curve_events']),short_accepted_fixture=True)
    check('SM10_forced_box_engineering_coexistence',forced)
    def guards():
        run('disabled_history_blocks',changes={'output7/dcycle':0,'output8/dcycle':0})
        require(not (work/'disabled_history_blocks/shock_events').exists(),'inactive history blocks enabled tracking')
        run('reject_mc_complete',full=(True,False),changes={'particles/particle_type':'lagrangian_mc',
            'particles/pusher':'lagrangian_mc'},reject='Shock tracking requires')
        run('reject_duplicate',full=(True,False),changes={'output9/file_type':'prtcl_thermo_history',
            'output9/mode':'encounters','output9/dcycle':5},reject='only one active encounter')
        run('reject_disabled_detector',full=(True,False),changes={'shock_tracking/enabled':'false'},reject='enabled detector')
        run('reject_wrong_mode',changes={'output2/mode':'complete_typo'},reject='unknown output mode')
        checkpoint=checkpoints(cases['active'])[37]
        run('reject_disabled_restart',restart=checkpoint,overrides=['output7/dcycle=0','output8/dcycle=0'],
            reject='explicit encounter output block')
        text=(reference/'planar/input.athinput').read_text()
        text=re.sub(r'(?m)^lineage\s*=.*\n','',text)
        no_lineage=work/'without_lineage.athinput';no_lineage.write_text(text)
        d=run('sampled_without_lineage',source=no_lineage)
        require(not (d/'shock_events').exists(),'sampled mode fabricated a history lineage')
        run('reject_curve_history_without_lineage',full=(False,True),source=no_lineage,reject='lineage')
        return dict(disabled_blocks_no_work=True,MC_history_rejected=True,duplicate_rejected=True,
                    disabled_detector_rejected=True,removed_history_restart_rejected=True)
    check('SM11_CFG_opt_in_and_restart_guards',guards)
    def mc_sample():
        d=run('mc_sample',changes={'particles/particle_type':'lagrangian_mc',
             'particles/pusher':'lagrangian_mc','output5/dcycle':5})
        x=read_run(d);require(calls(d)==(len(x['particle_times']),len(x['curve_times'])),'MC sample mode tracked every step')
        require(not x['particle_events'] and not x['curve_events'],'MC sample mode created histories')
        from shock_output_check import mesh_frames,mesh_sample
        frames=mesh_frames(d,'shock_mesh')
        for (g,n),o in x['particles'].items():
            for field in SHOCK:
                a,b=o['fields'][field],mesh_sample(frames[g],o,field)
                require(a==b or a!=a and b!=b,'MC point-sampling differs from native mesh')
        return dict(samples=len(x['particles']),point_samples_only=True,no_continuous_history_claim=True)
    check('SM12_MC_sampled_properties_without_encounter_claim',mc_sample)
    def overflow():
        d=run('overflow',full=(True,False),changes={'output7/dcycle':100,'output1/dcycle':106,
             'shock_tracking/event_capacity':11},reject='unflushed encounter buffer overflow')
        checkpoint=checkpoints(d)[106];folder=work/'overflow-packet';folder.mkdir()
        packet=pending_packet(checkpoint,False,folder);require(packet['count']==4,'overflow checkpoint lost pending events')
        shared=work/'overflow-journal';copy_prefix(d,shared,106)
        resumed=run('overflow_resume',restart=checkpoint,
            overrides=['shock_tracking/directory='+str(shared),'shock_tracking/event_capacity=4096'])
        x=read_run(resumed,shock_directory=str(shared));old=read_run(reference/'planar')
        require(normalized(x['particle_events'])==normalized(old['particle_events']),'capacity recovery lost/duplicated events')
        return dict(bounded_overflow_fails=True,pending_before_failure=4,recovered_events=len(x['particle_events']))
    check('SM13_IO_bounded_queue_and_capacity_recovery',overflow)
    def coarse():
        d=run('coarse',full=(True,True),changes={'output2/dcycle':107});x=read_run(d)
        require(set(x['curve_times'])=={0,107},'coarse sample clocks differ')
        require(all(o['fields']['shock_region']==0 for o in x['curve_nodes'].values()),'coarse fixture should miss curve passages')
        require(len(x['curve_events'])==16,'complete history missed between-snapshot passages')
        cases['coarse']=d
        return dict(snapshot_clocks=2,sampled_curve_inside_flags=0,complete_curve_passages=16)
    check('SM14_EVT_complete_passages_between_saved_snapshots',coarse)
    def bits_zero():
        d=run('bits_zero',full=(True,True),changes={'material_curves/curve_id_bits':0,'curve1/id':1})
        x=read_run(d);require(x['curve_id_layout']==(32,0,True) and len(x['curve_events'])==16,'zero-bit history layout lost')
        return dict(full_32bit_namespace=True,events=16)
    check('SM15_IO_full_width_curve_ID_layout',bits_zero)
    def pending_corrupt():
        checkpoint=cases['pending_checkpoint'];original=Path(str(checkpoint)+'.curves.shock.rank0').read_bytes()
        pos=original.find(shock_binary.MAGIC);require(pos>=0,'missing pending packet')
        payload=bytearray(original[pos:-8]);offset=96
        for _ in range(3):
            count=int.from_bytes(payload[offset:offset+4],'little');offset+=4+count
        payload[offset:offset+8]=(0).to_bytes(8,'little')  # predates first valid inside evidence
        payload[-8:]=shock_binary.digest(payload[:-8]).to_bytes(8,'little')
        changed=bytearray(original);changed[pos:-8]=payload
        changed[-8:]=shock_binary.digest(changed[:-8]).to_bytes(8,'little')
        folder=work/'corrupt-checkpoint';folder.mkdir()
        for source in checkpoint.parent.glob(checkpoint.name+'*'):shutil.copy2(source,folder/source.name)
        target=folder/checkpoint.name;Path(str(target)+'.curves.shock.rank0').write_bytes(changed)
        shared=work/'corrupt-journal';copy_prefix(cases['pending'],shared,80)
        run('reject_pending_semantics',restart=target,overrides=['shock_tracking/directory='+str(shared)],
            reject='buffered observation clock')
        return dict(repaired_checksums_do_not_hide_bad_pending_evidence=True)
    check('SM16_RST_semantic_validation_of_pending_payload',pending_corrupt)
    def oracle_dt():
        options={'output2/file_type':'crv','output2/id':'sampled_curve','output2/dcycle':None,
            'output2/dt':.05,'output2/variables':'default,vmag,cs,mach,shock',
            'problem/flow':'acceleration','problem/a1':.3,'problem/a2':.2}
        d=run('oracle_dt','amr_shock',changes=options,source=ROOT/'inputs/tests/shocks-curves/smooth.athinput')
        x=read_run(d);times=list(x['curve_times'].values())
        require(len(times)==9 and times[0]==0 and times[-1]==.4,'native dt snapshot cadence changed')
        require(calls(d)==(0,9),'dt snapshots detected between saved clocks')
        maximum=0
        import math
        for o in x['curve_nodes'].values():
            t=o['time'];v=(.2+.3*t,-.1+.2*t,0);pressure=(4-max(map(abs,v)))**2/(5/3)
            expected=dict(density=1,pressure=pressure,temperature=pressure,entropy=math.log(pressure),
                internal_energy=pressure/(2/3),v1=v[0],v2=v[1],v3=0,vmag=math.hypot(*v),
                sound_speed=math.sqrt(5*pressure/3),mach=math.hypot(*v)/math.sqrt(5*pressure/3))
            maximum=max(maximum,max(abs(o['fields'][key]-value) for key,value in expected.items()))
        require(maximum<=2e-12,'independent CT thermo oracle failed')
        cases['oracle']=d
        return dict(saved_times=times,max_thermo_error=maximum,unchanged_CT_budget=2e-12)
    check('SM17_NUM_independent_accelerating_thermo_and_dt_oracle',oracle_dt)
    def final_only():
        d=run('final_only',full=(True,True),changes={'output7/dcycle':-1,'output8/dcycle':-1})
        x=read_run(d);histories(x,read_run(reference/'planar'))
        require(len(list((d/'shock_events').glob('*.cycle*.thp')))==2 and
                len(list((d/'shock_events/curves').glob('*.cycle*.crv')))==2,
                'final-only cadence wrote intermediate event batches')
        return dict(complete_history=True,completed_batch_clocks=[0,107])
    check('SM18_IO_native_initial_final_only_flush_cadence',final_only)
    def tail_guard():
        seed=cases['active'];same=[p for p in sorted((seed/'rst').glob('*.rst'))
                                  if material_curves.read(str(p)+'.curves')['cycle']==37]
        require(len(same)>=2,'fixture lacks scheduled and terminal same-cycle checkpoints')
        ordinary=same[0];dirty=work/'tail-journal';copy_prefix(seed,dirty,37)
        before={str(p.relative_to(dirty)):sha(p) for p in dirty.rglob('*') if p.is_file()}
        run('reject_nonprefix_tail',restart=ordinary,overrides=['time/nlim=-1','shock_tracking/directory='+str(dirty)],
            reject='non-prefix committed tail')
        require(before=={str(p.relative_to(dirty)):sha(p) for p in dirty.rglob('*') if p.is_file()},
                'rejected restart mutated existing histories')
        packets=work/'ordinary-packets';packets.mkdir()
        limits={int(c):pending_packet(ordinary,c,packets)['previous_flush_generation'] for c in (False,True)}
        clean=work/'ordinary-journal';copy_prefix(seed,clean,37)
        for path in list(clean.rglob('*')):
            if path.suffix in ('.thp','.crv'):
                packet=shock_binary.read(path)
                if packet['generation']>limits[packet['carrier']]:path.unlink()
        resumed=run('ordinary_resume',restart=ordinary,overrides=['time/nlim=-1','shock_tracking/directory='+str(clean)])
        histories(read_run(resumed,shock_directory=str(clean)),read_run(reference/'planar'))
        return dict(orphan_tail_rejected_without_mutation=True,ordinary_checkpoint_clean_prefix_exact=True)
    check('SM19_RST_reject_forked_tail_and_resume_clean_ordinary_prefix',tail_guard)
    def store_unit():
        binary=a.executable.resolve().parents[1]/'shock_history_test';results=[]
        for mode,error in [('',None),('overflow','unflushed encounter buffer overflow'),
                           ('duplicate','duplicate buffered event identity'),('corrupt','checksum/schema mismatch')]:
            d=work/('unit_'+(mode or 'pass'));d.mkdir()
            command=[str(binary),str(d)]+([mode] if mode else [])
            with (d/'run.log').open('w') as log:
                code=subprocess.run(command,cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=30).returncode
            report['runs'].append(dict(name='native_unit_'+mode,command=command,exit=code,
                executable=str(binary),executable_sha256=sha(binary)));save()
            require((code!=0 and error in (d/'run.log').read_text()) if error else code==0,'native store unit failed '+mode)
            results.append(mode or 'roundtrip')
        node=shock_binary.read(next((work/'unit_pass').glob('*.crv')))['events'][0]
        require(node['node_id']==-1 and node['sequence']==9007199254740993,'typed identity/sequence lost bits')
        return dict(native_unit_modes=results,destroyed_ID_pending_roundtrip=True,
                    native_live_outflow_qualification=False,negative_one_and_large_sequence_exact=True)
    check('SM20_IO_native_pending_store_and_destroyed_identity_protocol',store_unit)
    def shared_directories():
        for name,folder in [('shared_root','.'),('shared_thp_directory','prtcl_thermo_history')]:
            d=run(name,full=(True,True),changes={'shock_tracking/directory':folder})
            x=read_run(d,shock_directory=folder);old=read_run(reference/'planar')
            histories(x,old);exact_samples(x,old)
            require(len(x['curve_times'])==23,'mixed directory lost snapshots')
        return dict(shared_curve_directory=True,shared_particle_directory=True,
                    recognized_magic_dispatch=True,unknown_magic_still_strict=True)
    check('SM21_IO_custom_roots_with_shared_snapshot_event_extensions',shared_directories)


    save();return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
