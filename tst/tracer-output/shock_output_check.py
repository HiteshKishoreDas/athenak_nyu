"""Native THP/BIN shock-output checks on accepted planar and live-AMR fixtures.

Run with the supplied vis/.venv Python. No new physical parameters or tolerances.
All native commands, inputs, output files and rejected fixtures are retained.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
import numpy as np
from bin_convert import read_binary
from common import edit_input,provenance,sha
from tracer_history import read_run
from tracer_thermo_history import read_particles,read_curves
from shock_history import read_trace as particle_trace
from curve_shock_history import read_trace as curve_trace
from check import normalized,geometry_by_cycle

FIELDS=['shock_valid','shock_region','shock_mach_pressure','shock_mach_temperature',
        'shock_mach_pressure_valid','shock_mach_temperature_valid']
BASE='rho,p,t,s,v1,v2,v3'


def require(condition,message):
    if not condition:raise ValueError(message)


def equal(a,b):
    return a==b or (isinstance(a,float) and isinstance(b,float) and math.isnan(a) and math.isnan(b))


def mesh_frames(directory,output_id):
    result={}
    for path in sorted((Path(directory)/'bin').glob('*.'+output_id+'.*.bin')):
        frame=read_binary(path);cycle=frame['cycle']
        if cycle in result:
            for name,values in frame['mb_data'].items():
                require(np.array_equal(values,result[cycle]['mb_data'][name],equal_nan=True),'same-cycle BIN differs')
        result[cycle]=frame
    require(result,'missing BIN stream '+output_id)
    return result


def mesh_sample(frame,observation,name):
    # Single-rank full-volume native BIN keeps the complete global-GID block order.
    # Bound checks independently verify that each stored sample belongs to that block.
    gid=observation['gid'];values=frame['mb_data'][name]
    require(0<=gid<len(values),'sample GID outside native BIN')
    bounds=frame['mb_geometry'][gid];indices=[]
    for axis,n in enumerate(values[gid].shape[::-1]):
        low,high=bounds[2*axis:2*axis+2]
        x=observation['x'][axis]
        domain_low,domain_high=frame[f'x{axis+1}min'],frame[f'x{axis+1}max']
        if n>1:x=domain_low+(x-domain_low)%(domain_high-domain_low)
        index=0 if n==1 else math.floor((x-low)/((high-low)/n))
        require(0<=index<n,'sample outside native leaf bounds')
        indices.append(index)
    i,j,k=indices
    return float(values[gid][k,j,i])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--prescribed-executable',type=Path,required=True)
    parser.add_argument('--reference-runs',type=Path,required=True,help='check.py accepted-input runs')
    parser.add_argument('--workdir',type=Path,required=True)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    reference=args.reference_runs.resolve()
    report=dict(checks=[],runs=[],runner_sha256=sha(__file__),scope='local Serial/double',scientific_approval='pending')
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save()
    def run(name,fixture='planar',fields=True,mesh=True,curve_only=False,changes=None,restart=None,overrides=(),reject=None):
        directory=work/name;directory.mkdir()
        binary=(args.executable if fixture in ('planar','forced') else args.prescribed_executable).resolve()
        source=(ROOT/'inputs/tests/shocks-curves/forced.athinput' if fixture=='forced' else reference/fixture/'input.athinput').read_text()
        if curve_only:
            source=re.sub(r'(?ms)^<(?:particles|tracer_seed[^>]*)>.*?(?=^<|\Z)','',source)
        options={'output2/variables':BASE+(',shock' if fields else ''),
            'output2/source':'curves' if curve_only else 'both','shock_tracking/trace':'true' if fields else 'false',
            'output5/file_type':'bin','output5/variable':'shock','output5/id':'shock_mesh',
            'output5/data_precision':'real','output5/dcycle':1 if mesh else 0,'output5/dt':-1,
            'output6/file_type':'bin','output6/variable':'hydro_w','output6/id':'fluid',
            'output6/data_precision':'real','output6/dcycle':1,'output6/dt':-1}
        if fixture=='forced':
            options.update({'output2/dcycle':10,'output5/dcycle':10 if mesh else 0,'output6/dcycle':10})
        if fixture=='amr_shock':options['curve1/id']=255
        options.update(changes or {})
        path=directory/'input.athinput';path.write_text(edit_input(source,options))
        if restart:
            archive=directory/'checkpoint_used';archive.mkdir()
            for p in restart.parent.glob(restart.name+'*'):shutil.copy2(p,archive/p.name)
        command=[str(binary),'-r' if restart else '-i',str(restart or path),*overrides]
        with (directory/'run.log').open('w') as log:
            code=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=120).returncode
        report['runs'].append(dict(name=name,command=command,exit=code,input_used=not bool(restart),input_sha256=sha(path),restart_sha256=sha(restart) if restart else None,**provenance(binary)));save()
        require((code!=0 and reject in (directory/'run.log').read_text()) if reject else code==0,'unexpected native outcome: '+name)
        return directory
    cases={}
    def parity(fixture):
        control=run(fixture+'_control',fixture,fields=False,mesh=False)
        directory=run(fixture,fixture);cases[fixture]=directory
        a,b=read_run(control),read_run(directory)
        for key in ('particles','curve_nodes'):
            actual=copy.deepcopy(b[key])
            for o in actual.values():o['fields']={k:v for k,v in o['fields'].items() if k not in FIELDS}
            require(normalized(actual)==normalized(a[key]),'shock sampling changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(a[key])==normalized(b[key]),'shock sampling changed '+key)
        require(geometry_by_cycle(a)==geometry_by_cycle(b),'shock output changed curve geometry')
        af,bf=mesh_frames(control,'fluid'),mesh_frames(directory,'fluid')
        require(af.keys()==bf.keys(),'fluid output clocks differ')
        for generation in af:
            for name in af[generation]['mb_data']:
                require(np.array_equal(af[generation]['mb_data'][name],bf[generation]['mb_data'][name]),'shock output changed fluid arrays')
        frames=mesh_frames(directory,'shock_mesh');counts={};missing=0
        for carrier,key,reader,folder in [('particle','particles',particle_trace,directory/'shock_events'),('curve','curve_nodes',curve_trace,directory/'shock_events/curves')]:
            traces={}
            for path in folder.glob('trace*.rank*.csv'):
                cycle=int(re.match(r'trace([0-9]+)',path.name)[1])
                if cycle in b['particle_times' if carrier=='particle' else 'curve_times']:
                    traces.update({(cycle,p):o for p,o in reader(path).items()})
            for identity,o in b[key].items():
                raw=traces[identity]
                expected=[int(raw['valid']),int(raw['inside']),raw['mach_pressure'],raw['mach_temperature'],int(raw['valid'] and math.isfinite(raw['mach_pressure'])),int(raw['valid'] and math.isfinite(raw['mach_temperature']))]
                for field,value in zip(FIELDS,expected):
                    require(equal(o['fields'][field],value),'THP differs from synchronized raw trace: '+field)
                    require(equal(mesh_sample(frames[o['cycle']],o,field),value),'BIN differs from THP/raw sample: '+field)
                missing+=int(math.isnan(expected[2]))+int(math.isnan(expected[3]))
            counts[carrier+'_samples']=len(b[key])
        require(missing>0,'fixture did not exercise unavailable Mach estimates')
        if fixture=='amr_shock':require(all(p<0 for _,p in b['curve_nodes']),'negative-ID AMR fixture ineffective')
        return dict(**counts,mesh_cycles=len(frames),unavailable_estimates=missing,all_fields_exact=True,fluid_and_histories_unchanged=True)
    check('SHOUT01_IO03_planar_THP_BIN_trace_exact',lambda:parity('planar'))
    check('SHOUT02_AMR04_IO03_negative_AMR_THP_BIN_trace_exact',lambda:parity('amr_shock'))
    def curve_only():
        directory=run('curve_only',curve_only=True)
        a=mesh_frames(cases['planar'],'shock_mesh');b=mesh_frames(directory,'shock_mesh')
        require(a.keys()==b.keys(),'curve-only mesh clocks differ')
        for generation in a:
            for field in FIELDS:require(np.array_equal(a[generation]['mb_data'][field],b[generation]['mb_data'][field],equal_nan=True),'curve-only detector differs from particle cache')
        data=read_run(directory);require(not data['particles'] and data['curve_nodes'],'curve-only source missing')
        return dict(mesh_fields_exact_without_particles=True,curve_samples=len(data['curve_nodes']))
    check('SHOUT03_curve_only_mesh_and_THP',curve_only)
    def restart():
        shared=work/'restart-journal'
        seed=run('restart_seed','amr_shock',changes={'time/nlim':13,'output1/dcycle':-1,'shock_tracking/directory':str(shared)})
        cp=next((seed/'rst').glob('*.rst'))
        resumed=run('restart_resume','amr_shock',restart=cp,overrides=('time/nlim=-1',))
        a=read_run(seed,shock_directory=str(shared));b=read_run(resumed,shock_directory=str(shared));expected=read_run(cases['amr_shock'])
        zero=run('restart_zero_step','amr_shock',restart=cp,overrides=('time/nlim=13',))
        z=read_run(zero,shock_directory=str(shared))
        for key in ('particles','curve_nodes'):
            before={identity:o for identity,o in a[key].items() if identity[0]==13}
            require(normalized(z[key])==normalized(before),'zero-step restart changed shock properties: '+key)
        before=mesh_frames(seed,'shock_mesh')[13];after=mesh_frames(zero,'shock_mesh')[13]
        for name in FIELDS:
            require(np.array_equal(before['mb_data'][name],after['mb_data'][name],equal_nan=True),
                    'zero-step restart changed mesh diagnostics: '+name)
        for key in ('particles','curve_nodes','particle_times','curve_times'):
            require(normalized(a[key]|b[key])==normalized(expected[key]),'restart changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(b[key])==normalized(expected[key]),'restart changed '+key)
        return dict(active_AMR_checkpoint_cycle=13,full_shock_properties_and_histories_exact=True,zero_step_properties_and_mesh_exact=True)
    check('SHOUT04_RST02_active_restart',restart)
    def guards():
        for label,change,reason in [
            ('disabled',{'shock_tracking/enabled':'false'},'shock fields require enabled'),
            ('ghosts',{'output5/ghost_zones':'true'},'active zones'),
            ('seed_target',{'tracer_seed_initial/target':'shock_region'},'not seed targets')]:
            run('guard_'+label,changes=dict(change,**{'time/nlim':0}),reject=reason)
        return dict(explicit_rejections=3)
    check('SHOUT05_CFG_guards',guards)
    def float32():
        directory=run('float32',changes={'output5/data_precision':'float32','time/nlim':4})
        actual=mesh_frames(directory,'shock_mesh');reference=mesh_frames(cases['planar'],'shock_mesh')
        for cycle,frame in actual.items():
            for name in FIELDS:
                require(np.asarray(frame['mb_data'][name]).dtype==np.float32,'native float32 output not exercised')
                expected=np.asarray(reference[cycle]['mb_data'][name],dtype=np.float32)
                require(np.array_equal(frame['mb_data'][name],expected,equal_nan=True),'float32 output differs from exact native cast')
        return dict(mesh_cycles=len(actual),native_float32_cast_exact=True)
    check('SHOUT06_IO04_native_float32_precision',float32)
    def empty():
        directory=run('empty',changes={'time/nlim':2,'material_curves/ncurves':0,
                                      'tracer_seed_initial/count_per_event':0})
        data=read_run(directory)
        require(not data['particles'] and not data['curve_nodes'],'empty population fixture ineffective')
        require(data['curve_id_layout']==(32,8,True) and data['curve_times']==data['particle_times'],
                'empty THP lost clocks or explicit layout')
        frames=mesh_frames(directory,'shock_mesh');expected=mesh_frames(cases['planar'],'shock_mesh')
        for cycle,frame in frames.items():
            for name in FIELDS:
                require(np.array_equal(frame['mb_data'][name],expected[cycle]['mb_data'][name],equal_nan=True),
                        'empty carrier populations changed detector mesh')
        return dict(empty_population_clocks=len(data['curve_times']),explicit_empty_layout=True,mesh_exact=True)
    check('SHOUT07_IO02_empty_populations',empty)
    def native_v2():
        directory=cases['amr_shock'];curves,_=read_curves(directory/'curve_thermo_history')
        total=0
        for path in (directory/'curve_thermo_history').glob('*.thp'):
            raw,_=read_particles(path)
            for key,o in raw.items():
                row=curves[key]
                require(o['tag']==row['node_id'] and o['seed_id']==row['curve_id'] and o['gid']==row['gid'],'THP integer slots lost IDs')
                require(o['fields']['birth']==row['birth'] and o['fields']['a']==row['a'],'THP birth/material label changed')
                total+=1
        return dict(existing_THP_v2_reader_records=total,negative_IDs_exact=True)
    check('THP01_IO04_existing_reader_and_exact_signed_slots',native_v2)
    check('SHOUT08_native_forced_3D_THP_BIN_trace_exact',lambda:parity('forced'))
    save();print(json.dumps(report['checks'],indent=2));return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
