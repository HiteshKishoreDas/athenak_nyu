"""Native configurable-ID layout checks with retained inputs, journals and restart fixtures.

Reuses accepted smooth and planar/AMR fixtures. No physical tolerance tuning.
"""
import argparse
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import edit_input,provenance,sha
from tracer_history import read_run
from plot_curves import select_curve
import material_curves
from tracer_thermo_history import read_curves
from check import normalized,geometry_by_cycle
from id_check import seal


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--prescribed-executable',type=Path,required=True)
    parser.add_argument('--reference-runs',type=Path,required=True,help='default-eight-bit check.py output')
    parser.add_argument('--workdir',type=Path,required=True)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    reference=args.reference_runs.resolve()
    report=dict(checks=[],runs=[],runner_sha256=sha(__file__),scope='local Serial/double',scientific_approval='pending')
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def require(test,message):
        if not test:raise ValueError(message)
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save()
    def run(name,bits=0,curve=1,fixture='smooth',changes=None,restart=None,overrides=(),reject=None):
        directory=work/name;directory.mkdir()
        binary=(args.executable if fixture=='planar' else args.prescribed_executable).resolve()
        source=(reference/fixture/'input.athinput' if fixture in ('planar','amr_shock') else
                ROOT/'inputs/tests/shocks-curves/smooth.athinput')
        options={'material_curves/curve_id_bits':bits,'curve1/id':curve}
        if fixture=='smooth':
            options.update({'time/nlim':4,'shock_tracking/enabled':'true','curve1/count':2,
                'material_curves/adaptive':'false','output3/file_type':'prtcl_thermo_history',
                'output3/source':'curves','output3/variables':'rho,p','output3/dcycle':1,
                'output1/dt':-1,'output1/dcycle':-1})
        options.update(changes or {})
        input=directory/'input.athinput';input.write_text(edit_input(source.read_text(),options))
        if restart:
            archived=directory/'checkpoint_used';archived.mkdir()
            for path in restart.parent.glob(restart.name+'*'):shutil.copy2(path,archived/path.name)
        command=[str(binary),'-r' if restart else '-i',str(restart or input),*overrides]
        with (directory/'run.log').open('w') as stream:
            code=subprocess.run(command,cwd=directory,stdout=stream,stderr=subprocess.STDOUT,timeout=90).returncode
        report['runs'].append(dict(name=name,command=command,exit=code,input_used=not bool(restart),
            input_sha256=sha(input),restart_sha256=sha(restart) if restart else None,**provenance(binary)));save()
        require((code!=0 and reject in (directory/'run.log').read_text()) if reject else code==0,name+' unexpected native result')
        return directory
    def layouts():
        result={}
        for bits in (0,1,4,8,16,30):
            cid=(1<<bits)-1 if bits else 1
            directory=run('bits_'+str(bits),bits,cid)
            data=read_run(directory);nodes={p for g,p in data['curve_nodes']}
            require(data['curve_id_layout']==(32,bits,True),'missing or wrong declared layout')
            require(len(nodes)==2 and all(material_curves.node_curve_id(p,curve_id_bits=bits,signed=True)==cid for p in nodes),'wrong identity decoding')
            require(select_curve(data,node_id=min(nodes))['curve_id']==cid,'CLI selection inferred wrong curve')
            require(not data['curve_events'] and len(data['geometry'])>0,'empty-journal fixture missing')
            require('node ID bits: curve='+str(bits) in (directory/'run.log').read_text(),'startup split summary missing')
            result[str(bits)]=dict(curve_id=cid,node_ids=sorted(nodes))
        return result
    check('VB01_declared_layouts_0_1_4_8_16_30_empty_events',layouts)
    def guards():
        for label,value,reason in [('negative','-1','unsigned decimal integer'),('text','garbage','unsigned decimal integer'),
                                   ('fraction','1.5','unsigned decimal integer'),('overflow','4294967296','curve_id_bits must lie'),
                                   ('31',31,'curve_id_bits must lie'),('32',32,'curve_id_bits must lie')]:
            run('bad_bits_'+label,bits=value,reject=reason)
        run('zero_multiple',changes={'material_curves/ncurves':2},reject='ncurves exceeds curve_id_bits namespace')
        run('zero_wrong_curve',curve=2,reject='curve IDs must be unique within curve_id_bits namespace')
        return dict(rejected_inputs=8)
    check('VB02_checked_input_and_single_curve_guards',guards)
    zero={}
    def remapped(data,bits=0,cid=1):
        data=copy.deepcopy(data)
        def node(old):
            if not old:return 0
            value=((cid<<(32-bits)) if bits else 0)|(old&((1<<24)-1))
            return value if value<1<<31 else value-(1<<32)
        data['curve_nodes']={(g,node(p)):dict(o,node_id=node(p),curve_id=cid) for (g,p),o in data['curve_nodes'].items()}
        for key in ('curve_events','curve_suspended'):
            for e in data[key]:
                e['particle_id']=e['node_id']=node(e['node_id']);e['curve_id']=cid
                e['identity']=(*e['identity'][:2],e['node_id'],e['identity'][3])
        for state in data['geometry']:
            for curve in state['curves']:
                curve['id']=cid
                for n in curve['nodes']:
                    for key in ('id','prev','next','parent_left','parent_right'):n[key]=node(n[key])
        return data
    def physical(fixture):
        directory=run('zero_'+fixture,fixture=fixture)
        data=read_run(directory);baseline=read_run(reference/fixture)
        require(baseline['curve_id_layout'] in ((32,8),(32,8,True)),'reference must use default eight curve bits')
        expected=remapped(baseline)
        for key in ('particles','curve_nodes','particle_times','curve_times','particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(data[key])==normalized(expected[key]),fixture+' changed '+key)
        require(geometry_by_cycle(data)==geometry_by_cycle(expected),fixture+' geometry changed')
        require(data['curve_id_layout']==(32,0,True),'zero layout not preserved')
        zero[fixture]=directory
        return dict(curve_samples=len(data['curve_nodes']),curve_events=len(data['curve_events']),all_records_exact_after_ID_remap=True)
    check('VB03_zero_bit_planar_particles_and_real_encounters',lambda:physical('planar'))
    check('VB04_zero_bit_live_AMR_births_and_histories',lambda:physical('amr_shock'))
    def restart():
        shared=work/'restart-journal'
        seed=run('restart_seed',fixture='amr_shock',changes={'time/nlim':13,'output1/dcycle':-1,'shock_tracking/directory':str(shared)})
        checkpoint=next((seed/'rst').glob('*.rst'))
        resumed=run('restart_resume',fixture='amr_shock',restart=checkpoint,overrides=('time/nlim=-1',))
        a=read_run(seed,shock_directory=str(shared));b=read_run(resumed,shock_directory=str(shared))
        expected=read_run(zero['amr_shock'])
        for key in ('particles','curve_nodes','particle_times','curve_times'):
            require(normalized(a[key]|b[key])==normalized(expected[key]),'zero-bit restart changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(b[key])==normalized(expected[key]),'zero-bit restart changed '+key)
        run('restart_changed_bits',fixture='amr_shock',restart=checkpoint,
            overrides=('material_curves/curve_id_bits=8','time/nlim=-1'),reject='restart curve_id_bits mismatch')
        return dict(checkpoint_cycle=13,exact_continuation=True,changed_split_rejected=True)
    check('VB05_zero_bit_active_AMR_restart_split_guard',restart)
    def stream_mismatch():
        target=work/'mixed_geometry';target.mkdir()
        src=zero['planar']
        path=next(src.glob('*.curves.*.txt'));shutil.copy2(path,target/path.name)
        shutil.copytree(reference/'planar/curve_thermo_history',target/'curve_thermo_history')
        try:read_run(target)
        except ValueError as exc:require('layout differs' in str(exc),'wrong geometry/thermo mismatch rejection')
        else:raise ValueError('mixed geometry/thermo layouts accepted')
        target=work/'mixed_empty_journal';target.mkdir()
        path=next((work/'bits_0').glob('*.curves.*.txt'));shutil.copy2(path,target/path.name)
        shutil.copytree(work/'bits_8/shock_events',target/'shock_events')
        try:read_run(target)
        except ValueError as exc:require('layout differs' in str(exc),'wrong empty-journal mismatch rejection')
        else:raise ValueError('empty journal layout not checked')
        folder=work/'mixed_thermo_stream';folder.mkdir()
        for directory,cycle in [('bits_0',0),('bits_8',1)]:
            source=next((work/directory/'curve_thermo_history').glob(f'*.cycle{cycle}.rank0.csv'))
            shutil.copy2(source,folder/source.name)
        try:read_curves(folder)
        except ValueError as exc:require('ID encoding changed' in str(exc),'wrong within-stream mismatch rejection')
        else:raise ValueError('changing thermo layout accepted')
        return dict(geometry_thermo_empty_journal_and_within_stream_mismatch_rejected=True)
    check('VB06_cross_stream_layout_consistency',stream_mismatch)
    def compatibility():
        plain=run('plain_v4',bits=8,curve=255,changes={'shock_tracking/enabled':'false'})
        plain_cp=next((plain/'rst').glob('*.rst'));plain_side=Path(str(plain_cp)+'.curves')
        shutil.copy2(plain_side,plain/'original-v4.curves')
        body=plain_side.read_bytes();lines=body[:body.rfind(b'commit ')].decode().splitlines()
        converted=[]
        for line in lines:
            f=line.split()
            if f[0]=='id_bits':continue
            if f[0]=='MATERIAL_CURVES':line='MATERIAL_CURVES 3 ASCII 64 U32 I64 I32'
            elif f[0]=='curve_columns':line=line.replace('id:i32','id:u32').replace('next_sequence:u64','next_sequence:u32')
            elif f[0]=='node_columns':
                for name in ('id','prev','next','parent_left','parent_right'):line=line.replace(name+':i32',name+':u32')
            elif f[0]=='node':
                for j in range(1,6):f[j]=str(int(f[j])&((1<<32)-1))
                line=' '.join(f)
            converted.append(line)
        lines=converted
        plain_side.write_bytes(seal('\n'.join(lines)+'\n'))
        run('plain_v3_resume',restart=plain_cp,overrides=('time/nlim=5',))
        return dict(default_v3_geometry_restart=True)
    check('VB07_fixed_8bit_geometry_compatibility',compatibility)
    def corrupt_layout_metadata():
        seed=work/'restart_seed';cp=next((seed/'rst').glob('*.rst'))
        artifacts=work/'invalid_layout_artifacts';shutil.copytree(seed/'rst',artifacts)
        checkpoint=artifacts/cp.name;side=Path(str(checkpoint)+'.curves')
        raw=side.read_bytes();body=raw[:raw.rfind(b'commit ')].decode()
        require('id_bits 0\n' in body,'metadata fixture ineffective')
        side.write_bytes(seal(body.replace('id_bits 0\n','id_bits 31\n',1)))
        run('bad_checkpoint_layout',fixture='amr_shock',restart=checkpoint,
            overrides=('time/nlim=-1',),reject='invalid curve_id_bits metadata')
        folder=work/'invalid_thermo_layout';folder.mkdir()
        source=next((work/'bits_0/curve_thermo_history').glob('*.csv'))
        lines=source.read_text().splitlines();head=lines[0].split();head[-1]='31'
        lines[0]=' '.join(head);(folder/source.name).write_text('\n'.join(lines)+'\n')
        try:read_curves(folder)
        except ValueError as exc:require('curve_id_bits metadata' in str(exc),'wrong thermo metadata rejection')
        else:raise ValueError('invalid thermo layout accepted')
        return dict(repaired_checksum_checkpoint_and_thermo_layout_rejected=True)
    check('VB08_corrupt_layout_metadata',corrupt_layout_metadata)
    def signed_shock_restart():
        directory=run('signed_amr_shock',bits=8,curve=255,fixture='amr_shock')
        data=read_run(directory);expected=remapped(read_run(reference/'amr_shock'),bits=8,cid=255)
        for key in ('particles','curve_nodes','particle_times','curve_times','particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(data[key])==normalized(expected[key]),'signed AMR changed '+key)
        require(geometry_by_cycle(data)==geometry_by_cycle(expected),'signed AMR geometry changed')
        require(all(p<0 for g,p in data['curve_nodes']) and data['curve_events'], 'negative-ID shock fixture ineffective')
        shared=work/'signed-restart-journal'
        seed=run('signed_restart_seed',bits=8,curve=255,fixture='amr_shock',changes={'time/nlim':13,'output1/dcycle':-1,'shock_tracking/directory':str(shared)})
        cp=next((seed/'rst').glob('*.rst'))
        resumed=run('signed_restart_resume',bits=8,curve=255,fixture='amr_shock',restart=cp,overrides=('time/nlim=-1',))
        a=read_run(seed,shock_directory=str(shared));b=read_run(resumed,shock_directory=str(shared))
        for key in ('particles','curve_nodes','particle_times','curve_times'):
            require(normalized(a[key]|b[key])==normalized(data[key]),'negative-ID restart changed '+key)
        for key in ('particle_events','curve_events','particle_suspended','curve_suspended'):
            require(normalized(b[key])==normalized(data[key]),'negative-ID restart changed '+key)
        require(select_curve(data,node_id=min(p for g,p in data['curve_nodes']))['curve_id']==255,'negative CLI decoding failed')
        return dict(curve_samples=len(data['curve_nodes']),curve_events=len(data['curve_events']),negative_IDs=True,all_records_exact_after_ID_remap=True,restart_exact=True)
    check('VB09_signed_negative_ID_live_AMR_events_and_restart',signed_shock_restart)
    save();print(json.dumps(report['checks'],indent=2));return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
