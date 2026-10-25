"""Native compact-ID boundaries and legacy analysis checks; retain fixtures/logs.

Uses the accepted smooth-flow fixture. Counter/geometry edits are deliberately
synthetic persistence tests, not physical initial conditions or science results.
"""
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import sys
import traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import edit_input,provenance,sha
import material_curves
from material_curves import node_curve_id
from tracer_thermo_history import read_curves


def seal(body):
    checksum=material_curves.OFFSET
    data=body.encode('ascii')
    for byte in data:checksum=((checksum^byte)*material_curves.PRIME)&material_curves.MASK64
    return data+f'commit {checksum} {len(data)}\n'.encode()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True,help='material_curve_test build')
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--historical-run',type=Path,help='optional pre-compact-ID output')
    parser.add_argument('--curve-bits',type=int,choices=(0,8),default=8)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    binary=args.executable.resolve()
    limit=(1<<(32-args.curve_bits))-1
    curve_id=255 if args.curve_bits else 1
    report=dict(checks=[],runs=[],runner_sha256=sha(__file__),scope='local Serial/double',
                curve_id_bits=args.curve_bits,fixture='synthetic counter/link boundaries on accepted uniform flow',scientific_approval='pending')
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save()
    def require(test,message):
        if not test:raise ValueError(message)
    def run(name,changes=None,restart=None,overrides=(),reject=None):
        directory=work/name;directory.mkdir()
        options={'shock_tracking/enabled':'false','time/nlim':0,'curve1/id':curve_id,'material_curves/curve_id_bits':args.curve_bits,'curve1/count':2,
            'curve1/x1':-.4,'curve1/end1':.4,'curve1/x2':0,'curve1/end2':0,
            'material_curves/adaptive':'true','material_curves/max_spacing':.6,
            'material_curves/spacing_cells':1000000,'material_curves/max_angle':3,
            'material_curves/chord_error':1000,'material_curves/output_dcycle':1,
            'output1/dt':-1,'output1/dcycle':-1,'output3/file_type':'prtcl_thermo_history',
            'output3/source':'curves','output3/variables':'rho,p','output3/dcycle':1}
        options.update(changes or {})
        input=directory/'input.athinput';input.write_text(edit_input((ROOT/'inputs/tests/shocks-curves/smooth.athinput').read_text(),options))
        if restart:
            artifacts=directory/'checkpoint_used';artifacts.mkdir()
            for source in (restart,Path(str(restart)+'.curves')):
                shutil.copy2(source,artifacts/source.name)
        command=[str(binary),'-r' if restart else '-i',str(restart or input),*overrides]
        with (directory/'run.log').open('w') as stream:
            code=subprocess.run(command,cwd=directory,stdout=stream,stderr=subprocess.STDOUT,timeout=30).returncode
        report['runs'].append(dict(name=name,command=command,exit=code,input_used=not bool(restart),
            input_sha256=sha(input),restart_sha256=sha(restart) if restart else None,
            curve_sidecar_sha256=sha(str(restart)+'.curves') if restart else None,**provenance(binary)));save()
        require((code!=0 and reject in (directory/'run.log').read_text()) if reject else code==0,name+' unexpected exit; inspect run.log')
        return directory
    seed=run('seed')
    checkpoint=next((seed/'rst').glob('*.rst'))
    sidecar=Path(str(checkpoint)+'.curves')
    original=sidecar.read_bytes();state=material_curves.read(sidecar)
    require(state['version']==5 and len(state['curves'][0]['nodes'])==3,'seed not compact/adapted')

    def patch(counter,mode):
        rows=original[:original.rfind(b'commit ')].decode().splitlines()
        for index,line in enumerate(rows):
            f=line.split()
            if f[0]=='curve':f[9]=str(counter)
            elif f[0]=='node':
                if mode=='single' and f[8]=='0':f[10]='-0.3'
                elif mode=='batch':f[10]=str(2*float(f[10]))
            else:continue
            rows[index]=' '.join(f)
        sidecar.write_bytes(seal('\n'.join(rows)+'\n'))
        material_curves.read(sidecar,fluid_path=checkpoint)
    def last_id():
        patch(limit,'single')
        directory=run('last_id',restart=checkpoint,overrides=('time/nlim=1',))
        result=material_curves.read(directory/'material.curves.final.txt');curve=result['curves'][0]
        ids=[n['id'] for n in curve['nodes']]
        require(len(ids)==4 and -1 in ids,'last 32-bit bit pattern not allocated as signed -1')
        require(curve['next_sequence']==limit+1,'exhausted allocator sentinel missing')
        for n in curve['nodes']:
            require(node_curve_id(n['id'],curve_id_bits=args.curve_bits,signed=True)==curve_id,'wrong high byte')
        report['last_id_state']=str(directory/'material.curves.final.txt')
        return dict(node_ids=ids,last_id=-1,next_sequence=curve['next_sequence'])
    check('ID01_last_sequence_signed_minus1_and_links',last_id)
    def batch():
        patch(limit,'batch')
        directory=run('batch_limit',restart=checkpoint,overrides=('time/nlim=1',))
        curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
        require(len(curve['nodes'])==3 and curve['next_sequence']==limit and curve['reason']==2,
                'over-capacity split batch allocated IDs or missed quality flag')
        return dict(node_count=3,next_sequence=limit,reason=2)
    check('ID02_atomic_split_batch_exhaustion',batch)
    def exhausted():
        patch(limit+1,'batch')
        directory=run('exhausted',restart=checkpoint,overrides=('time/nlim=1',))
        curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
        require(len(curve['nodes'])==3 and curve['next_sequence']==limit+1 and curve['reason']==2,
                'exhausted counter wrapped/reused an identity')
        return dict(no_wrap_or_reuse=True)
    check('ID03_exhausted_allocator_no_wrap',exhausted)
    def inputs():
        for value in (0,256,4294967297):
            run('invalid_curve_'+str(value),changes={'curve1/id':value},reject='curve IDs must be unique within curve_id_bits namespace')
        run('too_many_curves',changes={'material_curves/ncurves':256},reject='ncurves exceeds curve_id_bits namespace')
        return dict(out_of_range_values=[0,256,4294967297],curve_count_limit=curve_id)
    check('ID04_config_reject_before_narrowing',inputs)
    def corrupt_links():
        rows=original[:original.rfind(b'commit ')].decode().splitlines()
        for index,line in enumerate(rows):
            f=line.split()
            if f[0]=='node':f[2]=str(1<<32);rows[index]=' '.join(f);break
        sidecar.write_bytes(seal('\n'.join(rows)+'\n'))
        run('wide_link',restart=checkpoint,overrides=('time/nlim=1',),reject='invalid exact-width integer')
        return dict(UINT32_overflow_rejected=True)
    check('ID05_persistent_link_overflow_rejected',corrupt_links)
    def legacy():
        sidecar.write_bytes(original.replace(b'MATERIAL_CURVES 5 ASCII',b'MATERIAL_CURVES 2 ASCII',1))
        run('old_checkpoint',restart=checkpoint,overrides=('time/nlim=1',),reject='restart requires v3/v4/v5 compact IDs')
        if args.historical_run:
            geometry=next(args.historical_run.glob('*.curves.*.txt'))
            old=material_curves.read(geometry)
            rows,_=read_curves(args.historical_run/'curve_thermo_history')
            require(old['version']==2 and rows and all(p>1<<32 for g,p in rows),'historical readers lost 64-bit IDs')
        return dict(old_restart_rejected=True,historical_analysis=bool(args.historical_run))
    check('ID06_legacy_analysis_and_explicit_restart_rejection',legacy)
    if args.curve_bits==0:
        def sign_crossing():
            patch(1<<31,'single')
            directory=run('sign_crossing',restart=checkpoint,overrides=('time/nlim=1',))
            curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
            ids=[n['id'] for n in curve['nodes']]
            require(len(ids)==4 and -(1<<31) in ids and curve['next_sequence']==(1<<31)+1,
                    'sign crossing lost allocation, ancestry or links')
            report['sign_crossing_state']=str(directory/'material.curves.final.txt')
            return dict(node_ids=ids,next_sequence=curve['next_sequence'])
        check('ID07_INT32_MIN_insertion_with_positive_parents',sign_crossing)
        def legacy_v4_full_sequence():
            checkpoint=next((work/'last_id/rst').glob('*.rst'))
            path=Path(str(checkpoint)+'.curves');original=path.read_bytes()
            rows=original[:original.rfind(b'commit ')].decode().splitlines()
            for i,line in enumerate(rows):
                f=line.split()
                if f[0]=='MATERIAL_CURVES':rows[i]='MATERIAL_CURVES 4 ASCII 64 U32 I64 I32'
                elif f[0]=='curve_columns':rows[i]=line.replace('id:i32','id:u32')
                elif f[0]=='node_columns':
                    for name in ('id','prev','next','parent_left','parent_right'):line=line.replace(name+':i32',name+':u32')
                    rows[i]=line
                elif f[0]=='node':
                    for j in range(1,6):f[j]=str(int(f[j])&((1<<32)-1))
                    rows[i]=' '.join(f)
            try:
                path.write_bytes(seal('\n'.join(rows)+'\n'))
                legacy=material_curves.read(path,fluid_path=checkpoint)
                require(legacy['id_layout']==(32,0) and legacy['curves'][0]['next_sequence']==1<<32,
                        'legacy v4 fixture missing wide exhausted counter')
                require((1<<32)-1 in [n['id'] for n in legacy['curves'][0]['nodes']],
                        'legacy v4 fixture missing UINT32_MAX')
                directory=run('legacy_v4_full_sequence',restart=checkpoint,overrides=('time/nlim=2',))
                curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
                require(-1 in [n['id'] for n in curve['nodes']] and curve['next_sequence']==1<<32,
                        'legacy v4 full-width identity/counter changed on restart')
                return dict(legacy_uint32_max=-1,next_sequence=curve['next_sequence'],geometry_restart_exact_ID_bits=True)
            finally:path.write_bytes(original)
        check('ID08_legacy_unsigned_v4_B0_full_width_restart',legacy_v4_full_sequence)
    sidecar.write_bytes(original);save();print(json.dumps(report['checks'],indent=2))
    return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
