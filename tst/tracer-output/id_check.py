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

LIMIT=(1<<24)-1


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
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    binary=args.executable.resolve()
    report=dict(checks=[],runs=[],runner_sha256=sha(__file__),scope='local Serial/double',
                fixture='synthetic counter/link boundaries on accepted uniform flow',scientific_approval='pending')
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:row=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:row=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(row);save()
    def require(test,message):
        if not test:raise ValueError(message)
    def run(name,changes=None,restart=None,overrides=(),reject=None):
        directory=work/name;directory.mkdir()
        options={'shock_tracking/enabled':'false','time/nlim':0,'curve1/id':255,'curve1/count':2,
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
    require(state['version']==3 and len(state['curves'][0]['nodes'])==3,'seed not compact/adapted')

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
        patch(LIMIT,'single')
        directory=run('last_id',restart=checkpoint,overrides=('time/nlim=1',))
        result=material_curves.read(directory/'material.curves.final.txt');curve=result['curves'][0]
        ids=[n['id'] for n in curve['nodes']]
        require(len(ids)==4 and max(ids)==(1<<32)-1,'last uint32 ID not allocated exactly')
        require(curve['next_sequence']==1<<24,'exhausted allocator sentinel missing')
        for n in curve['nodes']:
            require(node_curve_id(n['id'],24)==255,'wrong high byte')
        report['last_id_state']=str(directory/'material.curves.final.txt')
        return dict(node_ids=ids,last_id=max(ids),next_sequence=curve['next_sequence'])
    check('ID01_last_sequence_UINT32_MAX_and_links',last_id)
    def batch():
        patch(LIMIT,'batch')
        directory=run('batch_limit',restart=checkpoint,overrides=('time/nlim=1',))
        curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
        require(len(curve['nodes'])==3 and curve['next_sequence']==LIMIT and curve['reason']==2,
                'over-capacity split batch allocated IDs or missed quality flag')
        return dict(node_count=3,next_sequence=LIMIT,reason=2)
    check('ID02_atomic_split_batch_exhaustion',batch)
    def exhausted():
        patch(1<<24,'batch')
        directory=run('exhausted',restart=checkpoint,overrides=('time/nlim=1',))
        curve=material_curves.read(directory/'material.curves.final.txt')['curves'][0]
        require(len(curve['nodes'])==3 and curve['next_sequence']==1<<24 and curve['reason']==2,
                'exhausted counter wrapped/reused an identity')
        return dict(no_wrap_or_reuse=True)
    check('ID03_exhausted_allocator_no_wrap',exhausted)
    def inputs():
        for value in (0,256,4294967297):
            run('invalid_curve_'+str(value),changes={'curve1/id':value},reject='curve IDs must be unique in [1,255]')
        run('too_many_curves',changes={'material_curves/ncurves':256},reject='ncurves must lie in [0,255]')
        return dict(out_of_range_values=[0,256,4294967297],curve_count_limit=255)
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
        sidecar.write_bytes(original.replace(b'MATERIAL_CURVES 3 ASCII',b'MATERIAL_CURVES 2 ASCII',1))
        run('old_checkpoint',restart=checkpoint,overrides=('time/nlim=1',),reject='restart requires v3 curve8-seq24 IDs')
        if args.historical_run:
            geometry=next(args.historical_run.glob('*.curves.*.txt'))
            old=material_curves.read(geometry)
            rows,_=read_curves(args.historical_run/'curve_thermo_history')
            require(old['version']==2 and rows and all(p>1<<32 for g,p in rows),'historical readers lost 64-bit IDs')
        return dict(old_restart_rejected=True,historical_analysis=bool(args.historical_run))
    check('ID06_legacy_analysis_and_explicit_restart_rejection',legacy)
    sidecar.write_bytes(original);save();print(json.dumps(report['checks'],indent=2))
    return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
