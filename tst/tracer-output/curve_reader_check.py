"""Check the curve plotting script against existing text/CSV/THP simulation data.

No simulations are launched. Use the existing scientific Python environment.
--runs is the output of tst/tracer-output/check.py; --output must be new.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'vis/python'))
from tracer_history import read_run
from plot_curves import select_curve,curve_points,plot_curve


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--three-d-run',type=Path)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    runs=args.runs.resolve()
    report=dict(checks=[],plots=[],simulation_runs=0,inputs=str(runs),
                sources={str(p.relative_to(ROOT)):sha(p) for p in
                    [ROOT/'vis/python/plot_curves.py',ROOT/'vis/python/tracer_history.py',Path(__file__)]})
    def save():
        (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def check(name,action):
        try:result=dict(id=name,status='pass',**(action() or {}))
        except Exception as exc:result=dict(id=name,status='fail',error=repr(exc),traceback=traceback.format_exc())
        report['checks'].append(result);save()
    def render(selected,name):
        result=plot_curve(selected,out/(name+'.png'))
        blob=(out/(name+'.png')).read_bytes()
        assert blob[:8]==b'\x89PNG\r\n\x1a\n' and len(blob)>10000
        width,height=struct.unpack('>II',blob[16:24]);assert width>=1000 and height>=700
        report['plots'].append(result);save();return result
    def rejected(action,needle):
        try:action()
        except ValueError as exc:
            assert needle in str(exc),str(exc)
        else:raise AssertionError('invalid input was accepted')

    old=read_run(runs/'planar_control',include_particles=False)
    new=read_run(runs/'planar',include_particles=False)
    a,b=select_curve(old),select_curve(new)
    def parity():
        assert a['curve_id']==b['curve_id']==7
        assert a['history']==b['history'] and len(b['history'])==108
        assert a['profile']==b['profile'] and len(b['profile'])==16
        x=render(a,'control');y=render(b,'scheduled')
        assert {k:v for k,v in x.items() if k!='output'}=={k:v for k,v in y.items() if k!='output'}
        return dict(node_samples=108,profile_nodes=16,control_scheduled_equal=True)
    check('CR01_existing_outputs_scheduled_and_control',parity)

    def newborn():
        data=read_run(runs/'amr_shock',include_particles=False)
        events=[e for key in ('curve_events','curve_suspended') for e in data[key]
                if e['birth']>0 and e['flags']&1]
        event=min(events,key=lambda e:(e['birth'],e['node_id']))
        selected=select_curve(data,node_id=event['node_id'])
        assert selected['history'][0]['time']==event['birth']
        assert all(o['time']>=o['birth']==event['birth'] for o in selected['history'])
        result=render(selected,'amr-newborn')
        assert result['node_birth']==result['first_sample_time']
        return result
    check('CR02_AMR_birth_no_parent_history',newborn)

    def large():
        directory=runs/'large_ids';node=((255<<24)|1)-(1<<32)
        command=[sys.executable,str(ROOT/'vis/python/plot_curves.py'),str(directory),
                 '--node-id',str(node),'--output',str(out/'large-id.png')]
        result=subprocess.run(command,capture_output=True,text=True,timeout=60)
        (out/'large-id-command.json').write_text(json.dumps(dict(command=command,exit=result.returncode,
            stdout=result.stdout,stderr=result.stderr),indent=2)+'\n')
        assert result.returncode==0,result.stderr
        row=json.loads(result.stdout);assert row['node_id']==node and row['curve_id']==255
        report['plots'].append(row)
        return dict(exact_id=row['node_id'],cli_exit=0)
    check('CR03_signed32_CLI',large)

    def geometry_only():
        directory=out/'geometry-only';directory.mkdir()
        for path in sorted((runs/'planar').glob('*.curves.*.txt'))[:3]:
            shutil.copy2(path,directory/path.name)
        data=read_run(directory,include_particles=False)
        selected=select_curve(data);assert selected['frames'] and not selected['history']
        render(selected,'geometry-only');return dict(optional_thermo_and_shocks_absent=True)
    check('CR04_geometry_only',geometry_only)

    def thermo_only():
        directory=out/'thermo-only/curve_thermo_history';directory.mkdir(parents=True)
        for path in sorted((runs/'planar/curve_thermo_history').glob('*.thp'))[:3]:
            shutil.copy2(path,directory/path.name)
        selected=select_curve(read_run(directory.parent,include_particles=False))
        assert not selected['frames'] and len(selected['history'])==3
        render(selected,'thermo-only');return dict(no_invented_geometry=True)
    check('CR05_thermo_only',thermo_only)

    def particle_skip():
        directory=out/'geometry-only';folder=directory/'prtcl_thermo_history';folder.mkdir()
        (folder/'bad.thp').write_bytes(b'invalid unrelated particle stream')
        assert read_run(directory,include_particles=False)['geometry']
        rejected(lambda:read_run(directory),'truncated thermo file header')
        full=read_run(runs/'planar')
        assert len(full['particles'])==13824
        return dict(curve_mode_skips_particles=True,default_still_loads_and_validates_particles=True)
    check('CR06_particle_IO_opt_out',particle_skip)

    def winding():
        frame=dict(domain=[0,1,0,2,0,3],curve=dict(closed=True,winding=[1,-1,1],
            nodes=[dict(x=[.1,.2,.3]),dict(x=[.3,.7,1]),dict(x=[.8,1.2,2])]))
        points=curve_points(frame)
        assert points[-1]==[1.1,-1.8,3.3] and points[-1]!=points[0]
        frame['curve']['closed']=False
        assert len(curve_points(frame))==3
        return dict(winding_aware_3d_closure=True)
    check('CR07_winding_closure',winding)

    def invalid():
        rejected(lambda:select_curve(new,curve_id=0),'curve ID')
        rejected(lambda:select_curve(new,curve_id=8),'not found')
        rejected(lambda:select_curve(new,node_id=2**64),'valid declared curve node ID')
        rejected(lambda:select_curve(new,curve_id=7,node_id=(8<<24)|1),'different curve')
        rejected(lambda:select_curve(new,node_id=(7<<24)|999999),'not found')
        rejected(lambda:select_curve(new,field='not_saved'),'available fields')
        many=copy.deepcopy(new)
        node=copy.deepcopy(next(iter(many['curve_nodes'].values())))
        node.update(curve_id=8,node_id=(8<<24)|1)
        many['curve_nodes'][0,node['node_id']]=node
        rejected(lambda:select_curve(many),'choose --curve-id')
        return dict(invalid_or_ambiguous_selections_rejected=7)
    check('CR08_selection_guards',invalid)

    def overwrite():
        before=sha(out/'scheduled.png')
        rejected(lambda:plot_curve(b,out/'scheduled.png'),'already exists')
        assert sha(out/'scheduled.png')==before
        return dict(existing_plot_preserved=True)
    check('CR09_output_preservation',overwrite)

    def malformed():
        directory=out/'corrupt';directory.mkdir()
        path=next((runs/'planar').glob('*.curves.*.txt'))
        (directory/path.name).write_bytes(path.read_bytes()[:-8])
        rejected(lambda:read_run(directory,include_particles=False),'truncated')
        return dict(corrupt_geometry_rejected=True)
    check('CR10_corrupt_curve_input',malformed)

    if args.three_d_run:
        def three_d():
            selected=select_curve(read_run(args.three_d_run,include_particles=False))
            result=render(selected,'three-d')
            assert result['geometry_3d']
            return dict(source=str(args.three_d_run.resolve()),stored_3d_geometry=True)
        check('CR11_existing_3d_run',three_d)
    save();print(json.dumps(report['checks'],indent=2))
    return int(any(c['status']!='pass' for c in report['checks']))


if __name__=='__main__':raise SystemExit(main())
