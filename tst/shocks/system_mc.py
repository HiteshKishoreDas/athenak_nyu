"""Bounded CPU regressions, launch evidence and local performance measurements."""
import argparse
import importlib
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tst/shocks-mc'))
from common import provenance, run, uniform_changes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve(); work.mkdir(parents=True, exist_ok=False)
    report = args.report.resolve(); report.parent.mkdir(parents=True, exist_ok=True)
    results = {'final':provenance(args.executable), 'baseline':provenance(args.baseline),
               'mpi':provenance(args.mpi), 'regressions':[], 'performance':[],
               'scientific_signoff':'pending', 'remote_cluster_used':False}
    def save():
        report.write_text(json.dumps(results, indent=2)+'\n')
    # Use the existing hydro regression and its original unchanged thresholds.
    sys.path.insert(0,str(ROOT/'tst'))
    sys.path.insert(0,str(ROOT/'vis/python'))
    module = importlib.import_module('scripts.hydro.hydro_linwave')
    module._int = ['rk2']; module._recon = ['plm']; module._flux = ['hllc']
    import os
    for name,exe in [('baseline',args.baseline),('final',args.executable)]:
        directory=work/('linwave_'+name);(directory/'build/src').mkdir(parents=True)
        commands=[]
        def existing_run(input_name, options):
            # Existing test uses the removed parameter vflow; vx0 has the
            # same documented meaning (background parallel to wavevector).
            options=[x.replace('problem/vflow=', 'problem/vx0=') for x in options]
            cmd=[str(exe.resolve()),'-i',str(ROOT/'inputs'/input_name)]+options
            commands.append(cmd)
            completed=subprocess.run(cmd,cwd=directory/'build/src',capture_output=True,
                                     text=True,timeout=60)
            with (directory/'run.log').open('a') as log:
                log.write(completed.stdout+completed.stderr)
            if completed.returncode:raise RuntimeError('legacy wave run failed')
        module.athena.run=existing_run
        module.run()
        previous=Path.cwd();os.chdir(directory)
        try:passed=bool(module.analyze())
        finally:os.chdir(previous)
        results['regressions'].append({'name':name,'six_existing_wave_runs':True,
            'unchanged_thresholds_passed':passed,'commands':commands,
            'error_table':(directory/'build/src/hydro_lin_wave-errs.dat').read_text()})
        save()
        assert passed
    assert results['regressions'][0]['error_table']==results['regressions'][1]['error_table']
    results['wave_error_tables_exact_equal']=True
    # Local small-job medians include startup/output; not a production-scale budget.
    for nx in (64,256):
        measurements={}
        for mode,count,enabled in [('off',1024,False),('full',1024,True),
                                   ('empty_off',0,False),('mesh_only_empty',0,True)]:
            samples=[];runs=[]
            for trial in range(3):
                changes=uniform_changes(count=count)
                changes.update({'mesh/nx1':nx,'meshblock/nx1':nx,
                    'shock_tracking/enabled':str(enabled).lower(),
                    'shock_tracking/trace':'false','output2/dcycle':100000,
                    'output3/dt':-1})
                directory=work/f'perf_{nx}_{mode}_{trial}'
                changes['shock_tracking/directory']=str(directory/'shock_events')
                start=time.monotonic();item=run(args.executable,directory,changes)
                elapsed=time.monotonic()-start
                assert item['exit']==0
                samples.append(elapsed);runs.append(item)
            measurements[mode]={'whole_process_seconds':samples,
                'median_seconds':statistics.median(samples),'runs':runs}
        measurements['full_off_ratio']=measurements['full']['median_seconds']/measurements['off']['median_seconds']
        results['performance'].append({'mesh':[nx,4,1], 'measurements':measurements,
            'qualification':'local whole-process medians only; no accepted cluster budget'})
        save()
    # MPI singleton and actual two-rank launcher are separate evidence categories.
    changes=uniform_changes(count=128,blocks=4)
    changes.update({'shock_tracking/enabled':'true','shock_tracking/trace':'false',
                    'time/tlim':.02,'output1/single_file_per_rank':'true'})
    item=run(args.mpi,work/'mpi_singleton',changes)
    results['mpi_singleton']=item;save();assert item['exit']==0
    from common import edit_input
    directory=work/'mpi_two_ranks';directory.mkdir()
    changes['shock_tracking/directory']=str(directory/'shock_events')
    inp=directory/'input.athinput'
    inp.write_text(edit_input((ROOT/'inputs/tests/shock_mc.athinput').read_text(),changes))
    cmd=['/opt/homebrew/bin/mpiexec','-n','2',str(args.mpi.resolve()),'-i',str(inp)]
    try:
        completed=subprocess.run(cmd,cwd=directory,capture_output=True,text=True,timeout=30)
        (directory/'run.log').write_text(completed.stdout+completed.stderr)
        results['mpi_two_ranks']={'command':cmd,'exit':completed.returncode,
                                 'runtime_tested':completed.returncode==0}
    except subprocess.TimeoutExpired as exc:
        results['mpi_two_ranks']={'command':cmd,'exit':'timeout','runtime_tested':False,
                                 'diagnostic':str(exc)}
    results['status']='passed_cpu_checks' if results.get('wave_error_tables_exact_equal') else 'failed'
    save();print(results['status'])


if __name__=='__main__':
    main()
