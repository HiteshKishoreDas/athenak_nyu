"""Old paired checkpoint collision and deterministic new-basename replay evidence."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from common import ROOT, provenance, run, sha
sys.path.insert(0,str(ROOT/'vis/python'))
from shock_history import read_history


def canonical(events):
    return json.dumps(events,sort_keys=True,allow_nan=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    report=args.report.resolve();report.parent.mkdir(parents=True,exist_ok=True)
    results=dict(provenance(args.executable),runs=[])
    seed=work/'seed'
    item=run(args.executable,seed,{'time/nlim':1,'output1/dt':1,
        'tracer_seed_initial/count_per_event':64,'shock_tracking/directory':str(seed/'shock_events')})
    item['label']='seed';results['runs'].append(item)
    if item['exit']!=0:raise RuntimeError('seed process failed')
    checkpoint=sorted((seed/'rst').glob('*.rst'))[-1]
    def execute(label,basename=None):
        overrides=['time/nlim=-1','output1/dt=1']
        if basename:overrides.append('job/basename='+basename)
        command=[str(args.executable.resolve()),'-r',str(checkpoint)]+overrides
        started=time.monotonic();log=seed/(label+'.log')
        with log.open('w') as file:
            code=subprocess.run(command,cwd=seed,stdout=file,stderr=subprocess.STDOUT,timeout=30).returncode
        result={'label':label,'command':command,'exit':code,'elapsed_seconds':time.monotonic()-started,
                'log_sha256':sha(log),'directory':str(seed)}
        results['runs'].append(result)
        return result,log.read_text()
    later,text=execute('later_checkpoint');assert later['exit']==0
    before={str(p):sha(p) for p in (seed/'rst').iterdir() if p.is_file()}
    events=canonical(read_history(seed/'shock_events'))
    collision,text=execute('old_checkpoint_collision')
    assert collision['exit']!=0 and 'existing paired checkpoint' in text
    assert all(sha(p)==digest for p,digest in before.items())
    assert canonical(read_history(seed/'shock_events'))==events
    replay,text=execute('new_basename_replay','mc_replayed');assert replay['exit']==0
    assert canonical(read_history(seed/'shock_events'))==events
    results['status']='passed_immutable_pair_and_canonical_replay'
    results['checks']={'test_ids':['RST04','RST07','IO08'],'old_checkpoint':str(checkpoint),
        'old_pairs_byte_unchanged':True,'canonical_event_history_unchanged':True,
        'new_basename_replay_supported':True,
        'legacy_thermo_scope':'legacy append-only thermo may append identical rows; no claim its whole file remains unchanged'}
    report.write_text(json.dumps(results,indent=2)+'\n')
    target=report.parent/'replay_logs';target.mkdir(exist_ok=True)
    import shutil
    for p in seed.glob('*.log'):shutil.copy2(p,target/p.name)
    shutil.copy2(seed/'input.athinput',target/'input.athinput')
    print(json.dumps({'status':results['status'],'report':str(report)}))


if __name__=='__main__':main()
