"""Run the pinned MC baseline; never asserts deterministic physical particle paths."""
import argparse
import json
from pathlib import Path
import shutil
import time

from common import ROOT, moment_oracle, run, sha, thermo, uniform_changes

SEEDS = (11, 29, 47, 83, 127, 211, 307, 401)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve(); work.mkdir(parents=True, exist_ok=False)
    report = args.report.resolve(); report.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    result = {'scope':'baseline MC tracer qualification only, shock history unimplemented',
              'executable':str(args.executable.resolve()), 'executable_sha256':sha(args.executable),
              'input_sha256':sha(ROOT/'inputs/tests/shock_mc.athinput'),
              'predeclared_seeds':SEEDS, 'runs':[], 'uniform_moments':[], 'limitations':[],
              'scientific_signoff':'pending'}
    def save():
        result['elapsed_seconds'] = time.monotonic()-started
        report.write_text(json.dumps(result, indent=2)+'\n')
    def execute(label, changes=None, restart=None, overrides=()):
        if time.monotonic()-started > 300:
            raise TimeoutError('baseline300s budget exhausted')
        item = run(args.executable, work/label, changes, restart, overrides)
        item['label'] = label; result['runs'].append(item); save()
        target = report.parent/'runs'/label; target.mkdir(parents=True, exist_ok=True)
        for name in ('input.athinput','run.log'):
            path = work/label/name
            if path.exists():shutil.copy2(path,target/name)
        (target/'manifest.json').write_text(json.dumps(item,indent=2)+'\n')
        if item['exit'] != 0:raise RuntimeError(label+': '+str(item['exit']))
        return work/label
    try:
        for seed in SEEDS:
            directory = execute('uniform_seed'+str(seed), uniform_changes(seed))
            observations,times = thermo(directory)
            measurement = moment_oracle(observations,times)
            measurement['seed'] = seed
            measurement['legacy_identical_terminal_replay_rows'] = observations.identical_replay_rows
            result['uniform_moments'].append(measurement); save()
        changes = uniform_changes(count=1024,blocks=4)
        directory = execute('uniform_fourblocks',changes)
        observations,times = thermo(directory)
        first,last = min(times),max(times)
        initial = {tag:o for (cycle,tag),o in observations.items() if cycle==first}
        final = {tag:o for (cycle,tag),o in observations.items() if cycle==last}
        result['local_migration'] = {'count':len(final), 'identity_preserved':initial.keys()==final.keys(),
            'changed_gid':sum(initial[tag]['gid'] != final[tag]['gid'] for tag in initial),
            'moments':moment_oracle(observations,times)}
        seedchanges=uniform_changes(count=1024,blocks=4)
        seedchanges.update({'time/nlim':17,'output1/dt':1})
        directory=execute('checkpoint',seedchanges)
        checkpoint=sorted((directory/'rst').glob('*.rst'))[-1]
        continued=execute('continuation',restart=checkpoint,
            overrides=('time/nlim=-1','job/basename=mc_continued','output1/dt=-1'))
        restored,restored_times=thermo(continued)
        # The same four-block fluid/tracer run uses exactly the same fixed seed.
        result['restart']={'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),
            'first_restored_cycle':min(restored_times),
            'exact_observations_equal':all(observations[key]==o for key,o in restored.items()),
            'observations_compared':len(restored)}
        zero=execute('zero',uniform_changes(count=0))
        empty,zero_times=thermo(zero)
        result['zero']={'observations':len(empty),'passed':len(empty)==0}
        passed=all(x['mean_pass'] and x['variance_pass'] for x in result['uniform_moments'])
        passed=passed and result['local_migration']['identity_preserved']
        passed=passed and result['local_migration']['changed_gid']>0
        passed=passed and result['restart']['exact_observations_equal'] and result['zero']['passed']
        result['status']='passed_selected_baseline_checks' if passed else 'baseline_acceptance_failure'
    except Exception as exc:
        result['status']='failed';result['error']=repr(exc)
    save();print(json.dumps({'status':result['status'],'report':str(report)}))


if __name__ == '__main__':
    main()
