#!/usr/bin/env python3
"""Join one verified SGS restart continuation to its diagnostic prefix.

Cumulative ledgers are never rebased. Native files stay untouched. Statistical
screening is operational evidence, not published acceptance or seed uncertainty.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

import verify
from paper_plots import statistics, write_report
from paper_test import cell_diagnostics


INVENTORIES = ('thermal-E', 'sgs-E', 'represented-E', 'force-work',
               'sgs-production-integral', 'sgs-dissipation-integral',
               'sgs-pressure-work-integral')


def join(prefix, continuation, origin):
    """Require an exact native boundary; reject reset counters and overlaps."""
    boundary = origin['time']
    if not prefix or not continuation:
        raise ValueError('empty segment')
    if any(b['time'] <= a['time'] for rows in (prefix, continuation)
           for a, b in zip(rows, rows[1:])):
        raise ValueError('segment time is not strictly increasing')
    if abs(prefix[-1]['time']-boundary) > 1e-12:
        raise ValueError('prefix does not end at pinned checkpoint')
    tolerance = verify.LIMITS['history_vs_restart_scaled_absolute']
    for key in INVENTORIES:
        if verify.scaled_diff(prefix[-1][key], origin[key]) > tolerance:
            raise ValueError('restart boundary mismatch: '+key)
    if continuation[0]['time'] < boundary:
        raise ValueError('continuation predates pinned checkpoint')
    if continuation[0]['time'] == boundary:
        for key in INVENTORIES:
            if verify.scaled_diff(continuation[0][key], origin[key]) > tolerance:
                raise ValueError('duplicate boundary mismatch: '+key)
    boundary_row = dict(prefix[-1], **origin)
    return [*prefix[:-1], boundary_row, *[row for row in continuation
                                  if row['time'] > boundary]]


def stationarity(data, result, protocol):
    if not result['complete_paper_window']:
        return dict(status='incomplete', protocol=protocol)
    metrics = {}
    for key in protocol['quantities']:
        stat = result['statistics'][key]
        early, late = result['early_late'][key]
        difference = abs(late-early)/abs(stat['mean']) if stat['mean'] else None
        drift = stat['fractional_linear_drift']
        passed = (drift is not None and difference is not None and
                  abs(drift) <= protocol['maximum_absolute_fractional_linear_drift'] and
                  difference <= protocol['maximum_early_late_relative_difference'])
        metrics[key] = dict(fractional_linear_drift=drift,
                            early_late_relative_difference=difference,
                            screen_passed=passed)
    imbalance = result['production_dissipation_relative_imbalance']
    passed = (all(row['screen_passed'] for row in metrics.values()) and
              imbalance is not None and
              imbalance <= protocol['maximum_production_dissipation_imbalance'])
    T = data['protocol']['T']
    cs = statistics([s['time'] for s in data['samples']],
                    [s['cs2_mean'] for s in data['samples']], 2*T, 10*T)
    return dict(status='pass' if passed else 'fail', protocol=protocol, metrics=metrics,
                production_dissipation_relative_imbalance=imbalance,
                sound_speed_squared=cs,
                qualification='Operational screen only; correlated blocks, one seed')


def run(args):
    prefix = json.loads(args.prefix_diagnostics.read_text())
    prefix.pop('final', None)
    if prefix.get('numerical_status') != 'pass':
        raise ValueError('prefix did not pass numerical verification')
    histpath, = args.prefix_run.glob('*/*.hydro.hst')
    hist = verify.history(histpath)
    first = json.loads((histpath.parent/'history-origin.json').read_text())
    if first['time'] < hist[0]['time']:
        matching, = [h for h in prefix['history']
                     if abs(h['time']-first['time']) < 1e-10]
        hist.insert(0, dict(matching, **first))
    prefix_results = json.loads((args.prefix_run/'results.json').read_text())
    if len(prefix_results) != 1 or prefix_results[0]['status'] != 'pass':
        raise ValueError('prefix native verification did not pass')
    passivity = json.loads((args.prefix_run/'passivity.json').read_text())
    if not (passivity['passed'] and passivity['exact_forcing_metadata'] and
            passivity['exact_rng']):
        raise ValueError('prefix fields/forcing were not validated as passive')
    results = json.loads((args.continuation/'results.json').read_text())
    if len(results) != 1 or results[0]['status'] != 'pass':
        raise ValueError('continuation has not passed native verification')
    newpath, = args.continuation.glob('*/*.hydro.hst')
    origin = json.loads((newpath.parent/'history-origin.json').read_text())
    state = verify.restart(args.checkpoint)
    cfg = state['config']
    for key, parameter in [('n', 'mesh/nx1'), ('gamma', 'hydro/gamma'),
                           ('amplitude', 'turb_driving/accel_rms'),
                           ('seed', 'turb_driving/rseed'),
                           ('sol_fraction', 'turb_driving/sol_fraction'),
                           ('T', 'turb_driving/tcorr')]:
        error = verify.scaled_diff(float(prefix['protocol'][key]), float(cfg[parameter]))
        if error > 1e-12:
            raise ValueError('prefix protocol disagrees with checkpoint: '+key)
    if prefix['protocol']['reconstruction'] != cfg['hydro/reconstruct']:
        raise ValueError('prefix reconstruction disagrees with checkpoint')
    if 'cfl' in prefix['protocol']:
        error = verify.scaled_diff(prefix['protocol']['cfl'], float(cfg['time/cfl_number']))
        if error > 1e-12:
            raise ValueError('prefix CFL disagrees with checkpoint')
    for key, value in [('time', state['time']), ('thermal-E', state['thermal']),
                       ('sgs-E', state['totals'][-1]),
                       ('represented-E', sum(state['totals'][v] for v in (4, -1)))]:
        if verify.scaled_diff(origin[key], value) > 1e-12:
            raise ValueError('continuation origin disagrees with checkpoint: '+key)
    for label, parameter in [('force-work', 'turb_driving/integrated_work'),
                             ('sgs-production-integral', 'sgs/integrated_production'),
                             ('sgs-dissipation-integral', 'sgs/integrated_dissipation'),
                             ('sgs-pressure-work-integral',
                              'sgs/integrated_pressure_work')]:
        if verify.scaled_diff(origin[label], float(cfg[parameter])) > 1e-12:
            raise ValueError('counter origin disagrees with checkpoint: '+label)
    samples = [s for s in prefix['samples'] if s['time'] >= first['time']-1e-10]
    if abs(samples[0]['time']-first['time']) > 1e-10:
        raise ValueError('prefix diagnostic sample missing at heating origin')
    if abs(samples[-1]['time']-origin['time']) > 1e-12:
        raise ValueError('prefix samples do not end at pinned checkpoint')
    if samples[0]['time'] > 2*float(cfg['turb_driving/tcorr']):
        raise ValueError('missing diagnostics before 2T')
    for key, value in [('thermal', state['thermal']), ('K', state['totals'][-1]),
                       ('represented', sum(state['totals'][v] for v in (4, -1)))]:
        if verify.scaled_diff(samples[-1][key], value) > 1e-12:
            raise ValueError('prefix sample does not match pinned checkpoint: '+key)
    history = join(hist, verify.history(newpath), origin)
    n = int(cfg['mesh/nx1'])
    gamma = float(cfg['hydro/gamma'])
    T = float(cfg['turb_driving/tcorr'])
    checkpoints = sorted((newpath.parent/'rst').glob('*.rst'))
    del state
    for path in checkpoints:
        state = verify.restart(path)
        for key in cfg:
            if (key.startswith(('mesh/', 'meshblock/', 'hydro/', 'sgs/', 'turb_driving/'))
                    or key in ('time/integrator', 'time/evolution', 'time/cfl_number')):
                if '/integrated_' not in key and key != 'sgs/fallbacks':
                    if state['config'][key] != cfg[key]:
                        raise ValueError('changed restart physics: '+key)
        if state['time'] <= origin['time']:
            continue
        values = np.array(list(state['field'].values()))
        sample, mach = cell_diagnostics(values, gamma, 1/n, float(cfg['sgs/ceps']))
        sample.update(time=state['time'], kinetic=state['kinetic'],
                      thermal=state['thermal'],
                      represented=sum(state['totals'][v] for v in (4, -1)))
        samples.append(sample)
        del state
    protocol = dict(prefix['protocol'])
    protocol.update(end_T=samples[-1]['time']/T, T=T, n=n, gamma=gamma,
                    amplitude=float(cfg['turb_driving/accel_rms']),
                    seed=int(cfg['turb_driving/rseed']),
                    reconstruction=cfg['hydro/reconstruct'],
                    cfl=float(cfg['time/cfl_number']), actual_restart_config=cfg)
    protocol['departures'] = [d for d in protocol['departures']
                              if not d.startswith('No resolved pressure-work')]
    protocol['departures'].append(
        'Discrete pressure-work ledger includes spatial/RK effects')
    data = dict(protocol=protocol, samples=samples, history=history,
                final=dict(density=values[:, 0].tolist(), sgs_mach=mach.tolist()),
                numerical_status='pass')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'diagnostics.json').write_text(json.dumps(data, indent=2)+'\n')
    result = write_report(args.output)
    screen = stationarity(data, result, json.loads(args.screen.read_text()))
    (args.output/'stationarity.json').write_text(json.dumps(screen, indent=2)+'\n')
    assembled = args.output/'assembled'
    assembled.mkdir()
    keys = list(hist[0])
    header = '# Assembled verified native history; no counter rebasing\n# '+ ' '.join(
        f'[{i}]={key}' for i, key in enumerate(keys))+'\n'
    (assembled/'assembled.hydro.hst').write_text(header + ''.join(
        ' '.join(f'{row[key]:.17e}' for key in keys)+'\n' for row in history))
    (assembled/'rst').mkdir()
    (assembled/'rst/initial.rst').symlink_to(args.checkpoint.resolve())
    files = [args.prefix_diagnostics, histpath, newpath, args.checkpoint, args.screen,
             args.prefix_run/'passivity.json', args.prefix_run/'results.json',
             histpath.parent/'history-origin.json', newpath.parent/'history-origin.json',
             args.continuation/'results.json',
             args.continuation/'manifest.json', Path(__file__).resolve(),
             Path(__file__).with_name('paper_plots.py'),
             Path(__file__).with_name('paper_test.py'),
             Path(__file__).with_name('verify.py')]
    (args.output/'provenance.json').write_text(json.dumps(dict(
        files=[dict(path=str(p.resolve()), sha256=verify.sha(p)) for p in files],
        checkpoints=[dict(path=str(p.resolve()), sha256=verify.sha(p))
                     for p in checkpoints], boundary=origin,
        prefix_passivity=passivity, rebased_counters=False,
        analysis_source_commit=subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=verify.SOURCE, text=True).strip()), indent=2)+'\n')
    print(json.dumps(dict(output=str(args.output), screen=screen['status'],
                          complete=result['complete_paper_window'])), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('prefix-diagnostics', 'prefix-run', 'continuation', 'checkpoint',
                 'screen', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    run(p.parse_args())


if __name__ == '__main__':
    main()
