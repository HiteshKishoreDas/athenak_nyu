#!/usr/bin/env python3
"""Signed periodic-box heating budget from applied counters and native history.

The residual includes spatial/RK effects. It is not an LLF-only diagnostic.
Fractions use compression-corrected heating, not forcing power or total heat.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from paper_plots import window
from verify import history, parameters

PRESSURE = 'sgs-pressure-work-integral'


def budget(hist, start, stop):
    times = [row['time'] for row in hist]
    series = {}
    for label, key in (('thermal', 'thermal-E'), ('compression', PRESSURE),
                       ('sgs', 'sgs-dissipation-integral'),
                       ('production', 'sgs-production-integral'),
                       ('K', 'sgs-E'), ('forcing', 'force-work'),
                       ('represented', 'represented-E')):
        t, y = window(times, [row[key] for row in hist], start, stop)
        series[label] = y-y[0]
    series['numerical'] = series['thermal']-series['compression']-series['sgs']
    series['corrected'] = series['thermal']-series['compression']
    sgs, numerical = series['sgs'][-1], series['numerical'][-1]
    denominator = sgs+numerical
    summary = {key: float(y[-1]) for key, y in series.items()}
    summary.update(start=float(start), stop=float(stop),
                   sgs_fraction=float(sgs/denominator) if denominator > 0 else None,
                   numerical_fraction=float(numerical/denominator)
                   if denominator > 0 else None,
                   forcing_budget_residual=float(
                       series['represented'][-1]-series['forcing'][-1]),
                   sgs_budget_residual=float(
                       series['K'][-1]-series['production'][-1]+sgs))
    return t, series, summary


def self_check():
    rows = []
    for t in (0, .5, 1):
        rows.append({'time': t, 'thermal-E': 100+2.5*t, PRESSURE: -.5*t,
                     'sgs-dissipation-integral': 2*t,
                     'sgs-production-integral': 3*t, 'sgs-E': 10+t,
                     'force-work': 4*t, 'represented-E': 110+4*t})
    t, series, result = budget(rows, .25, .75)
    assert np.allclose(t, [.25, .5, .75])
    assert np.isclose(result['sgs'], 1)
    assert np.isclose(result['compression'], -.25)
    assert np.isclose(result['numerical'], .5)
    assert np.isclose(result['sgs_fraction'], 2/3)
    assert result['forcing_budget_residual'] == result['sgs_budget_residual'] == 0
    for rate in (-.5, -1.):
        for row in rows:
            row['thermal-E'] = 100+rate*row['time']
        _, _, result = budget(rows, .25, .75)
        assert result['sgs_fraction'] is result['numerical_fraction'] is None
    print('Signed heating budget, clipped endpoints and invalid fractions: PASS')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--start-T', type=float, default=2.)
    parser.add_argument('--end-T', type=float, default=4.)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_check()
        return
    if not args.runs or args.output is None:
        parser.error('--runs and --output are required')
    args.output.mkdir(parents=True, exist_ok=False)
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
    summaries = []
    for run in args.runs:
        paths = list(run.glob('*/*.hydro.hst'))
        if len(paths) != 1:
            raise ValueError(f'need exactly one native hydro history in {run}')
        hist = history(paths[0])
        origin = paths[0].parent/'history-origin.json'
        if origin.exists():
            row = json.loads(origin.read_text())
            if row['time'] > hist[0]['time']:
                raise ValueError('restart origin is later than first stored history')
            if row['time'] < hist[0]['time']:
                hist.insert(0, row)
        checkpoints = sorted((paths[0].parent/'rst').glob('*.rst'))
        with checkpoints[0].open('rb') as handle:
            header = handle.read(1024*1024)
        if b'<par_end>\n' not in header:
            raise ValueError('native parameter header exceeds reader cap')
        cfg = parameters(header.split(b'<par_end>\n', 1)[0].decode())
        turnover = float(cfg['turb_driving/tcorr'])
        t, series, result = budget(hist, args.start_T*turnover, args.end_T*turnover)
        label = f'{cfg["mesh/nx1"]}³ {cfg["hydro/reconstruct"].upper()}'
        label += f' CFL {cfg["time/cfl_number"]}'
        result.update(label=label, run=str(run.resolve()),
                      window_T=[args.start_T, args.end_T])
        summaries.append(result)
        line, = axes[0, 0].plot(t/turnover, series['sgs'], label=label)
        color = line.get_color()
        axes[0, 0].plot(t/turnover, series['numerical'], '--', color=color)
        axes[0, 1].plot(t/turnover, series['compression'], color=color, label=label)
        increments = np.diff(series['corrected'])
        fraction = np.divide(np.diff(series['sgs']), increments,
                             out=np.full_like(increments, np.nan), where=increments > 0)
        axes[1, 0].plot((t[:-1]+t[1:])/(2*turnover), fraction, color=color)
        axes[1, 1].plot(t/turnover, series['K'], color=color)
        axes[2, 0].plot(t/turnover, series['production'], color=color)
        axes[2, 0].plot(t/turnover, series['sgs'], '--', color=color)
    x = np.arange(len(summaries))
    fractions = [np.nan if row['sgs_fraction'] is None else row['sgs_fraction']
                 for row in summaries]
    numerical = [np.nan if row['numerical_fraction'] is None
                 else row['numerical_fraction'] for row in summaries]
    axes[2, 1].bar(x-.18, fractions, width=.36, label='SGS')
    axes[2, 1].bar(x+.18, numerical, width=.36, label='Numerical residual')
    axes[2, 1].set_xticks(x, [row['label'] for row in summaries], rotation=20)
    axes[2, 1].axhline(.75, color='black', ls=':', label='Paper SGS ≈0.75')
    axes[1, 0].axhline(.75, color='black', ls=':')
    titles = ['Cumulative heating: solid SGS, dashed numerical',
              'Signed resolved compression heating',
              'SGS share in each stored time interval', 'Change in SGS energy',
              'Cumulative SGS transfer: solid production, dashed dissipation',
              'Whole-window heating fractions']
    for index, (ax, title) in enumerate(zip(axes.flat, titles)):
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=.25)
        if index != 5:
            ax.set_xlabel('Time / T')
        ax.set_ylabel('Fraction' if index in (2, 5) else 'Energy (code units)')
    axes[0, 0].legend(fontsize=8)
    axes[2, 1].legend(fontsize=8)
    window_label = ('2T–10T averages' if args.start_T == 2 and args.end_T == 10
                    else 'short-window pilots')
    fig.suptitle('Compression-corrected heating budget — '+window_label+'\n'
                 'Numerical residual includes spatial/RK effects; signed values retained',
                 fontsize=13)
    fig.savefig(args.output/'heating-comparison.png', dpi=170)
    fig.savefig(args.output/'heating-comparison.pdf')
    plt.close(fig)
    report = dict(window_T=[args.start_T, args.end_T], cases=summaries,
                  definitions=__doc__, reference='https://arxiv.org/html/1010.4492#S4.SS2',
                  caveats=['One seed per case; stationarity requires separate trend assessment',
                           'Pressure work uses the discrete reconstructed momentum flux',
                           'Residual includes spatial/RK and shock-associated heating',
                           'Forcing equivalence to paper remains unproven',
                           'Endpoints use linear interpolation of native history',
                           'Interval fractions require positive corrected heating'])
    (args.output/'heating-summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
