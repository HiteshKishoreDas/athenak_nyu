"""Run small local MC tests and export measured plots; no scientific pass thresholds.

Requires installed NumPy/Matplotlib and AthenaK's existing binary reader (h5py).
Run from any directory; --workdir must be new unless --reuse is explicit.
"""
import argparse
import csv
import html
import json
import math
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys

from common import ROOT, provenance, run, sha, thermo, uniform_changes, moment_oracle
from integration import canonical, traces
from qualification import rh_states

sys.path.insert(0, str(ROOT / 'vis/python'))
from shock_history import read_history, read_suspended
from bin_convert import read_binary
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reuse', action='store_true', help='reuse completed test reports')
    args = parser.parse_args()
    work, output = args.workdir.resolve(), args.output.resolve()
    binary = args.executable.resolve()
    if not args.reuse:
        work.mkdir(parents=True, exist_ok=False)
    output.mkdir(parents=True, exist_ok=False)
    commands = []
    for suite in ('integration', 'qualification'):
        report = work / (suite + '-report.json')
        command = [sys.executable, str(ROOT / 'tst/shocks-mc' / (suite + '.py')),
                   '--executable', str(binary), '--workdir', str(work / suite),
                   '--report', str(report)]
        if not args.reuse:
            subprocess.run(command, check=True, timeout=330)
        commands.append(command)
    engineering = json.loads((work / 'integration-report.json').read_text())
    physical = json.loads((work / 'qualification-report.json').read_text())
    assert engineering['status'] == 'passed_supported_engineering_checks'
    assert len(physical['runs']) == 21 and all(r['exit'] == 0 for r in physical['runs'])
    assert all('error' not in m for m in physical['measurements'])
    identity = provenance(binary)
    for report in (engineering, physical):
        assert report['executable_sha256'] == identity['executable_sha256']
        assert report['compiled_source_fingerprint'] == identity['compiled_source_fingerprint']
    for name in ('integration-report.json', 'qualification-report.json'):
        shutil.copy2(work / name, output / name)
    shutil.copytree(work / 'runs', output / 'runs')
    figures, measurements = [], {m['label']: m for m in physical['measurements']}
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.facecolor': 'white'})
    pdf = PdfPages(output / 'all-plots.pdf')

    def save(fig, name, caption):
        fig.tight_layout()
        fig.savefig(output / (name + '.png'), dpi=170)
        pdf.savefig(fig)
        plt.close(fig)
        figures.append({'name': name, 'caption': caption})

    def execute(label, changes):
        changes['shock_tracking/directory'] = str(work / label / 'shock_events')
        item = run(binary, work / label, changes)
        assert item['exit'] == 0, item
        target = output / 'runs' / label
        target.mkdir()
        for name in ('input.athinput', 'run.log'):
            shutil.copy2(work / label / name, target / name)
        return item

    # Full native grid snapshots: never expand a 1D slice into a fake 2D field.
    snapshot_run = execute('visual_snapshots', {
        'output4/file_type': 'bin', 'output4/variable': 'hydro_w', 'output4/dt': .025})
    snapshots = [read_binary(p) for p in sorted((work / 'visual_snapshots/bin').glob('*.bin'))]
    assert len(snapshots) >= 10
    assert all(d['n_mbs'] == 1 and d['mb_data']['dens'][0].shape == (1, 4, 128)
               for d in snapshots)
    data = output / 'data'
    data.mkdir()
    shutil.copytree(work / 'visual_snapshots/bin', data / 'native-2d-snapshots')
    observed = traces(work / 'visual_snapshots/shock_events')
    selected = [snapshots[i] for i in (0, 3, 6, len(snapshots) - 1)]
    x = (np.arange(128) + .5) / 128
    for field, title in [('dens', 'Density'), ('pressure', 'Pressure'),
                         ('temperature', 'Temperature proxy P / density')]:
        arrays = [(d['mb_data']['dens'][0][0] if field == 'dens' else
                   .4 * d['mb_data']['eint'][0][0] if field == 'pressure' else
                   .4 * d['mb_data']['eint'][0][0] / d['mb_data']['dens'][0][0])
                  for d in selected]
        fig, axes = plt.subplots(4, 1, figsize=(11, 8), sharex=True, sharey=True)
        for ax, d, a in zip(axes, selected, arrays):
            im = ax.imshow(a, origin='lower', extent=(0, 1, 0, .1), aspect='auto',
                           vmin=min(v.min() for v in arrays),
                           vmax=max(v.max() for v in arrays), cmap='viridis')
            hits = [o for (g, p), o in observed.items()
                    if g == d['cycle'] and o['inside']]
            ax.scatter([o['x'][0] for o in hits], [o['x'][1] for o in hits],
                       s=12, color='white', edgecolors='black', linewidths=.4)
            ax.axvline(.3 + .5 * d['time'], color='tomato', ls='--')
            ax.set_ylabel('y')
            ax.set_title(f"t={d['time']:.4f}, generation={d['cycle']}", fontsize=10)
            fig.colorbar(im, ax=ax, pad=.01, label=title)
        axes[0].set_title(title + ': actual 128 × 4 grid; dashed = analytic front')
        axes[-1].set_xlabel('x (y stretched to display all four cell rows)')
        save(fig, '2d_' + field, 'Native 2D fields. White dots are sampled particles '
             'inside the unreduced shock region, not a reconstruction of the full mesh mask.')
    last = snapshots[-1]
    rho = last['mb_data']['dens'][0][0, 0]
    pressure = .4 * last['mb_data']['eint'][0][0, 0]
    velocity = last['mb_data']['velx'][0][0, 0]
    states = rh_states()
    front = .3 + .5 * last['time']
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, values, left, right, name in zip(axes.flat,
            (rho, pressure, velocity, pressure / rho),
            (states['problem/dl'], states['problem/pl'], states['problem/ul'],
             states['problem/pl'] / states['problem/dl']),
            (1, 1, states['problem/ur'], 1), ('Density', 'Pressure', 'Velocity x', 'P / density')):
        ax.plot(x, values, label='numerical cell centers')
        ax.plot([0, front, front, 1], [left, left, right, right], 'k--', label='RH plateaus')
        ax.set(xlabel='x', ylabel=name)
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle(f'Independent Mach 2 Rankine–Hugoniot reference at t={last["time"]:.3f}')
    save(fig, 'rh_profiles', 'Analytic discontinuity and plateaus versus resolved fluid profiles; '
         'finite shock thickness is visible. No strength acceptance tolerance is inferred.')
    fig, ax = plt.subplots(figsize=(10, 5))
    tt = [d['time'] for d in snapshots]
    rr = [d['mb_data']['dens'][0][0, 0] for d in snapshots]
    ax.pcolormesh(x, tt, rr, shading='nearest', cmap='viridis')
    ax.plot(.3 + .5 * np.array(tt), tt, 'r--', label='known Eulerian front')
    hits = [o for o in observed.values() if o['inside']]
    ax.scatter([o['x'][0] for o in hits], [o['time'] for o in hits], s=1, c='white', alpha=.25)
    ax.set(xlabel='x', ylabel='code time', title='Density evolution and sampled shock-region membership')
    ax.legend()
    save(fig, 'shock_space_time', 'White dots are synchronized sampled mask hits; '
         'the red line uses the prescribed shock speed, not a fitted particle crossing time.')
    midpoint = (states['problem/dl'] + 1) / 2
    fitted = []
    for row in rr:
        indices = np.flatnonzero((row[:-1] >= midpoint) & (row[1:] < midpoint))
        assert len(indices) == 1, 'density midpoint is ambiguous'
        i = indices[0]
        fitted.append(x[i] + (midpoint - row[i]) / (row[i + 1] - row[i]) / 128)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(tt, fitted, 'o-', label='interpolated density midpoint')
    axes[0].plot(tt, .3 + .5 * np.array(tt), 'k--', label='analytic front')
    axes[0].legend(fontsize=9)
    axes[0].set(xlabel='code time', ylabel='Eulerian front x')
    error_cells = (np.array(fitted) - (.3 + .5 * np.array(tt))) * 128
    axes[1].plot(tt, error_cells, 'o-')
    axes[1].axhline(0, color='black', ls='--')
    axes[1].set(xlabel='code time', ylabel='midpoint error / cell width')
    save(fig, 'front_position', 'Density-midpoint interpolation is an Eulerian diagnostic, '
         'not an encounter timing estimator. Initial cell placement contributes an offset.')
    with (data / 'front-position.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'density_midpoint_x', 'analytic_front_x', 'error_cells'])
        writer.writerows(zip(tt, fitted, .3 + .5 * np.array(tt), error_cells))

    # Sampled MC parent histories, with known failing parents explicitly selected.
    histories = {}
    labels = ('seed11', 'mach1.2', 'mach1.5', 'mach3', 'mach5', 'boost2',
              'stationary', 'reflected', 'nx64', 'nx256')
    for label in labels:
        source = work / 'qualification' / label / 'shock_events'
        events, suspended = read_history(source), read_suspended(source)
        rows = traces(source)
        extras = sorted(p for p, count in measurements[label]['particles_with_multiple_episodes'].items())
        ids = sorted({0, 1, 2} | {int(p) for p in extras} |
                     {e['particle_id'] for e in suspended})
        histories[label] = (rows, events, suspended, ids)
        speed = 2.5 if label == 'boost2' else 0 if label == 'stationary' else -.5 if label == 'reflected' else .5
        origin = .7 if label == 'reflected' else .3
        fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
        exported = []
        for pid in ids:
            obs = [o for (g, p), o in sorted(rows.items()) if p == pid]
            t, xx = [o['time'] for o in obs], [o['x'][0] for o in obs]
            line, = axes[0].step(t, xx, where='post', lw=1, label='ID ' + str(pid))
            inside = [o for o in obs if o['inside']]
            axes[0].scatter([o['time'] for o in inside], [o['x'][0] for o in inside],
                            s=9, color=line.get_color())
            axes[1].step(t, np.array(xx) - (origin + speed * np.array(t)),
                         where='post', color=line.get_color(), lw=1)
            exported.extend([dict(particle_id=pid, **o) for o in obs])
        end = max(o['time'] for o in rows.values())
        axes[0].plot([0, end], [origin, origin + speed * end], 'k--', label='analytic front')
        axes[0].legend(ncol=4, fontsize=9)
        axes[0].set_ylabel('sampled parent-cell x')
        axes[0].set_title(label + ': MC jumps; dots = inside unreduced shock region')
        axes[1].axhline(0, color='black', ls='--')
        axes[1].set(xlabel='code time', ylabel='parent x − analytic front x')
        save(fig, label + '_parents', 'MC positions are cell centers. Lines display discrete '
             'sample sequences, not continuous material trajectories or exact crossing times.')
        with (data / (label + '-selected.csv')).open('w') as f:
            writer = csv.writer(f)
            writer.writerow(['particle_id', 'generation', 'time', 'x', 'y', 'valid', 'inside',
                             'phase', 'sequence', 'rho', 'pressure', 'mach_pressure', 'mach_temperature'])
            for o in exported:
                writer.writerow([o['particle_id'], o['generation'], o['time'], *o['x'][:2],
                    int(o['valid']), int(o['inside']), o['phase'], o['sequence'], o['rho'],
                    o['pressure'], o['mach_pressure'], o['mach_temperature']])
    for label in ('seed11', 'mach5', 'boost2'):
        rows, events, suspended, ids = histories[label]
        fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
        for pid in ids[:6]:
            obs = [o for (g, p), o in sorted(rows.items()) if p == pid]
            t = [o['time'] for o in obs]
            axes[0].step(t, [o['phase'] for o in obs], where='post', lw=1, label='ID ' + str(pid))
            for ax, key in zip(axes[1:], ('mach_pressure', 'mach_temperature')):
                hits = [o for o in obs if o['inside'] and math.isfinite(o[key])]
                ax.scatter([o['time'] for o in hits], [o[key] for o in hits], s=9)
        axes[0].set_yticks([0, 1, 2], ['outside', 'active', 'pending exit'])
        axes[0].legend(ncol=3, fontsize=9)
        axes[0].set_title(label + ': encounter state and all inside strength samples')
        expected = measurements[label]['expected_input_mach']
        for ax, text in zip(axes[1:], ('Pressure Mach', 'Temperature Mach')):
            ax.axhline(expected, color='black', ls='--', label='input Mach')
            ax.set_ylabel(text)
        axes[-1].set_xlabel('code time')
        save(fig, label + '_strength_time', 'Strength is based on fixed detector endpoints; '
             'these are not certified upstream/downstream plateaus. Gaps mean no inside sample.')
    names = [m['label'] for m in physical['measurements']]
    fig, ax = plt.subplots(figsize=(13, 5))
    complete = [measurements[n]['completed_events'] for n in names]
    partial = [measurements[n]['partial_events'] for n in names]
    ax.bar(names, complete, label='complete episodes', color='#337eb8')
    ax.bar(names, partial, bottom=complete, label='partial episodes', color='#ee9955')
    ax.plot(names, [measurements[n]['emitted_particles'] for n in names], 'ko', ms=4,
            label='unique particles with events')
    ax.plot(names, [measurements[n]['terminal_suspended'] for n in names], 'rx',
            label='suspended at cutoff')
    ax.axhline(128, color='gray', ls='--', label='128 seeded particles')
    ax.tick_params(axis='x', rotation=55)
    ax.set(ylabel='count', title='All 21 fixtures: episodes versus unique particles')
    ax.legend(ncol=2, fontsize=9)
    save(fig, 'episode_counts', 'Execution succeeded for all cases. Extra episodes and '
         'censored observations remain physical limitations, not scientific passes.')
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for ax, field in zip(axes, ('pressure', 'temperature')):
        arrays = [[e[field + '_mach'] - measurements[n]['expected_input_mach']
                   for e in measurements[n]['event_evidence']
                   if math.isfinite(e[field + '_mach'])] for n in names]
        ax.boxplot(arrays, tick_labels=names, showfliers=True)
        ax.axhline(0, color='black', ls='--')
        ax.set_ylabel(field.capitalize() + ' Mach − input Mach')
        ax.grid(axis='y', alpha=.2)
    axes[0].set_title('First-inside event strength: bias and spread against independent input Mach')
    axes[1].tick_params(axis='x', rotation=55)
    save(fig, 'mach_errors', 'Every emitted event contributes, including weak secondary lobes. '
         'No outlier removal, calibrated endpoint selection or relaxed tolerance is applied.')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for label in ('nx64', 'seed11', 'nx256', 'cfl0.15', 'cfl0.45'):
        for ax, key in zip(axes, ('entry_region_bracket', 'exit_region_bracket')):
            widths = sorted(e[key][1] - e[key][0] for e in measurements[label]['event_evidence']
                            if all(math.isfinite(t) for t in e[key]))
            ax.step(widths, np.arange(1, len(widths) + 1) / len(widths), label=label)
            ax.set(xlabel=key.split('_')[0] + ' bracket width (code time)', ylabel='cumulative fraction')
            ax.grid(alpha=.2)
    axes[0].legend(fontsize=9)
    fig.suptitle('Resolution and cadence change recorded timing brackets')
    save(fig, 'timing_bracket_widths', 'Brackets bound adjacent region-membership observations. '
         'Their widths are not errors against an invented exact MC crossing time.')
    fig, axes = plt.subplots(3, 1, figsize=(12, 11))
    for ax, label in zip(axes, ('seed11', 'mach5', 'boost2')):
        rows, events, suspended, ids = histories[label]
        selected_ids = sorted(set(ids) | set(range(16)))
        index = {pid: i for i, pid in enumerate(selected_ids)}
        for e in events:
            if e['particle_id'] not in index:
                continue
            y = index[e['particle_id']]
            color = '#cc3333' if e['sequence'] else '#337eb8'
            ax.plot([e['first']['time'], e['last']['time']], [y, y], color=color, lw=3)
            for left, right in (('pre', 'first'), ('last', 'post')):
                a, b = e[left]['time'], e[right]['time']
                if math.isfinite(a) and math.isfinite(b):
                    ax.plot([a, b], [y, y], color=color, lw=8, alpha=.3)
            ax.scatter(e['first']['time'], y, color=color, s=12)
        for e in suspended:
            if e['particle_id'] in index:
                ax.scatter(e['last']['time'], index[e['particle_id']], marker='x', c='black', s=35)
        ax.set_yticks(list(index.values()), list(index.keys()), fontsize=8)
        ax.set(xlabel='code time', ylabel='particle ID', title=label + ': blue first episode; red rearmed episode; × suspended')
        ax.grid(axis='x', alpha=.2)
    save(fig, 'encounter_brackets', 'Thick pale ends show recorded entry/exit brackets; '
         'narrow solid segments join first/last inside observations, not exact physical residence.')

    observations, times = thermo(work / 'integration/uniform_on')
    oracle = moment_oracle(observations, times)
    assert oracle['mean_pass'] and oracle['variance_pass']
    generations = sorted(times)
    ids = sorted(p for g, p in observations if g == generations[0])
    hops = []
    for pid in ids:
        xx = [observations[g, pid]['x'][0] for g in generations]
        hops.append(round(sum((b - a) % 1 for a, b in zip(xx, xx[1:])) * 64))
    probabilities = [(times[b] - times[a]) * 64 for a, b in zip(generations, generations[1:])]
    distribution = np.array([1.])
    for p in probabilities:
        distribution = np.convolve(distribution, [1 - p, p])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].hist(hops, bins=np.arange(min(hops) - .5, max(hops) + 1.5),
                 density=True, alpha=.65, label='1024 MC tracers')
    axes[0].plot(np.arange(len(distribution)), distribution, 'k.-', label='independent Bernoulli sum')
    axes[0].set(xlim=(min(hops) - 2, max(hops) + 2), xlabel='net positive cell hops', ylabel='probability')
    axes[0].legend(fontsize=9)
    axes[1].bar(['mean', 'variance'], [oracle['mean_standard_errors'], oracle['variance_standard_errors']])
    axes[1].axhline(6, color='gray', ls='--')
    axes[1].axhline(-6, color='gray', ls='--')
    axes[1].set(ylabel='deviation / predeclared standard error', title='Existing ±6 SE diagnostic')
    save(fig, 'uniform_transport', 'Physical flux probabilities define this independent '
         'ensemble reference. A fixed seed diagnostic does not prove RNG independence at research scales.')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for label in ('uniform_on', 'contact'):
        rows = traces(work / 'integration' / label / 'shock_events')
        first, last_g = min(g for g, p in rows), max(g for g, p in rows)
        for g, ls in ((first, '--'), (last_g, '-')):
            sample = sorted((o for (generation, p), o in rows.items() if generation == g),
                            key=lambda o: o['x'][0])
            axes[0].plot([o['x'][0] for o in sample], [o['rho'] for o in sample],
                         ls, label=f'{label}, generation {g}')
        assert not read_history(work / 'integration' / label / 'shock_events')
    axes[0].set(xlabel='sampled parent x', ylabel='density', title='Uniform flow and advected contact')
    axes[0].legend(fontsize=8)
    axes[1].bar(['uniform', 'contact', 'zero particles'], [0, 0, 0])
    axes[1].set(ylim=(0, 1), ylabel='emitted shock episodes', title='Negative controls: zero false detections')
    for i in range(3):
        axes[1].text(i, .04, '0', ha='center')
    save(fig, 'negative_controls', 'Contact samples show a density change at constant '
         'pressure and velocity. Uniform/contact/empty fixtures emitted no events.')
    reference, _ = thermo(work / 'integration/planar')
    events = read_history(work / 'integration/planar/shock_events')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    restart_results = []
    for ax, phase in zip(axes, ('outside', 'active', 'pending')):
        restored, _ = thermo(work / 'integration' / ('resume_' + phase))
        assert all(reference[key] == value for key, value in restored.items())
        replay = read_history(work / 'integration' / ('seed_' + phase) / 'shock_events')
        assert canonical(events) == canonical(replay)
        delta_x = [abs(value['x'][0] - reference[key]['x'][0]) for key, value in restored.items()]
        delta_rho = [abs(value['fields']['density'] - reference[key]['fields']['density'])
                     for key, value in restored.items()]
        ax.bar(['max |Δx|', 'max |Δρ|'], [max(delta_x) / 1e-14, max(delta_rho) / 1e-14])
        ax.set(ylim=(0, 1), title=phase + ' checkpoint')
        ax.text(.5, .6, f'{len(restored):,} exact rows\n{len(replay)} exact events\n0 duplicate identities',
                ha='center', fontsize=10, transform=ax.transAxes)
        restart_results.append({'phase': phase, 'rows_exact': len(restored), 'events_exact': len(replay)})
    axes[0].set_ylabel('absolute difference / 1e−14 code units')
    save(fig, 'restart_agreement', 'Replay from outside, active and pending-exit states '
         'preserved exact same-partition particle rows and event identities. No cross-rank claim.')

    # Run timings sequentially, after all suites; rotate mode order across repeats.
    benchmark = []
    modes = ('off', 'tracking', 'trace')
    for nx in (64, 256):
        for trial in range(3):
            for mode in modes[trial:] + modes[:trial]:
                changes = uniform_changes(count=1024)
                changes.update({'mesh/nx1': nx, 'meshblock/nx1': nx, 'time/tlim': .1,
                    'shock_tracking/enabled': str(mode != 'off').lower(),
                    'shock_tracking/trace': str(mode == 'trace').lower(),
                    'output2/dcycle': 100000, 'output3/dt': -1})
                label = f'benchmark_{nx}_{mode}_{trial}'
                item = execute(label, changes)
                log = (work / label / 'run.log').read_text()
                item.update(mesh=[nx, 4, 1], particles=1024, mode=mode, trial=trial)
                item['solver_cpu_seconds'] = float(re.search(r'cpu time used\s*=\s*(\S+)', log)[1])
                item['zone_cycles_per_cpu_second'] = float(re.search(r'zone-cycles/cpu_second\s*=\s*(\S+)', log)[1])
                match = re.search(r'diagnostic_seconds_this_process=(\S+)', log)
                item['diagnostic_seconds'] = float(match[1]) if match else 0
                files = list((work / label / 'shock_events').glob('*'))
                item['shock_output_bytes'] = sum(p.stat().st_size for p in files if p.is_file())
                item['shock_output_files'] = sum(p.is_file() for p in files)
                benchmark.append(item)
                print(f'Finished {label}: {item["elapsed_seconds"]:.3f}s', flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    summary = []
    for j, mode in enumerate(modes):
        grouped = [[r for r in benchmark if r['mesh'][0] == nx and r['mode'] == mode]
                   for nx in (64, 256)]
        medians = [statistics.median(r['elapsed_seconds'] for r in group) for group in grouped]
        positions = np.arange(2) + (j - 1) * .23
        axes[0].bar(positions, medians, width=.23, label=mode)
        for pos, group in zip(positions, grouped):
            axes[0].scatter([pos] * 3, [r['elapsed_seconds'] for r in group], s=12, c='black')
        ratios = [median / statistics.median(r['elapsed_seconds'] for r in benchmark
                  if r['mesh'][0] == nx and r['mode'] == 'off') for nx, median in zip((64, 256), medians)]
        axes[1].bar(positions, ratios, width=.23, label=mode)
        for nx, median, ratio in zip((64, 256), medians, ratios):
            summary.append({'nx': nx, 'mode': mode, 'median_whole_process_seconds': median,
                            'ratio_to_tracking_off': ratio})
    for ax in axes:
        ax.set_xticks([0, 1], ['64 × 4 cells', '256 × 4 cells'])
        ax.legend(fontsize=9)
    axes[0].set(ylabel='whole-process wall seconds', title='1024 tracers: median of 3 runs; dots = repeats')
    axes[1].axhline(1, color='gray', ls='--')
    axes[1].set(ylabel='runtime / tracking-off runtime', title='Same fluid/tracers; only diagnostics differ')
    save(fig, 'runtime_cost', 'Release Serial CPU, uniform flow, tlim=0.1. Startup, '
         'durable empty event chunks and output are included. These are laptop costs, not cluster budgets.')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for j, mode in enumerate(modes):
        groups = [[r for r in benchmark if r['mesh'][0] == nx and r['mode'] == mode] for nx in (64, 256)]
        positions = np.arange(2) + (j - 1) * .23
        for ax, key, scale in zip(axes, ('shock_output_bytes', 'shock_output_files'), (1e6, 1)):
            ax.bar(positions, [statistics.median(r[key] for r in group) / scale for group in groups],
                   width=.23, label=mode)
            ax.set_xticks([0, 1], ['64 × 4 cells', '256 × 4 cells'])
            ax.legend(fontsize=9)
    axes[0].set(ylabel='shock diagnostics / MB (decimal)', title='Bytes actually written')
    axes[1].set(ylabel='shock diagnostic file count', title='Immutable chunks and optional traces')
    save(fig, 'output_cost', 'Actual directory sizes/counts exclude fluid and thermo output. '
         'Trace mode records every observation/transfer and is intended for small qualification jobs.')
    pdf.close()
    result = dict(identity, platform=platform.platform(), python=sys.version,
                  numpy=np.__version__, matplotlib=matplotlib.__version__, commands=commands,
                  raw_workdir=str(work), snapshot_run=snapshot_run, independent_transport=oracle,
                  restart_results=restart_results, benchmark_runs=benchmark, benchmark_summary=summary,
                  figures=figures, scientific_signoff='pending',
                  scope='Serial CPU double; planar thin 2D. One thin 3D heldout. No MPI/GPU/cluster timing.',
                  known_limits={n: measurements[n] for n in ('mach5', 'boost2', 'mach1.5')})
    (output / 'report.json').write_text(json.dumps(result, indent=2, allow_nan=True) + '\n')
    cards = ''.join('<article><h2>' + html.escape(f['name'].replace('_', ' ')) + '</h2>'
                    '<a href="' + f['name'] + '.png"><img loading="lazy" src="' + f['name'] + '.png"></a>'
                    '<p>' + html.escape(f['caption']) + '</p></article>' for f in figures)
    document = ('<!doctype html><html lang="en"><meta charset="utf-8">'
                '<title>AthenaK local shock checks</title><style>'
                'body{font:17px system-ui;max-width:1200px;margin:32px auto;padding:0 20px;background:#fafafa;color:#182434}'
                'img{width:100%}article{background:white;padding:20px;margin:28px 0;border:1px solid #ddd;border-radius:10px}'
                'a{color:#1967a4}p{line-height:1.5}code{overflow-wrap:anywhere}</style>'
                '<h1>AthenaK: laptop shock checks</h1><p>' + str(len(figures)) + ' plots from fresh local tests. '
                'Click a plot for full resolution. <a href="all-plots.pdf">Download all plots as PDF</a>. '
                '<a href="report.json">Measurements and source identity</a>.</p>'
                '<p>All 20 engineering processes and 21 physical fixture executions succeeded. '
                '<strong>Scientific acceptance remains pending.</strong> Mach 5 has 134 episodes / 128 unique particles; '
                'boosted Mach 2 has 132 / 127 plus one suspended; Mach 1.5 has 127 plus one suspended. '
                'These are preserved failures/limits. The tracer makes stochastic cell hops; plots do not '
                'assign exact material crossing times.</p><p>Source fingerprint: <code>' +
                identity['compiled_source_fingerprint'] + '</code>. Serial CPU, double, RK2, uniform Cartesian ideal gas.</p>'
                + cards + '</html>')
    (output / 'index.html').write_text(document)
    files = {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    (output / 'sha256.json').write_text(json.dumps(files, indent=2) + '\n')
    assert len(figures) == 28 and all((output / (f['name'] + '.png')).stat().st_size > 1000 for f in figures)
    print(json.dumps({'plots': len(figures), 'report': str(output / 'index.html'),
                      'runtime_medians': summary}, indent=2))


if __name__ == '__main__':
    main()
