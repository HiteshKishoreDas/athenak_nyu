"""Bounded decaying supersonic box and measured MC encounter statistics.

Use existing NumPy/Matplotlib, native binary/thermo readers and strict event reader.
This low-resolution pilot is not a calibrated physical shock census.
"""
import argparse
import csv
from collections import Counter
import html
import json
from pathlib import Path
import re
import subprocess
import sys
import time
import zipfile

from common import ROOT, edit_input, provenance, sha, thermo
sys.path.insert(0, str(ROOT / 'vis/python'))
from bin_convert import read_binary
from shock_history import read_history, read_suspended
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FLAGS = {'left_censored': 1, 'invalid_observation': 2, 'merged_flicker': 4,
         'no_plateau': 8, 'cadence_limited': 16, 'temperature_invalid': 32,
         'right_censored': 64, 'domain_exit': 128, 'direction_ambiguous': 256,
         'downstream_reentry': 512}
N = 4096
GAMMA = 1.4
TLIM = .2


def snapshots(directory):
    result = []
    for path in sorted((directory / 'bin').glob('*.bin')):
        d = read_binary(path)
        assert d['n_mbs'] == 1
        f = {k: v[0] for k, v in d['mb_data'].items()}
        assert f['dens'].shape == (32, 32, 32) and f['dens'].dtype == np.float64
        r = f['dens']
        ke = sum(f[k] ** 2 for k in ('mom1', 'mom2', 'mom3')) / (2 * r)
        pressure = (GAMMA - 1) * (f['ener'] - ke)
        assert np.isfinite(list(f.values())).all() and np.all(r > 0)
        assert np.all(pressure > 0)
        vel = [f[k] / r for k in ('mom1', 'mom2', 'mom3')]
        div = sum((np.roll(v, -1, axis=a) - np.roll(v, 1, axis=a)) * 16
                  for v, a in zip(vel, (2, 1, 0)))
        result.append({'path': path, 'time': d['time'], 'cycle': d['cycle'], 'fields': f,
                       'rho': r, 'pressure': pressure, 'div': div,
                       'mass': float(r.mean()), 'energy': float(f['ener'].mean()),
                       'kinetic': float(ke.mean()), 'thermal': float((f['ener']-ke).mean()),
                       'mach': float(np.sqrt(np.mean(2 * ke / (GAMMA * pressure))))})
    assert result[0]['time'] == 0 and result[-1]['time'] == TLIM
    return result


def analyze(directory, field):
    events = read_history(directory / 'shock_events')
    suspended = read_suspended(directory / 'shock_events')
    all_events = events + suspended
    ids = [e['identity'] for e in all_events]
    assert len(set(ids)) == len(ids)
    observations, times = thermo(directory)
    populations = {g: {p for generation, p in observations if generation == g} for g in times}
    assert all(p == set(range(N)) for p in populations.values())
    for e in all_events:
        assert e['particle_id'] in populations[min(times)]
        assert json.dumps(e['representative'], sort_keys=True) == json.dumps(e['first'], sort_keys=True)
    for pid in range(N):
        seq = sorted(e['sequence'] for e in all_events if e['particle_id'] == pid)
        assert seq == list(range(len(seq))), (pid, seq)
    counts = Counter(e['particle_id'] for e in events)
    begins = Counter(e['particle_id'] for e in all_events)
    flags = {name: sum(bool(e['flags'] & bit) for e in all_events) for name, bit in FLAGS.items()}
    pressure_mach = [e['first']['mach_pressure'] for e in events if np.isfinite(e['first']['mach_pressure'])]
    temp_mach = [e['first']['mach_temperature'] for e in events if np.isfinite(e['first']['mach_temperature'])]
    counters = re.findall(r'^Shock rank=0 (.+)$', (directory/'run.log').read_text(), re.M)
    assert len(counters) == 1
    counters = {k: float(v) if 'seconds' in k else int(v) for k, v in
                re.findall(r'(\w+)=([^ ]+)', counters[0])}
    assert counters['complete'] == sum(e['status'] == 0 for e in events)
    assert counters['partial'] == sum(e['status'] == 1 for e in events)
    assert counters['suspended'] == len(suspended)
    assert counters['sampled'] == (field[-1]['cycle'] + 1) * N
    mass_error = max(abs(d['mass'] - field[0]['mass']) for d in field)
    energy_error = max(abs(d['energy'] / field[0]['energy'] - 1) for d in field)
    assert mass_error <= 1e-12 and energy_error <= 1e-12
    summary = {'completed': counters['complete'], 'partial': counters['partial'],
               'suspended': len(suspended), 'particles': N,
               'particles_with_emitted_events': len(counts),
               'particles_with_any_recorded_episode': len(begins),
               'particles_with_no_recorded_episode': N-len(begins),
               'mean_emitted_episodes_per_particle': len(events)/N,
               'mean_recorded_starts_per_particle': len(all_events)/N,
               'maximum_emitted_per_particle': max(counts.values(), default=0),
               'emitted_count_histogram': dict(sorted(Counter(counts[p] for p in range(N)).items())),
               'started_count_histogram': dict(sorted(Counter(begins[p] for p in range(N)).items())),
               'flags_all_recorded_episodes': flags,
               'pressure_mach_quantiles': np.quantile(pressure_mach, [.1,.5,.9]).tolist() if pressure_mach else None,
               'temperature_mach_quantiles': np.quantile(temp_mach, [.1,.5,.9]).tolist() if temp_mach else None,
               'finite_pressure_estimates': len(pressure_mach), 'finite_temperature_estimates': len(temp_mach),
               'counters': counters, 'particle_ids_preserved': True, 'unique_contiguous_event_ids': True,
               'max_mass_absolute_error': mass_error, 'max_total_energy_relative_error': energy_error,
               'field_evolution': [{k:d[k] for k in ('time','cycle','mass','energy','kinetic','thermal','mach')} for d in field]}
    return summary, events, suspended


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reuse', action='store_true')
    args = p.parse_args()
    work, output, binary = args.workdir.resolve(), args.output.resolve(), args.executable.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {'test_id': 'PHY16', 'scientific_signoff': 'pending',
              'scope': 'Exploratory 32^3 periodic decaying ideal-gas box, initial rms Mach3; no forcing/cooling.',
              'definitions': {'mach': 'sqrt(volume mean(v^2/(gamma P/rho)))',
                              'strength': 'first valid inside spatial detector estimate; not plateau-certified',
                              'counts': 'numerical episodes, not certified unique physical shocks',
                              'timing': 'recorded prior-outside / first-inside brackets, not exact crossings'},
              'experimental_settings': {'width': 4, 'min_mach': 1.1, 'velocity_seed': 24680,
                                       'particle_seed': 11, 'seed_schedule_seed': 903, 'tlim': TLIM},
              **provenance(binary)}
    if not args.reuse:
        work.mkdir(parents=True, exist_ok=False)
        report['runs'] = []
        base = (ROOT/'inputs/tests/shock_mc_turb_decay.athinput').read_text()
        for name, changes in [('baseline', {}), ('off', {'shock_tracking/enabled':'false'}),
                              ('half_cfl', {'time/cfl_number':.15, 'shock_tracking/lineage':'exploratory-decay-m3-32-cfl015'})]:
            directory = work/name
            directory.mkdir()
            inp = directory/'input.athinput'
            inp.write_text(edit_input(base, changes))
            command = [str(binary), '-i', str(inp)]
            start = time.monotonic()
            with (directory/'run.log').open('w') as log:
                try:
                    status = subprocess.run(command, cwd=directory, stdout=log,
                                            stderr=subprocess.STDOUT, timeout=180).returncode
                except subprocess.TimeoutExpired:
                    status = 'timeout'
            report['runs'].append({'name': name, 'command': command, 'exit':status,
                                   'wall_seconds':time.monotonic()-start,
                                   'input_sha256':sha(inp), 'log_sha256':sha(directory/'run.log')})
            (work/'execution.json').write_text(json.dumps(report, indent=2)+'\n')
            assert status == 0, report['runs'][-1]
    else:
        old = json.loads((work/'execution.json').read_text())
        assert old['executable_sha256'] == report['executable_sha256']
        report['runs'] = old['runs']
    fields = {name:snapshots(work/name) for name in ('baseline','off','half_cfl')}
    initial = fields['baseline'][0]
    assert abs(initial['mach']-3) <= 1e-12
    assert np.max(abs(initial['pressure']-1/GAMMA)) <= 1e-12
    assert np.array_equal(initial['rho'], np.ones((32,32,32)))
    assert len(fields['baseline']) == len(fields['off'])
    for a,b in zip(fields['baseline'], fields['off']):
        assert (a['time'],a['cycle']) == (b['time'],b['cycle'])
        assert all(np.array_equal(a['fields'][k],b['fields'][k]) for k in a['fields'])
    assert thermo(work/'baseline') == thermo(work/'off')
    assert not (work/'off/shock_events').exists()
    report['non_interference'] = 'Exact equality of all saved double conserved volumes and sampled native particle records.'
    report['cases'], details = {}, {}
    for name in ('baseline','half_cfl'):
        summary, events, suspended = analyze(work/name, fields[name])
        report['cases'][name] = summary
        details[name] = (events,suspended)
    report['cfl_count_change_fraction'] = report['cases']['half_cfl']['completed']/report['cases']['baseline']['completed']-1
    report['limitations'] = ['One realization, short low-resolution decaying run; not steady driven turbulence.',
                             'CFL change affects the fluid solution, MC random draws and detector sampling jointly.',
                             'Episodes, first-inside strengths and reentries require scientific qualification.',
                             'All episodes retain no_plateau/cadence flags; open episodes remain right-censored.',
                             'No independent truth shock catalogue, recall, physical waiting-time law or HPC qualification.']
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False,'axes.spines.right':False})
    figures=[]
    def save(fig,name,description):
        fig.tight_layout()
        fig.savefig(output/(name+'.png'),dpi=150)
        plt.close(fig)
        figures.append({'name':name,'description':description})
    chosen=[fields['baseline'][i] for i in (0,2,4,len(fields['baseline'])-1)]
    fig, axes=plt.subplots(3,4,figsize=(14,10))
    for col,d in enumerate(chosen):
        for row,(key,label) in enumerate([('rho','log10 density'),('pressure','pressure'),('div','compression -div(v)')]):
            value=np.log10(d[key]) if row==0 else d[key] if row==1 else np.maximum(-d[key],0)
            planes = [np.log10(s[key][16]) if row == 0 else s[key][16] if row == 1 else np.maximum(-s[key][16], 0) for s in chosen]
            lo, hi = min(v.min() for v in planes), max(v.max() for v in planes)
            im=axes[row,col].imshow(value[16],origin='lower',extent=(0,1,0,1),vmin=lo,vmax=hi,cmap='magma' if row<2 else 'Blues')
            axes[row,col].set_title(f'{label}, t={d["time"]:.4f}')
            fig.colorbar(im,ax=axes[row,col],shrink=.75)
    save(fig,'box_slices','Real z=0.515625 slices of the 32^3 box. Each row shares a color scale across times; compression is not the detector mask.')
    fig, ax=plt.subplots(1,3,figsize=(14,4))
    for name,label in [('baseline','CFL 0.30'),('half_cfl','CFL 0.15')]:
        f=fields[name];t=[d['time'] for d in f]
        ax[0].plot(t,[d['mach'] for d in f],'-o',label=label)
        ax[1].plot(t,[d['kinetic'] for d in f],label=label+' kinetic')
        ax[1].plot(t,[d['thermal'] for d in f],'--',label=label+' thermal')
        ax[2].plot(t,[d['energy']/f[0]['energy']-1 for d in f],label=label)
    ax[0].axhline(1,color='k',ls=':');ax[0].set_ylabel('Volume rms local Mach')
    ax[1].set_ylabel('Energy per box volume');ax[2].set_ylabel('Relative total-energy change')
    for a in ax:a.set_xlabel('Code time');a.legend(fontsize=8)
    save(fig,'fluid_evolution','The box begins supersonic and heats as kinetic energy decays. This is not stationary turbulence. Conservation is checked independently at 1e-12.')
    fig, ax=plt.subplots(1,3,figsize=(14,4))
    bins=np.arange(0,TLIM+.001,.025)
    for name,label in [('baseline','CFL 0.30'),('half_cfl','CFL 0.15')]:
        e,s=details[name];c=Counter(x['particle_id'] for x in e)
        hist=Counter(c[i] for i in range(N));ax[0].plot(list(hist),[hist[k] for k in hist],'o',label=label)
        starts=np.array([x['first']['time'] for x in e+s]);ordered=np.sort(starts)
        ax[1].step(ordered,np.arange(1,len(ordered)+1)/N,where='post',label=label)
        rates,_=np.histogram(starts,bins);ax[2].stairs(rates/(N*np.diff(bins)),bins,label=label)
    ax[0].set_xlabel('Emitted episodes per particle');ax[0].set_ylabel('Particles (zero encounters included)')
    ax[1].set_xlabel('First inside observation time');ax[1].set_ylabel('Cumulative recorded starts / particle')
    ax[2].set_xlabel('First inside observation time');ax[2].set_ylabel('Recorded starts / particle / code time')
    for a in ax:a.legend()
    save(fig,'encounter_counts','Numerical episode counts, including all 4096 particles. Starts include terminal suspended episodes; emitted counts exclude those open episodes. CFL changes are a joint numerical sensitivity check.')
    fig, ax=plt.subplots(1,2,figsize=(11,4))
    for name,label in [('baseline','CFL 0.30'),('half_cfl','CFL 0.15')]:
        for a,key in zip(ax,('mach_pressure','mach_temperature')):
            values=sorted(e['first'][key] for e in details[name][0] if np.isfinite(e['first'][key]))
            a.step(values,np.arange(1,len(values)+1)/len(values),label=label)
            a.set_xlabel('First-inside '+key.replace('_',' '));a.set_ylabel('Emitted-event cumulative fraction');a.legend()
    save(fig,'strength_estimates','These are spatial detector estimates at the first inside sample, not Mach numbers from certified plateaus. Both estimators are shown; particle encounter counts are not a truth shock catalogue.')
    fig, ax=plt.subplots(1,2,figsize=(13,5))
    names=['merged_flicker','direction_ambiguous','downstream_reentry','invalid_observation','right_censored','temperature_invalid']
    for i,(name,label) in enumerate([('baseline','CFL 0.30'),('half_cfl','CFL 0.15')]):
        summary=report['cases'][name]
        ax[0].barh(np.arange(len(names))+i*.35,[summary['flags_all_recorded_episodes'][n] for n in names],height=.35,label=label)
        events=details[name][0]
        widths=sorted(e['first']['time']-e['pre']['time'] for e in events if e['pre']['valid'])
        ax[1].step(widths,np.arange(1,len(widths)+1)/len(widths),label=label)
    ax[0].set_yticks(np.arange(len(names))+.175,names);ax[0].set_xlabel('Recorded episodes with flag (overlap allowed)');ax[0].legend()
    ax[1].set_xlabel('Observed entry bracket width (code time)');ax[1].set_ylabel('Cumulative fraction');ax[1].legend()
    save(fig,'quality_and_timing','Flags overlap. All episodes additionally have no_plateau and cadence_limited. Entry brackets span recorded prior-outside to first-inside observations; gaps can include suppressed reentry. They are not exact crossing errors. Suspended episodes retain right censoring.')
    report['source_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip()
    report['source_sha256'] = {name:sha(ROOT/name) for name in [
        'src/shocks/detector.hpp', 'src/shocks/encounter.hpp',
        'src/srcterms/initial_perturbations.cpp', 'tst/shocks-mc/turbulence_pilot.py',
        'inputs/tests/shock_mc_turb_decay.athinput']}
    with (output/'episodes.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['case', 'outcome', 'status', 'particle_id', 'sequence', 'flags',
                         'entry_lower', 'first_inside', 'last_inside', 'post_observation',
                         'pressure_mach', 'temperature_mach'])
        for name, (events, suspended) in details.items():
            for kind, rows in [('emitted', events), ('suspended', suspended)]:
                for e in rows:
                    writer.writerow([name, kind, e['status'], e['particle_id'], e['sequence'], e['flags'],
                                     e['pre']['time'] if e['pre']['valid'] else '',
                                     e['first']['time'], e['last']['time'],
                                     e['post']['time'] if e['post']['valid'] else '',
                                     e['first']['mach_pressure'], e['first']['mach_temperature']])
    report['archive_scope'] = 'All canonical event/suspended chunks and inputs/logs; four baseline volumes and half-CFL endpoint volumes. Full local raw data are SHA256 inventoried separately.'
    report['figures']=figures
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    archive=output/'selected-raw-evidence.zip'
    volumes = {d['path'] for d in chosen}
    volumes.update(fields['half_cfl'][i]['path'] for i in (0,-1))
    raw_hashes = {}
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as z:
        for name in ('baseline','off','half_cfl'):
            for path in sorted((work/name).rglob('*')):
                if not path.is_file():
                    continue
                relative = str(path.relative_to(work))
                raw_hashes[relative] = sha(path)
                if path.parent.name == 'shock_events' or path in volumes or path.parent == work/name:
                    z.write(path,relative)
        z.write(work/'execution.json','execution.json')
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    (output/'raw-sha256.json').write_text(json.dumps(raw_hashes,indent=2)+'\n')
    text=f'<h1>Decaying Mach 3 box: exploratory particle shock statistics</h1><p>32³ cells, 4096 mass-seeded MC tracers, ideal gas gamma=1.4, periodic box, t=0.2. No forcing or cooling. Scientific sign-off pending.</p>'
    for name,s in report['cases'].items():
        text+=f'<p><b>{name}</b>: {s["completed"]} complete, {s["partial"]} partial, {s["suspended"]} suspended; {s["particles_with_emitted_events"]}/4096 particles have emitted episodes; mean {s["mean_emitted_episodes_per_particle"]:.3f} emitted episodes per particle.</p>'
    text+='<p>'+html.escape(report['non_interference'])+'</p><p>'+' '.join(html.escape(v) for v in report['limitations'])+'</p>'
    for f in figures:
        text+=f'<article><h2>{f["name"].replace("_"," ")}</h2><p>{html.escape(f["description"])}</p><a href="{f["name"]}.png"><img src="{f["name"]}.png"></a></article>'
    (output/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Supersonic box shock pilot</title><style>body{font:17px system-ui;max-width:1200px;margin:32px auto;padding:20px}p{line-height:1.5}img{width:100%}article{margin:40px 0}</style>'+text)
    (output/'sha256.json').write_text(json.dumps({str(p.relative_to(output)):sha(p) for p in sorted(output.iterdir()) if p.is_file()},indent=2)+'\n')
    print(json.dumps({'output':str(output),'counts':{n:{k:s[k] for k in ('completed','partial','suspended','particles_with_emitted_events')} for n,s in report['cases'].items()},'cfl_count_change_fraction':report['cfl_count_change_fraction']}))


if __name__=='__main__':
    main()
