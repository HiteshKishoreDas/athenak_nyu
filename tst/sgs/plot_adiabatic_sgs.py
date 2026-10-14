"""Render the saved compression results; no simulation is launched."""
import json
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True, help='completed diagnostic directory')
parser.add_argument('--output', type=Path, help='plot directory; defaults to input')
args = parser.parse_args()
root = args.output or args.input
root.mkdir(parents=True, exist_ok=True)
result = json.loads((args.input/'adiabatic-results.json').read_text())
profiles = result['profiles']
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
for ax, mode in zip(axes.flat, ['strong', 'compressible', 'near_isothermal', 'smooth']):
    for r in profiles:
        if r['mode'] == mode and r['cfl'] == .1:
            ax.plot(r['x'], r['entropy_ratio'], label=f"{r['grid']} cells")
    ax.axhline(1, color='black', ls='--', label='Constant entropy expectation')
    ax.set(xlim=(-.15, .15), xlabel='x / box length', ylabel=r'$K/(S_0\rho^{5/3})$',
           title=mode.replace('_', ' ').title())
    ax.legend(fontsize=8)
fig.suptitle('Native pressure-only SGS: shock entropy error versus smooth compression\n'
             'C1=C2=Cepsilon=Ckappa=0 for isolation; baseline model remains unchanged')
fig.savefig(root/'compression-profiles.png', dpi=160)
fig.savefig(root/'compression-profiles.pdf')
plt.close(fig)
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
for mode in ['weak', 'strong', 'compressible', 'near_isothermal', 'smooth']:
    rows = [r for r in profiles if r['mode'] == mode and r['cfl'] == .1]
    axes[0, 0].loglog([r['grid'] for r in rows], [r['mass_weighted_entropy_rms'] for r in rows],
                      'o-', label=mode.replace('_', ' '))
    axes[0, 1].plot([r['grid'] for r in rows], [r['total_K_over_adiabatic'] for r in rows],
                    'o-', label=mode.replace('_', ' '))
for mode in ['strong', 'compressible', 'near_isothermal']:
    rows = [r for r in profiles if r['mode'] == mode and r['grid'] == 256]
    axes[1, 0].plot([r['cfl'] for r in rows], [r['mass_weighted_entropy_rms'] for r in rows],
                    'o-', label=mode.replace('_', ' '))
d = result['decay']
axes[1, 1].loglog([r['cfl'] for r in d], [abs(r['relative_error']) for r in d], 'o-')
axes[0, 0].set(xlabel='Grid cells', ylabel='Mass-weighted entropy RMS error',
               title='Smooth error falls; shock error remains')
axes[0, 1].axhline(1, color='black', ls='--')
axes[0, 1].set(xlabel='Grid cells', ylabel='Mean K / adiabatic mean K',
               title='Whole-domain effect is an excess in these controls')
axes[1, 0].set(xlabel='CFL number', ylabel='Mass-weighted entropy RMS error',
               title='Spatial shock error persists with smaller timesteps')
axes[1, 1].set(xlabel='CFL number', ylabel='Relative native K error',
               title='Uniform decay versus ART analytic update')
for ax in axes.flat:
    ax.grid(alpha=.2)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)
fig.suptitle('Compression transport and decay integration are separate numerical questions')
fig.savefig(root/'compression-and-decay.png', dpi=160)
fig.savefig(root/'compression-and-decay.pdf')
plt.close(fig)
