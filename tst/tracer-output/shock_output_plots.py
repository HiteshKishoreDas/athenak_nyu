"""Plot actual native THP/BIN output checks; no reconstructed trajectories."""
import argparse
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from shock_output_check import FIELDS,mesh_frames,mesh_sample,equal
from check import normalized
from tracer_history import read_run
from tracer_thermo_history import read_curves
from shock_history import read_trace as particle_trace
from curve_shock_history import read_trace as curve_trace


def paint(ax,frame,field,**kwargs):
    image=None
    if 'vmin' not in kwargs:
        all_values=np.asarray(frame['mb_data'][field]);finite=all_values[np.isfinite(all_values)]
        if finite.size:
            low,high=float(np.min(finite)),float(np.max(finite))
            kwargs.update(vmin=low if low!=high else low-.1,vmax=high if low!=high else high+.1)
        else:kwargs.update(vmin=0,vmax=1)
    for bounds,values in zip(frame['mb_geometry'],frame['mb_data'][field]):
        low,high,y0,y1,_,_=bounds
        image=ax.pcolormesh(np.linspace(low,high,values.shape[2]+1),
            np.linspace(y0,y1,values.shape[1]+1),values[0],shading='flat',**kwargs)
    ax.set(xlabel=r'$x_1$ [code units]',ylabel=r'$x_2$ [code units]',
           xlim=(frame['x1min'],frame['x1max']),ylim=(frame['x2min'],frame['x2max']))
    return image


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True,help='passing shock_output_check.py workdir')
    parser.add_argument('--thermo-runs',type=Path,required=True,help='passing standard check.py workdir')
    parser.add_argument('--csv-reference',type=Path,help='optional retained CSV controls for matched-size plot')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.run/'report.json').read_text());assert all(c['status']=='pass' for c in report['checks'])
    figures=[]
    def save(fig,name,caption):
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
        figures.append(dict(file=name+'.png',caption=caption))
    directory=args.run/'planar';data=read_run(directory);frames=mesh_frames(directory,'shock_mesh')
    cycle=sorted(frames)[len(frames)//2];frame=frames[cycle]
    fig,axes=plt.subplots(2,1,figsize=(13,7))
    for ax,field,title in zip(axes,['shock_region','shock_mach_pressure'],
            ['Unreduced native shock-region membership','Native pressure-jump Mach estimate (dimensionless)']):
        image=paint(ax,frame,field,cmap='Oranges' if field=='shock_region' else 'viridis',
            **({'vmin':0,'vmax':1} if field=='shock_region' else {}))
        plt.colorbar(image,ax=ax,label='membership [0/1]' if field=='shock_region' else 'Mach estimate')
        particles=[o for (g,_),o in data['particles'].items() if g==cycle]
        ax.scatter([o['x'][0] for o in particles],[o['x'][1] for o in particles],s=8,c='black',alpha=.5,label='actual particle samples')
        geometry=next(s for s in data['geometry'] if s['cycle']==cycle)['curves'][0]
        ax.plot([n['x'][0] for n in geometry['nodes']],[n['x'][1] for n in geometry['nodes']],'.-',c='dodgerblue',label='stored material-curve connectivity')
        ax.set_title(title);ax.legend(loc='upper right',fontsize=8)
    fig.suptitle(f'Existing BIN + THP readers · cycle {cycle}, time {frame["time"]:.6g}')
    save(fig,'native-shock-map','Actual saved planar-shock mesh and same-cycle particle/curve positions. The connected blue line uses stored curve geometry/connectivity; no particle trajectory is invented. White gaps in the Mach panel are unavailable NaN estimates. Region membership is the unreduced detector result. SHOUT01 checks all six properties against raw traces and native cells at every saved sample, beyond the single snapshot shown.')

    fig,axes=plt.subplots(1,2,figsize=(12,4.8));finite_total=0
    for key,reader,where,label,color in [('particles',particle_trace,directory/'shock_events','particles','tab:orange'),('curve_nodes',curve_trace,directory/'shock_events/curves','curve nodes','tab:blue')]:
        raw={}
        for path in where.glob('trace*.rank*.csv'):
            generation=int(path.name.split('.')[0][5:]);raw.update({(generation,p):o for p,o in reader(path).items()})
        rows=[(o,raw[k]) for k,o in data[key].items() if math.isfinite(o['fields']['shock_mach_pressure'])]
        x=[r['mach_pressure'] for _,r in rows];y=[o['fields']['shock_mach_pressure'] for o,_ in rows]
        axes[0].scatter(x,y,s=7,alpha=.4,label=label,c=color)
        delta=[mesh_sample(frames[o['cycle']],o,'shock_mach_pressure')-o['fields']['shock_mach_pressure'] for o,_ in rows]
        assert all(v==0 for v in delta)
        axes[1].scatter([o['time'] for o,_ in rows],delta,s=6,alpha=.4,label=label,c=color)
        finite_total+=len(rows)
    axes[0].set(xlabel='raw trace pressure-jump Mach',ylabel='THP pressure-jump Mach',title=f'Finite estimates: {finite_total} exact matches')
    axes[1].set(xlabel='saved sample time [code units]',ylabel='BIN minus THP Mach',title='Zero difference at actual containing cells')
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    save(fig,'native-property-parity','Each dot is an actual finite pressure-jump estimate. Left: raw trace versus THP. Right: native Real-precision BIN minus THP at the stored sample position/leaf. All plotted differences are exactly zero. SHOUT01 also checks membership, validity, temperature estimates and every unavailable NaN. Fluid Mach remains its separate existing quantity; no tolerance or physical crossing time is fitted.')

    before=mesh_frames(args.run/'restart_seed','shock_mesh')[13]
    after=mesh_frames(args.run/'restart_zero_step','shock_mesh')[13]
    fig,axes=plt.subplots(1,3,figsize=(15,4.6),gridspec_kw={'width_ratios':[1,1,1.2]})
    for ax,record,title in [(axes[0],before,'Saved cycle 13'),(axes[1],after,'Zero-step restart, same cycle/time')]:
        paint(ax,record,'shock_region',cmap='Oranges',vmin=0,vmax=1);ax.set_title(title)
    differences=[]
    for field in FIELDS:
        a,b=np.asarray(before['mb_data'][field]),np.asarray(after['mb_data'][field])
        differences.append(int(np.count_nonzero(~((a==b)|(np.isnan(a)&np.isnan(b))))))
    assert not any(differences)
    axes[2].barh(range(6),differences);axes[2].set_yticks(range(6),[f.replace('shock_','') for f in FIELDS],fontsize=8)
    axes[2].set(xlim=(0,1),xlabel='mismatched native cells',title='All six arrays: zero mismatches')
    for i,n in enumerate(differences):axes[2].text(.02,i,str(n),va='center')
    save(fig,'zero-step-cache-restart','Actual active-AMR checkpoint at cycle 13, before and after a restart taking zero fluid timesteps. The detector cache is rebuilt from synchronized restart primitives without advancing encounters/counters/journals. All six native BIN arrays and particle/curve THP properties match exactly; NaNs are compared as unavailable values. This specifically covers the restart case found during independent timing review.')

    rows=[]
    if args.csv_reference:
        for name in ['planar','amr_shock']:
            old=args.csv_reference/name/'curve_thermo_history';new=args.thermo_runs/name/'curve_thermo_history'
            a,_=read_curves(old);b,_=read_curves(new);assert normalized(a)==normalized(b)
            op=list(old.glob('*.csv'));np_=list(new.glob('*.thp'))
            rows.append(dict(case=name,samples=len(b),csv_bytes=sum(p.stat().st_size for p in op),thp_bytes=sum(p.stat().st_size for p in np_),csv_files=len(op),thp_files=len(np_)))
        fig,axes=plt.subplots(1,2,figsize=(11,4.6));x=np.arange(len(rows))
        axes[0].bar(x-.18,[r['csv_bytes']/1024 for r in rows],width=.36,label='retained CSV')
        axes[0].bar(x+.18,[r['thp_bytes']/1024 for r in rows],width=.36,label='native THP v2')
        axes[0].set(xticks=x,xticklabels=[r['case'] for r in rows],ylabel='stored curve thermo size [KiB]',title='Matched exact node/thermo records');axes[0].legend()
        axes[1].bar(x-.18,[r['csv_files'] for r in rows],width=.36,label='CSV');axes[1].bar(x+.18,[r['thp_files'] for r in rows],width=.36,label='THP v2')
        axes[1].set(xticks=x,xticklabels=[r['case'] for r in rows],ylabel='curve thermo file count',title='Existing cycle/rank chunks retained');axes[1].legend()
        (out/'storage.json').write_text(json.dumps(rows,indent=2)+'\n')
        save(fig,'native-thermo-storage','Measured file sizes and file counts for exact matched planar/live-AMR node thermo records, checked against retained signed-ID CSV controls. Only transport changes; identities, births, clocks, positions and sampled fields match exactly. The existing cycle/rank chunking remains. These bytes are local fixture measurements, not a filesystem throughput, cluster-memory or runtime benchmark.')

    forced=args.run/'forced'
    if forced.is_dir():
        data=read_run(forced);frames=mesh_frames(forced,'shock_mesh');ordered=sorted(frames)
        times=[frames[c]['time'] for c in ordered]
        volume=[float(np.mean(np.asarray(frames[c]['mb_data']['shock_region']))) for c in ordered]
        fig,axes=plt.subplots(1,2,figsize=(12,4.6));axes[0].plot(times,volume,'o-')
        axes[0].set(xlabel='saved time [code units]',ylabel='mesh volume fraction',title='Unreduced region in the native 32³ forced box')
        for key,label in [('particles','regular particles'),('curve_nodes','curve nodes')]:
            fractions=[]
            for c in ordered:
                samples=[o['fields']['shock_region'] for (g,_),o in data[key].items() if g==c]
                fractions.append(float(np.mean(samples)))
            axes[1].plot(times,fractions,'o-',label=label)
        axes[1].set(xlabel='saved time [code units]',ylabel='sampled occupancy fraction',title='Actual THP flags on the same saved clocks');axes[1].legend()
        for ax in axes:ax.grid(alpha=.2)
        save(fig,'native-forced-box-properties','Short accepted 32³ supersonic forced-box engineering fixture with actual native BIN and THP shock properties. Mesh fractions weight uniform cell volumes; particle fractions count mass-seeded regular tracers; curve fractions count nodes on the seeded material curve. Every saved property is checked against raw traces and current native cells by SHOUT08. Full encounters remain separately recorded. This short run establishes output coexistence/passivity, not converged turbulent statistics or scientific approval.')
    (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n')
    print(json.dumps(dict(figures=len(figures),storage=rows),indent=2))


if __name__=='__main__':main()
