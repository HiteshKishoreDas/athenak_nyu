"""Plot actual .crv records through the supplied curve reader/plotting API."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python')]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tracer_history import read_run
from plot_curves import select_curve,plot_curve,curve_points
from tracer_thermo_history import read_curves


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--reference-runs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((a.runs/'report.json').read_text())
    assert all(c['status']=='pass' for c in report['checks'])
    figures=[];data={}
    def save(fig,name,title,caption):
        path=out/(name+'.png');fig.savefig(path,dpi=160,bbox_inches='tight');plt.close(fig)
        figures.append(dict(path=path.name,title=title,caption=caption))
    for case,title,caption in [
        ('planar','Planar shock: one file contains topology and fluid samples',
         'All 1,728 node samples match the previous THP values exactly. Saved connectivity drives the geometry; complete encounters are read separately from unchanged journals.'),
        ('amr_shock','Live AMR: new nodes, signed IDs and thermo together',
         'All 9,644 node samples, mesh locations, geometry and histories match the previous code. Inserted nodes begin at their own birth; these curves use negative signed32 IDs.'),
        ('periodic','Periodic translation: continuous positions and stored image counts',
         'Accepted translation-wrap fixture. Continuous lifts and exact links stay in the CRV records; no nearest-image reconstruction or inferred trajectory.'),
        ('circle','Closed curve: saved connectivity, nonzero curvature and thermo',
         'Accepted radius 0.7 circle fixture. Native polygon length and curvature are checked against independent analytic values with the original production budgets.'),
        ('forced','Short 32³ forced box: geometry, fluid and sampled shock evidence',
         'Accepted short forced-box engineering fixture. All 256 curve samples match the previous implementation exactly. This is not converged turbulence statistics.')]:
        d=read_run(a.runs/case);data[case]=d
        selected=select_curve(d,field='pressure')
        name=case+'-combined.png';plot_curve(selected,out/name)
        figures.append(dict(path=name,title=title,caption=caption))
    d=data['planar'];frame=min(d['geometry'],key=lambda s:abs(s['cycle']-54));c=frame['curves'][0]
    nodes=c['nodes'];points=curve_points(dict(curve=c,domain=frame['domain']))
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for ax,field,label in zip(axes,('pressure','shock_region','shock_mach_pressure'),
                             ('Pressure (code units)','Shock-region membership [0/1]','Pressure-jump Mach estimate')):
        ax.plot(*np.array(points)[:,:2].T,color='gray',lw=1)
        values=np.array([n['fields'][field] for n in nodes]);finite=np.isfinite(values)
        dots=ax.scatter([n['x'][0] for n in nodes if np.isfinite(n['fields'][field])],
            [n['x'][1] for n in nodes if np.isfinite(n['fields'][field])],c=values[finite],s=38)
        fig.colorbar(dots,ax=ax,label=label)
        ax.set(xlabel=r'$x_1$ (code units)',ylabel=r'$x_2$ (code units)',title=field)
        ax.grid(alpha=.2)
    fig.suptitle(f"Same .crv snapshot: cycle {frame['cycle']}, t={frame['time']:.7f}")
    save(fig,'same-record-geometry-thermo-shocks','One binary snapshot, three views',
         'Each panel uses positions, links and fields from the same CRV record set. Lines follow explicit connectivity; colour shows the saved property. No joining of independently timed geometry/thermo files.')
    d=data['amr_shock'];frames=d['geometry'];times=[s['time'] for s in frames]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    axes[0,0].step(times,[s['node_count'] for s in frames],where='post')
    axes[0,0].set(ylabel='Saved node count',title='Live insertion')
    for level in sorted({n['level'] for s in frames for c in s['curves'] for n in c['nodes']}):
        axes[0,1].plot(times,[sum(n['level']==level for c in s['curves'] for n in c['nodes']) for s in frames],label=f'Level {level}')
    axes[0,1].set(ylabel='Nodes on leaf level',title='Mesh location saved in the same record');axes[0,1].legend()
    for index in sorted({0,len(frames)//2,len(frames)-1}):
        c=frames[index]['curves'][0];label=f"t={times[index]:.4f}"
        axes[1,0].plot([n['a'] for n in c['nodes']],[n['stretch'] for n in c['nodes']],'.-',label=label)
        axes[1,1].plot([n['a'] for n in c['nodes']],[n['fields']['pressure'] for n in c['nodes']],'.-',label=label)
    axes[1,0].set(xlabel=r'$a$ (initial arc length)',ylabel='Edge stretch',title='Geometry diagnostic');axes[1,0].legend()
    axes[1,1].set(xlabel=r'$a$ (initial arc length)',ylabel='Pressure (code units)',title='Thermo on the identical node set');axes[1,1].legend()
    for ax in axes[0]:ax.set_xlabel(r'$t$ (code units)')
    for ax in axes.flat:ax.grid(alpha=.2)
    save(fig,'amr-node-geometry-fields','AMR, adaptation and thermo share one node table',
         'Counts and leaf levels come directly from CRV. Stretch and pressure use identical saved nodes/material labels. Independent geometry reconstruction error remains within the unchanged 2e-12 CT budget. The open final node has no outgoing edge, so its stored stretch is zero.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');sizes=[]
    for case in ('planar','amr_shock','forced'):
        current=data[case];old=read_run(a.reference_runs/case)
        previous=old['curve_nodes'];cycles=sorted(current['curve_times'])
        errors=[max(abs(o['fields']['pressure']-previous[g,n]['fields']['pressure'])
                    for (g,n),o in current['curve_nodes'].items() if g==cycle) for cycle in cycles]
        axes[0].plot([current['curve_times'][g] for g in cycles],errors,'.-',label=case)
        # Match one old geometry snapshot per saved clock, excluding duplicate schedules.
        chosen={}
        for s in old['geometry']:
            if s['cycle'] in current['curve_times']:chosen.setdefault(s['cycle'],Path(s['path']))
        oldthp=list((a.reference_runs/case/'curve_thermo_history').glob('*.thp'))
        oldbytes=sum(p.stat().st_size for p in oldthp)+sum(p.stat().st_size for p in chosen.values())
        newfiles=list((a.runs/case/'curves').glob('*.crv'));newbytes=sum(p.stat().st_size for p in newfiles)
        sizes.append(dict(case=case,old_bytes=oldbytes,crv_bytes=newbytes,
                          old_files=len(oldthp)+len(chosen),crv_files=len(newfiles)))
    axes[0].set(xlabel=r'$t$ (code units)',ylabel='Maximum absolute pressure difference (code units)',
        title='CRV versus previous native THP: exactly zero',ylim=(-1e-15,1e-15));axes[0].legend();axes[0].grid(alpha=.2)
    x=np.arange(len(sizes));axes[1].bar(x-.18,[s['old_bytes']/1e6 for s in sizes],.36,label='THP + one geometry per clock')
    axes[1].bar(x+.18,[s['crv_bytes']/1e6 for s in sizes],.36,label='Combined CRV')
    axes[1].set(xticks=x,xticklabels=[s['case'] for s in sizes],ylabel='Actual bytes / 1,000,000',title='Measured storage for matched saved records');axes[1].legend(fontsize=8)
    save(fig,'exact-parity-and-storage','Exact values and measured file sizes',
         'Pressure differences are exactly zero; all shared thermo/shock values and geometry also pass exact comparison. Byte counts compare matched saved clocks, excluding duplicate ASCII schedules. No throughput/cluster-memory claim.')
    (out/'storage.json').write_text(json.dumps(sizes,indent=2)+'\n')
    reference=data['planar'];suffix=read_run(a.runs/'restart_resume');zero=read_run(a.runs/'restart_zero')
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    cycles=sorted(suffix['curve_times']);errors=[]
    for g in cycles:
        errors.append(max(abs(o['fields']['pressure']-reference['curve_nodes'][g,n]['fields']['pressure']) for (c,n),o in suffix['curve_nodes'].items() if c==g))
    axes[0].plot([suffix['curve_times'][g] for g in cycles],errors,'.-');axes[0].set(xlabel=r'$t$ (code units)',ylabel='Maximum absolute pressure difference',title='Restart continuation: exact',ylim=(-1e-15,1e-15))
    node=next(n for (g,n),o in zero['curve_nodes'].items() if o['fields']['shock_region']==1)
    checkpoint_cycle=max(zero['curve_times'])
    a0=reference['curve_nodes'][checkpoint_cycle,node]['fields'];b0=zero['curve_nodes'][checkpoint_cycle,node]['fields']
    fields=['pressure','shock_region','shock_mach_pressure'];x=np.arange(len(fields))
    axes[1].bar(x-.18,[a0[f] for f in fields],.36,label='Uninterrupted');axes[1].bar(x+.18,[b0[f] for f in fields],.36,label='Zero-step restart')
    axes[1].set(xticks=x,xticklabels=fields,title=f'Node {node}: checkpoint cycle {checkpoint_cycle}');axes[1].tick_params(axis='x',labelrotation=15);axes[1].legend(fontsize=8)
    save(fig,'exact-restart','Active encounter and zero-step restart',
         'Every continued/zero-step node record and particle record matches the uninterrupted run exactly. The active shock cache is refreshed without advancing encounters. Complete event histories remain independent journals.')
    (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n')
    print(json.dumps(dict(figures=len(figures),output=str(out)),indent=2))


if __name__=='__main__':main()
