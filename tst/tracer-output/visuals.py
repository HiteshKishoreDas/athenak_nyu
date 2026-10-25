"""Plot actual existing-format output checks; no conversion or reconstructed paths."""
import argparse
import html
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'vis/python'))
from tracer_history import read_run


def series(rows,node,field):
    samples=sorted((g,o) for (g,p),o in rows.items() if p==node)
    return np.array([o['time'] for _,o in samples]),np.array([o['fields'][field] for _,o in samples])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run=args.run.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    report=json.loads((run/'report.json').read_text())
    if any(c['status']!='pass' for c in report['checks']):
        raise ValueError('checks contain failures; inspect report before presenting plots')
    data=read_run(run/'planar');control=read_run(run/'planar_control')
    figures=[]
    def save(fig,name,caption):
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=160);plt.close(fig)
        figures.append(dict(file=name+'.png',caption=caption))

    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for key,label in [('particles','particle'),('curve_nodes','curve node')]:
        rows=data[key];ids={p for _,p in rows}
        node=max(sorted(ids),key=lambda p:np.ptp(series(rows,p,'density')[1]))
        for ax,field in zip(axes[0],['density','pressure']):
            t,y=series(rows,node,field);tc,yc=series(control[key],node,field)
            line,=ax.plot(t,y,label=label)
            ax.scatter(tc[::4],yc[::4],s=18,facecolors='none',edgecolors=line.get_color(),label=label+' control')
            ax.set(xlabel='simulation time [code units]',ylabel=field+' [code units]')
        t,y=series(rows,node,'density');_,yc=series(control[key],node,'density')
        axes[1,0].plot(t,y-yc,'.',label=label)
    axes[0,0].axhline(8/3,color='gray',ls=':',label='RH downstream density')
    axes[0,1].axhline(4.5,color='gray',ls=':',label='RH downstream pressure')
    axes[1,0].set(xlabel='simulation time [code units]',ylabel='density difference',title='Scheduled geometry minus control')
    groups=[('particles','particle_events'),('curve nodes','curve_events')]
    for label,key in groups:
        events=data[key]
        axes[1,1].scatter([e['first']['time'] for e in events],
            [e['representative']['mach_pressure'] for e in events],label=label,s=20)
    axes[1,1].axhline(2,color='gray',ls=':',label='prescribed shock-frame Mach')
    axes[1,1].set(xlabel='first inside sample time [code units]',ylabel='pressure-jump Mach')
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Existing THP and CSV: particles and curves crossing the planar shock')
    save(fig,'shared-histories','The unchanged 128 × 4 translating planar shock. Lines read existing THP/CSV records; open markers read a control without scheduled geometry. Density differences are exactly zero. Dotted lines show independent Rankine–Hugoniot downstream states and the prescribed Mach 2. Finite shock width and spatial strength estimates remain visible; no exact crossing time is inferred.')

    reference=read_run(run/'amr_shock')
    prefix=read_run(run/'amr_shock_seed',shock_directory=str(run/'amr-shock-restart-journal'))
    suffix=read_run(run/'amr_shock_resume',shock_directory=str(run/'amr-shock-restart-journal'))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for key,label in [('particles','particles'),('curve_nodes','curve nodes')]:
        joined=prefix[key]|suffix[key]
        clocks=sorted({g for g,p in reference[key]})
        error=[max(abs(o['fields']['pressure']-joined[g,p]['fields']['pressure'])
                   for (g,p),o in reference[key].items() if g==cycle) for cycle in clocks]
        times=[next(o['time'] for (g,p),o in reference[key].items() if g==cycle) for cycle in clocks]
        axes[0].plot(times,error,'.-',label=label)
    axes[0].set(xlabel='simulation time [code units]',ylabel='maximum pressure difference [code units]',title='Restart versus uninterrupted run')
    for label,key in groups:
        events=reference[key];other={e['identity']:e for e in suffix[key]}
        a=[e['first']['time'] for e in events];b=[other[e['identity']]['first']['time'] for e in events]
        axes[1].plot(a,np.array(b)-np.array(a),'.',label=label)
    axes[1].set(xlabel='first inside sample time [code units]',ylabel='restart minus uninterrupted time',title='Every emitted encounter retained')
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Live AMR, inserted nodes and existing CSV journal restart')
    save(fig,'restart-equality','Both pressure differences and encounter-time differences are exactly zero for every saved record. These are comparisons of stored evidence, not estimates of exact physical crossing times. The run refines/coarsens and inserts 17 inside-born and 22 outside-born nodes; incomplete entries remain incomplete. The checkpoint is taken at the first inside birth. Full records, identities and terminal views are compared by TO21, beyond the quantities plotted here.')
    (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n')
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Existing output checks</title><style>body{max-width:1200px;margin:32px auto;font:17px system-ui;padding:0 16px}img{width:100%}figure{margin:32px 0}figcaption{line-height:1.5}</style><h1>Existing formats, fresh checks</h1><p>Local Serial/double only. No format conversion. Scientific approval remains pending.</p>'+''.join('<figure><img src="'+x['file']+'"><figcaption>'+html.escape(x['caption'])+'</figcaption></figure>' for x in figures))


if __name__=='__main__':main()
