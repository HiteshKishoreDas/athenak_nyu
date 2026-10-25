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
import material_curves


def series(rows,node,field):
    samples=sorted((g,o) for (g,p),o in rows.items() if p==node)
    return np.array([o['time'] for _,o in samples]),np.array([o['fields'][field] for _,o in samples])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--id-boundaries',type=Path,help='retained id_check.py output')
    parser.add_argument('--variable-layout-run',type=Path,help='variable_id_check.py results')
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
    if args.id_boundaries:
        boundary=json.loads((args.id_boundaries/'report.json').read_text())
        if any(c['status']!='pass' for c in boundary['checks']):
            raise ValueError('ID boundary checks contain failures')
        state=material_curves.read(boundary['last_id_state'])
        curve=state['curves'][0];nodes=curve['nodes'];curve_bits=state['id_layout'][1]
        fig,axes=plt.subplots(2,1,figsize=(12,8),height_ratios=[1.2,1])
        axes[0].plot([n['a'] for n in nodes],[n['x'][0] for n in nodes],'o-')
        for n in nodes:
            sequence=n['id']&((1<<(32-curve_bits))-1)
            axes[0].annotate(f"sequence {sequence}\nID {n['id']}",
                (n['a'],n['x'][0]),xytext=(0,12),textcoords='offset points',ha='center',fontsize=9)
        axes[0].margins(x=.15,y=.35)
        axes[0].set(xlabel=r'immutable material label $a$ [code units]',
            ylabel=r'lifted $x_1$ [code units]',title='Native insertion: links follow material order, IDs retain creation order')
        bits=np.array([[(n['id']>>(31-bit))&1 for bit in range(32)] for n in nodes])
        axes[1].imshow(bits,cmap='Blues',vmin=0,vmax=1,aspect='auto',interpolation='nearest')
        if curve_bits:axes[1].axvline(curve_bits-.5,color='orange',lw=3)
        axes[1].set_yticks(range(len(nodes)),[str(n['id']) for n in nodes])
        axes[1].set_xticks(range(32),range(31,-1,-1),fontsize=8)
        axes[1].set(xlabel='bit position (blue = 1, white = 0)',ylabel='exact unsigned node ID',
            title=f"{curve_bits} curve bits (logical curve {curve['id']}) | {32-curve_bits} creation-sequence bits")
        fig.suptitle('32-bit curve IDs: maximum node ID and consistent neighbour links')
        save(fig,'compact-id-boundary','A synthetic counter-boundary test on the accepted uniform-flow fixture: the allocation counter is deliberately advanced to its last available sequence. Native insertion creates node 4294967295 between earlier nodes, with reciprocal prev/next links verified by the production reader. The upper line connects actual saved nodes in material order; it is not a reconstructed trajectory. The lower panel shows the exact 32 stored ID bits. Separate tests verify exhausted counters and over-capacity split batches never wrap or reuse IDs. This does not simulate millions of nodes or qualify their performance.')
    if args.variable_layout_run:
        variable=json.loads((args.variable_layout_run/'report.json').read_text())
        if any(c['status']!='pass' for c in variable['checks']):
            raise ValueError('variable layout checks contain failures')
        actual=next(c for c in variable['checks'] if c['id'].startswith('VB01'))
        bits=list(range(31))
        fig,axes=plt.subplots(2,1,figsize=(12,8),height_ratios=[1.4,1])
        axes[0].semilogy(bits,[max(1,(1<<b)-1) for b in bits],label='curve ID capacity')
        axes[0].semilogy(bits,[(1<<(32-b))-1 for b in bits],label='node creations per curve')
        tested=[0,1,4,8,16,30]
        axes[0].scatter(tested,[(1<<(32-b))-1 for b in tested],label='native layout checks',zorder=3)
        axes[0].axvline(8,color='gray',ls=':',label='unchanged default')
        axes[0].set(xlabel='input curve_id_bits',ylabel='identity capacity (log scale)',
            title='Mathematical ID capacities; not a memory or throughput qualification')
        axes[0].grid(alpha=.2);axes[0].legend()
        rows=[[b,32-b,actual[str(b)]['curve_id'],', '.join(map(str,actual[str(b)]['node_ids']))] for b in tested]
        axes[1].axis('off')
        table=axes[1].table(cellText=rows,colLabels=['curve bits','sequence bits','logical curve ID','actual stored node IDs'],
            colWidths=[.13,.16,.22,.49],loc='center')
        table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,1.8)
        axes[1].set_title('Actual records from six short native constant-flow checks')
        fig.suptitle('One optional input: curve bits versus node sequence capacity')
        save(fig,'configurable-id-capacity','Solid curves show the mathematical unsigned-32 identity namespace with zero reserved. With zero curve bits, logical curve ID 1 lives in metadata and all 32 ID bits are sequence. The table contains actual saved IDs from six native checks. Defaults remain 8/24. Separate native tests cover zero-bit real shock histories, live AMR births, exact restart and split-change rejection. These capacities do not qualify billions of resident nodes: max_nodes, MPI batch counts and curve IO caps still apply.')
    (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n')
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Existing output checks</title><style>body{max-width:1200px;margin:32px auto;font:17px system-ui;padding:0 16px}img{width:100%}figure{margin:32px 0}figcaption{line-height:1.5}</style><h1>Existing formats, fresh checks</h1><p>Local Serial/double only. No format conversion. Scientific approval remains pending.</p>'+''.join('<figure><img src="'+x['file']+'"><figcaption>'+html.escape(x['caption'])+'</figcaption></figure>' for x in figures))


if __name__=='__main__':main()
