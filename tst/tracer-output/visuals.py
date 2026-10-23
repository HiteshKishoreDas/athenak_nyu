"""Plot actual native/legacy records from check.py; no reconstructed crossings."""
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
import tracer_binary
import curve_shock_history


def sensor(rows):
    ids={p for _,p in rows}
    return max(ids,key=lambda p:np.ptp([v['fields']['density'] for (g,q),v in rows.items() if q==p]))


def series(rows,p,field):
    selected=sorted((g,v) for (g,q),v in rows.items() if q==p)
    return np.array([v['time'] for _,v in selected]),np.array([v['fields'][field] for _,v in selected])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); run=args.run.resolve(); out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    report=json.loads((run/'report.json').read_text())
    if any(c['status']!='pass' for c in report['checks']):
        raise ValueError('qualification contains failures; inspect report before presenting final figures')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.facecolor':'white','axes.grid':True,'grid.alpha':.18})
    figures=[]
    def save(fig,name,title,caption):
        fig.suptitle(title,fontsize=15)
        fig.tight_layout(rect=(0,.02,1,.95))
        fig.savefig(out/(name+'.png'),dpi=165)
        plt.close(fig)
        figures.append(dict(file=name+'.png',title=title,caption=caption))

    legacy,native=(read_run(run/('planar_'+m)) for m in ('legacy','binary'))
    pid,nid=sensor(native['particles']),sensor(native['curve_nodes'])
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for key,p,color,label in [('particles',pid,'#2677af','particle'),('curve_nodes',nid,'#c24d2d','curve node')]:
        for ax,field in [(axes[0,0],'density'),(axes[0,1],'pressure')]:
            t,v=series(native[key],p,field);lt,lv=series(legacy[key],p,field)
            ax.plot(t,v,color=color,label=label+' binary')
            ax.scatter(lt[::4],lv[::4],s=20,facecolors='none',edgecolors=color,label=label+' legacy')
            ax.set(xlabel='simulation time [code units]',ylabel=field+' [code units]')
        t,v=series(native[key],p,'density');_,lv=series(legacy[key],p,'density')
        axes[1,0].plot(t,v-lv,'.',color=color,label=label)
    axes[0,0].axhline(8/3,color='gray',ls=':',label='RH downstream density')
    axes[0,1].axhline(4.5,color='gray',ls=':',label='RH downstream pressure')
    axes[1,0].set(xlabel='simulation time [code units]',ylabel='binary − legacy density',title='Exact record comparison')
    steps=sorted({s['cycle'] for s in native['geometry']})
    for index,g in enumerate([steps[0],steps[len(steps)//2],steps[-1]]):
        a=next(s for s in native['geometry'] if s['cycle']==g)
        b=next(s for s in legacy['geometry'] if s['cycle']==g)
        x=np.array(a['curves'][0]['nodes'][0]['x'])
        points=np.array([n['x'] for n in a['curves'][0]['nodes']])
        old=np.array([n['x'] for n in b['curves'][0]['nodes']])
        axes[1,1].plot(points[:,0],points[:,1],'-',label=f"binary t={a['time']:.4f}")
        axes[1,1].scatter(old[:,0],old[:,1],s=17,facecolors='none',edgecolors=f'C{index}')
    axes[1,1].set(xlabel='lifted x₁ [code units]',ylabel='lifted x₂ [code units]',title='Same stored curve geometry')
    for ax in axes.flat:ax.legend(fontsize=8)
    save(fig,'01-shared-samples','Thermo and curve geometry: native binary versus legacy',
         f'The existing 128 × 4 translating planar-shock fixture, with 128 continuous particles and 16 curve nodes. '
         f'Particle {pid} and node {nid} are selected by their observed density variation. Lines use stored binary samples; '
         'open markers use legacy samples. Particle THP is unchanged in both modes. The zero difference panel compares every saved density value for these carriers. '
         'Dotted horizontal lines are the independent Rankine–Hugoniot downstream states for the prescribed Mach-2 fixture; '
         'finite numerical shock width remains visible. Geometry points are stored nodes, with edges joining their native connectivity.')

    fig,axes=plt.subplots(3,1,figsize=(12,11),gridspec_kw={'height_ratios':[2.4,1.3,1]})
    rows=[]
    for label,key in [('particle','particle_events'),('node','curve_events')]:
        rows += [(label,e) for e in native[key][:18]]
    for label,key in [('particle','particle_suspended'),('node','curve_suspended')]:
        rows += [(label+' censored',e) for e in native[key][:4]]
    labels=[]
    for j,(label,e) in enumerate(rows):
        color='#8d559c' if 'censored' in label else '#2677af' if label=='particle' else '#c24d2d'
        first,last=e['first']['time'],e['last']['time']
        axes[0].plot([first,last],[j,j],color=color,lw=4)
        axes[0].scatter([first,last],[j,j],s=12,color=color)
        if e['pre']['valid']:axes[0].plot([e['pre']['time'],first],[j,j],color=color,ls=':')
        if e['post']['valid']:axes[0].plot([last,e['post']['time']],[j,j],color=color,ls=':')
        labels.append(f"{label} {e['particle_id'] % (1<<32)} / event {e['sequence']}")
    axes[0].set_yticks(range(len(labels)),labels,fontsize=7)
    axes[0].set(xlabel='simulation time [code units]',title='Recorded inside intervals and entry/exit evidence brackets')
    for label,key,color in [('particles','particle_events','#2677af'),('curve nodes','curve_events','#c24d2d')]:
        e=native[key];old=legacy[key]
        axes[1].scatter([x['first']['time'] for x in e],[x['representative']['mach_pressure'] for x in e],
                        color=color,label=label+' binary',s=30)
        axes[1].scatter([x['first']['time'] for x in old],[x['representative']['mach_pressure'] for x in old],
                        edgecolors='black',facecolors='none',s=60)
    axes[1].axhline(2,ls=':',color='gray',label='prescribed upstream shock-frame Mach = 2')
    axes[1].set(xlabel='first inside sample time [code units]',ylabel='spatial pressure-jump Mach')
    axes[1].legend(fontsize=8)
    trace=[]
    for p in (run/'planar_binary/shock_events/curves').glob('trace*.tbin'):
        row=curve_shock_history.read_trace(p).get(nid)
        if row:trace.append(row)
    trace.sort(key=lambda o:o['time'])
    axes[2].step([o['time'] for o in trace],[int(o['inside']) for o in trace],where='post',label='raw shock-region flag')
    axes[2].scatter([o['time'] for o in trace],[o['phase'] for o in trace],s=11,label='encounter phase (0/1/2)')
    axes[2].set(xlabel='simulation time [code units]',ylabel='flag / phase',yticks=[0,1,2]);axes[2].legend(fontsize=8)
    save(fig,'02-encounters','Encounter histories retain evidence, flags and censored outcomes',
         f"Native output records {len(native['particle_events'])} particle and {len(native['curve_events'])} curve-node emitted events; "
         f"the latest terminal views separately contain {len(native['particle_suspended'])} particle and {len(native['curve_suspended'])} curve censored records. "
         'The raster shows a bounded selection, not unique shock surfaces. Solid intervals run from first to last stored inside evidence; '
         'dotted segments show available outside/inside timing brackets. No exact physical crossing time is inferred. '
         'Black rings are legacy Mach values over the corresponding binary points. The lower panel preserves raw membership separately from accepted-episode phase.')

    amr_old,amr=(read_run(run/('amr_'+m)) for m in ('legacy','binary'))
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    tables=[tracer_binary.read(p) for p in (run/'amr_binary').glob('*.tbin')]
    clocks=sorted({t['cycle']:(t['time'],int(t['metadata']['leaf_count'])) for t in tables}.values())
    axes[0,0].step([t for t,n in clocks],[n for t,n in clocks],where='post')
    axes[0,0].set(xlabel='simulation time [code units]',ylabel='mesh leaves',title='Live refinement and coarsening')
    for d,label,style in [(amr,'binary','-'),(amr_old,'legacy','o')]:
        unique={s['cycle']:s for s in d['geometry']}
        axes[0,1].plot([s['time'] for _,s in sorted(unique.items())],
                       [s['node_count'] for _,s in sorted(unique.items())],style,ms=3,label=label)
    axes[0,1].set(xlabel='simulation time [code units]',ylabel='curve nodes');axes[0,1].legend()
    first={}
    for (g,p),v in sorted(amr['curve_nodes'].items()):first.setdefault(p,v)
    births=[v for v in first.values() if v['birth']>0]
    axes[1,0].scatter([v['birth'] for v in births],[v['time'] for v in births],s=50,facecolors='none',edgecolors='#c24d2d')
    extent=max(v['time'] for v in births)
    axes[1,0].plot([0,extent],[0,extent],':',color='gray')
    axes[1,0].set(xlabel='actual node birth time',ylabel='first saved thermo time',title=f'{len(births)} inserted nodes: no inherited prehistory')
    ref=native['curve_nodes'];node=sensor(ref)
    a=read_run(run/'restart_seed',shock_directory=str(run/'restart-journal'))
    b=read_run(run/'restart_resume',shock_directory=str(run/'restart-journal'))
    t,v=series(ref,node,'pressure');axes[1,1].plot(t,v,label='uninterrupted native')
    for d,label,color in [(a,'before checkpoint','#c24d2d'),(b,'resumed','#2677af')]:
        t,v=series(d['curve_nodes'],node,'pressure');axes[1,1].scatter(t,v,s=12,facecolors='none',edgecolors=color,label=label)
    axes[1,1].set(xlabel='simulation time [code units]',ylabel='node pressure [code units]',title='Restart samples versus uninterrupted run');axes[1,1].legend(fontsize=8)
    save(fig,'03-amr-births-restart','AMR insertion and restart preserve the histories',
         'The existing prescribed translation fixture changes from 4 to 16 leaves, then back to 4. '
         'Binary and legacy node counts overlap. Every newly inserted node begins its saved thermo history at its actual birth under every-step sampling. '
         'The diagonal compares birth time against first saved time, rather than filling a parent trajectory. '
         'The pressure panel uses the separate planar-shock restart case: saved pre-checkpoint and resumed samples overlap uninterrupted output exactly. '
         'The test also compares complete event records and checks that zero-step replay leaves durable journal files unchanged.')

    fig,axes=plt.subplots(1,2,figsize=(12,5))
    categories=['particle thermo','curve thermo','geometry','events','optional traces']
    values={}
    for mode in ('legacy','binary'):
        d=run/('planar_'+mode);sums=dict.fromkeys(categories,0)
        for p in d.rglob('*'):
            if not p.is_file():continue
            if p.suffix=='.thp':key='particle thermo'
            elif 'curve_thermo_history' in p.parts and p.suffix in ('.csv','.tbin'):key='curve thermo'
            elif 'shock_events' in p.parts and p.suffix in ('.csv','.tbin'):
                key='optional traces' if 'trace' in p.name else 'events'
            elif p.suffix in ('.txt','.tbin') and ('curves' in p.parts or '.curves.' in p.name):key='geometry'
            else:continue
            sums[key]+=p.stat().st_size
        values[mode]=sums
    bottoms=np.zeros(2)
    for key in categories:
        amounts=np.array([values[m][key] for m in ('legacy','binary')])/1e6
        axes[0].bar(['legacy','binary'],amounts,bottom=bottoms,label=key);bottoms+=amounts
    axes[0].set(ylabel='actual output bytes [MB, decimal]',title='Measured storage for this small fixture');axes[0].legend(fontsize=8)
    durations={r['name']:r['seconds'] for r in report['runs']}
    axes[1].bar(['legacy','binary'],[durations['planar_'+m] for m in ('legacy','binary')],color=['#777','#2677af'])
    axes[1].set(ylabel='one native run wall time [seconds]',title='One local trial; not a performance qualification')
    save(fig,'04-measured-cost','Output size and local run cost: measurements, not a budget pass',
         'Bytes come from actual files, including unchanged particle THP output. Native tables repeat named schemas and metadata in each immutable chunk. '
         'This short, sparsely emitting shock fixture therefore produces larger binary output; binary does not guarantee smaller files. '
         'The runtime bars are single observations, with possible concurrent laptop activity and no uncertainty estimate. '
         'No RUSTY, Frontier, MPI, GPU, production scaling, memory budget or scientific-performance approval is claimed.')

    cadence=read_run(run/'cadence')
    fig,ax=plt.subplots(figsize=(12,3.7))
    for y,key,label in [(0,'particle_times','particle thermo'),(1,'curve_times','curve thermo')]:
        ax.scatter(sorted(cadence[key]),[y]*len(cadence[key]),s=65,label=label)
    cycles=sorted({s['cycle'] for s in cadence['geometry']})
    ax.scatter(cycles,[2]*len(cycles),s=65,label='geometry')
    ax.set(yticks=[0,1,2],yticklabels=['particle thermo','curve thermo','curve geometry'],
           xlabel='actual completed fluid cycle',xticks=range(11),ylim=(-.5,2.5))
    save(fig,'05-independent-clocks','Same output machinery, explicit independent sampling clocks',
         'Thermo output is requested every two cycles and geometry every three through ordinary output blocks. '
         'Both also receive the native terminal sample at cycle 10. Markers show actual saved records; no samples are invented on intervening cycles. '
         'Shock detection and durable event draining continue every completed timestep, independently of these diagnostic output clocks.')

    if (run/'amr_shock_binary').is_dir():
        coupled=read_run(run/'amr_shock_binary')
        observations={}
        for pattern in ('trace*.tbin','birthtrace*.tbin'):
            for path in (run/'amr_shock_binary/shock_events/curves').glob(pattern):
                for node,o in curve_shock_history.read_trace(path).items():
                    if node not in observations or o['time']<observations[node]['time']:
                        observations[node]=o
        born={node:o for node,o in observations.items() if o['birth']>0}
        inside={node:o for node,o in born.items() if o['inside']}
        fig,axes=plt.subplots(1,2,figsize=(12,5))
        for condition,color,label in [(True,'#c24d2d','born inside: entry unknown'),
                                      (False,'#2677af','born outside: clean start')]:
            selected=[(node,o) for node,o in born.items() if o['inside']==condition]
            labels={node:next(v['a'] for (g,q),v in coupled['curve_nodes'].items() if q==node)
                    for node,o in selected}
            axes[0].scatter([o['birth'] for node,o in selected],[labels[node] for node,o in selected],
                            color=color,s=35,label=label)
        axes[0].set(xlabel='actual insertion time [code units]',ylabel='immutable material label a',
                    title='Inside and outside births during live shock AMR');axes[0].legend(fontsize=8)
        node=min(inside);birth=inside[node]['birth']
        t,v=series(coupled['curve_nodes'],node,'pressure')
        axes[1].plot(t,v,label='uninterrupted binary')
        journal=str(run/'amr-shock-restart-journal')
        for folder,label,color in [('amr_shock_seed','before checkpoint','#c24d2d'),
                                   ('amr_shock_resume','resumed','#2677af')]:
            d=read_run(run/folder,shock_directory=journal)
            t,v=series(d['curve_nodes'],node,'pressure')
            axes[1].scatter(t,v,s=13,facecolors='none',edgecolors=color,label=label)
        axes[1].axvline(birth,color='gray',ls=':',label='actual birth')
        axes[1].axvspan(0,birth,color='gray',alpha=.1)
        axes[1].set(xlabel='simulation time [code units]',ylabel='newborn node pressure [code units]',
                    xlim=(0,max(coupled['curve_times'].values())),title='An inside-born node survives AMR and restart')
        axes[1].legend(fontsize=8)
        save(fig,'06-coupled-amr-shock-restart','Live AMR, shock encounters and curve thermo work together',
             f'The existing CSH live-shock AMR fixture creates {len(inside)} nodes inside the detected region and '
             f'{len(born)-len(inside)} outside it at positive time. Each earliest observation has sequence zero and actual time equal to birth; '
             'inside-born nodes retain the incomplete-entry flag. The points show the actual insertion times and immutable material labels. '
             'The pressure history starts at birth; the shaded pre-birth interval contains no fabricated parent samples. '
             'A checkpoint at an actual inside birth resumes both continuous particles and curves. Thermo, geometry, events and terminal views '
             'match the uninterrupted run exactly. This fixture records births in the normal completed-step trace; it does not manufacture separate birthtrace files. '
             'This is an interaction/persistence check, not AMR shock-strength convergence or scientific sign-off.')

    (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    body=''.join(f'<section><h2>{html.escape(f["title"])}</h2><p>{html.escape(f["caption"])}</p>'
                 f'<img src="{f["file"]}" alt="{html.escape(f["title"])}"></section>' for f in figures)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>AthenaK native tracer outputs</title>'
        '<style>body{font:17px/1.5 system-ui;max-width:1250px;margin:35px auto;padding:0 22px;background:#f4f5f7;color:#17212b}'
        'section{background:white;border-radius:12px;padding:22px;margin:24px 0}img{width:100%;height:auto}p{max-width:1100px}</style>'
        '<h1>AthenaK: native particle, curve and shock outputs</h1>'
        f'<p>{len(report["checks"])} output checks passed on local Serial CPU/double. Physical signals below come from saved native files. '
        'Human scientific sign-off and cluster qualification remain pending.</p>'+body)
    print(out/'index.html')


if __name__=='__main__':main()
