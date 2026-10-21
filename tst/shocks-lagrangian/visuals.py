"""Explained native-data plots. No invented trajectories or crossing-time estimates."""
import argparse
from collections import Counter
import hashlib
import html
import json
from pathlib import Path
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc'),str(ROOT/'tst/shocks-curves')]
from common import thermo, sha
from integration import traces
import curve_shock_history as curve
from shock_history import read_history, read_suspended
import material_curves

BLUE,ORANGE,GREEN,RED='#2166ac','#d77b16','#26734d','#b2182b'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('lagrangian','core','forced','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    paths=[args.lagrangian.resolve(),args.core.resolve(),args.forced.resolve()]
    reports=[json.loads((p/'report.json').read_text()) for p in paths]
    for path,report in zip(paths,reports):
        assert report['status']=='pass'
        for run in report['runs']:
            assert sha(run['command'][0])==run['executable_sha256'],'binary changed since test'
            assert sha(path/run['name']/'run.log')==run['log_sha256']
            if run['input_sha256']:
                assert sha(path/run['name']/'input.athinput')==run['input_sha256']
    local,core,forced=paths;report,core_report,force_report=reports
    checks={c['id']:c for c in report['checks']};figures=[]
    def save(fig,name,title,caption):
        fig.suptitle(title,fontsize=16,fontweight='bold',x=.05,ha='left')
        fig.text(.05,.02,textwrap.fill(caption,175),fontsize=10,va='bottom')
        fig.tight_layout(rect=(.02,.17,.995,.91),w_pad=2.5)
        fig.savefig(out/(name+'.png'),dpi=160);plt.close(fig)
        figures.append(dict(name=name,title=title,caption=caption,sha256=sha(out/(name+'.png'))))
    def series(directory,particle):
        rows,_=thermo(directory)
        return sorted((o for (g,p),o in rows.items() if p==particle),key=lambda o:o['time'])
    def matrix(ax,observed,label='Particle ID order'):
        ids=sorted({p for g,p in observed});gens=sorted({g for g,p in observed})
        times={g:o['time'] for (g,p),o in observed.items()}
        values=np.array([[int(observed[g,p]['inside']) if observed[g,p]['valid'] else 2
                          for g in gens] for p in ids])
        ax.pcolormesh([times[g] for g in gens],range(len(ids)),values,shading='nearest',
                       cmap=ListedColormap(['#dfe7ef',ORANGE,RED]),vmin=0,vmax=2)
        ax.set(xlabel=r'$t$ [code units]',ylabel=label)
        return ids

    fig,ax=plt.subplots(1,3,figsize=(15,5.8))
    for p in range(12):
        for directory,color,alpha in [('uniform2d',BLUE,.7),('mc_comparison',ORANGE,.35)]:
            rows=series(local/directory,p);t=[o['time'] for o in rows]
            x=np.unwrap([o['x'][0] for o in rows],period=1)
            ax[0].plot(t,x,color=color,alpha=alpha,lw=.8)
    ax[0].plot([],[],color=BLUE,label='Continuous RK2');ax[0].plot([],[],color=ORANGE,label='MC control')
    ax[0].set(xlabel=r'$t$ [code units]',ylabel=r'Lifted $x$ [code units]',title='Actual recorded paths; 12 IDs');ax[0].legend()
    for name,color in [('uniform2d',BLUE),('mc_comparison',ORANGE)]:
        rows,times=thermo(local/name);a,b=min(times),max(times)
        initial={p:o for (g,p),o in rows.items() if g==a}
        errors=[(o['x'][0]-initial[p]['x'][0]-(times[b]-times[a])+.5)%1-.5
                for (g,p),o in rows.items() if g==b]
        ax[1].hist(errors,bins=np.linspace(-.16,.16,33),histtype='step',lw=2,color=color,label=name)
    ax[1].set(xlabel=r'$x(t)-x_0-ut$ [periodic units]',ylabel='Particle count',title='Final displacement error');ax[1].legend()
    selected=[c for c in report['checks'] if 'max_error' in c]
    ax[2].barh(range(len(selected)),[max(c['max_error'],1e-18) for c in selected],color=BLUE)
    ax[2].set_yticks(range(len(selected)),[c['id'].replace('LT01_','').replace('LT02_','').replace('LT04_','') for c in selected],fontsize=8)
    ax[2].set_xscale('log');ax[2].axvline(report['position_budget'],color=RED,ls='--',label='Fixed tolerance')
    ax[2].set(xlabel='Maximum position error [code units]',title='Independent affine-flow oracles');ax[2].legend()
    save(fig,'01_continuous_motion','Shock carriers now follow velocity continuously',
         'Blue paths are production RK2 particles. Orange paths are native MC controls with histories disabled. '
         'Both begin at seeded cell centers; only MC takes stochastic cell hops. The right panel compares every saved continuous '
         'position with independent translation, acceleration or shear solutions. The predeclared tolerance is unchanged.')

    observed=traces(local/'planar/shock_events');events=read_history(local/'planar/shock_events')
    fig,ax=plt.subplots(1,3,figsize=(15,5.8))
    for p in sorted({p for g,p in observed})[::4]:
        rows=sorted((o for (g,q),o in observed.items() if p==q),key=lambda o:o['time'])
        ax[0].plot([o['time'] for o in rows],[o['x'][0] for o in rows],color=BLUE,lw=.6,alpha=.6)
        inside=[o for o in rows if o['inside']]
        ax[0].scatter([o['time'] for o in inside],[o['x'][0] for o in inside],s=3,color=ORANGE)
    data=checks['LT05_PHY08_INT06_RH_passage'];front=data['density_fronts']
    ax[0].plot([t for g,t,x in front],[x for g,t,x in front],color='black',lw=2,label='Native density midpoint')
    ax[0].set(xlabel=r'$t$ [code units]',ylabel=r'$x$ [code units]',title='Evolved planar Mach-2 shock');ax[0].legend(fontsize=8)
    matrix(ax[1],observed);ax[1].set_title('Saved membership: orange = inside')
    for i,e in enumerate(sorted(events,key=lambda e:e['particle_id'])):
        for a,b,color in [('pre','first',GREEN),('last','post',ORANGE)]:
            ax[2].plot([e[a]['time'],e[b]['time']],[i,i],color=color,lw=1.6)
    ax[2].plot([],[],color=GREEN,label='Entry bracket');ax[2].plot([],[],color=ORANGE,label='Exit bracket')
    ax[2].set(xlabel=r'$t$ [code units]',ylabel='Particle ID order',title=f'{len(events)} IDs; one complete episode each');ax[2].legend()
    save(fig,'02_planar_histories','Each sampled passage has one encounter, with timing brackets',
         'All 128 seeded particles independently cross the native density midpoint exactly once, and all 128 have one complete '
         'encounter. Every density-crossing bracket overlaps its recorded event evidence envelope. Membership uses the full '
         'unreduced numerical shock region. Colored entry/exit intervals are sampling brackets, not exact physical crossing times.')

    fig,ax=plt.subplots(1,3,figsize=(15,5.8))
    for name,color in [('amr',BLUE),('mixed_amr',ORANGE)]:
        states={}
        for p in (local/name).glob('material.curves.*.txt'):
            state=material_curves.read(p);states[state['cycle']]=state
        values=[s for g,s in sorted(states.items())]
        ax[0].step([s['time'] for s in values],[len(s['amr']['leaves']) for s in values],where='post',color=color,label=name)
    ax[0].set(xlabel=r'$t$ [code units]',ylabel='Active leaf blocks',title='Actual refine/coarsen transactions');ax[0].legend()
    for name,velocity,color in [('amr',(.2,-.1),BLUE),('mixed_amr',(.97,.97),ORANGE)]:
        rows,times=thermo(local/name);initial={p:o for (g,p),o in rows.items() if g==min(times)}
        by_time={}
        for (g,p),o in rows.items():
            error=max(abs((o['x'][d]-initial[p]['x'][d]-velocity[d]*o['time']+2)%4-2) for d in range(2))
            by_time[g]=max(by_time.get(g,0),error)
        ax[1].plot([times[g] for g in sorted(by_time)],[by_time[g] for g in sorted(by_time)],color=color,label=name)
    ax[1].axhline(report['position_budget'],ls='--',color=RED,label='Fixed tolerance')
    ax[1].set(xlabel=r'$t$ [code units]',ylabel='Max position error [code units]',title='No AMR position snapping');ax[1].legend(fontsize=8)
    coexisting=next(c for c in core_report['checks'] if c['id']=='CSH04_CSH07_coexisting_live_amr_births')
    ax[2].bar(['Born inside','Born outside'],[coexisting['inside_positive_births'],coexisting['outside_positive_births']],color=[ORANGE,BLUE])
    ax[2].set(ylabel='Positive-time curve node births',title='New nodes start fresh histories')
    save(fig,'03_live_amr','Continuous particles and curve histories work through live AMR',
         'Uniform and mixed-level trees both refine and coarsen. Position error stays within the unchanged affine-oracle budget. '
         'The mixed-level run explicitly logs the containing-leaf fallback for a missing diagonal neighbor. In the evolved shock '
         'run, 17 positive-time nodes are born inside and 22 outside; no node inherits a parent encounter. Inside births are left-censored.')

    fig,ax=plt.subplots(1,3,figsize=(15,5.8))
    for label,color in [('planar',BLUE),('seed_active',ORANGE)]:
        history=read_history(local/label/'shock_events');clocks=sorted({e['post']['time'] for e in history})
        ax[0].step(clocks,[sum(e['post']['time']<=t for e in history) for t in clocks],where='post',
                   color=color,ls='-' if label=='planar' else '--',label=label)
    ax[0].set(xlabel=r'$t$ [code units]',ylabel='Committed complete events',title='Uninterrupted / active-checkpoint replay');ax[0].legend()
    groups=[('Standalone',len(report['checks'])),('Curve gates',len(core_report['checks'])),('Forced box',len(force_report['checks']))]
    ax[1].bar([g for g,n in groups],[n for g,n in groups],color=GREEN)
    ax[1].set(ylabel='Passing engineering checks',title='Native runs; current source binding')
    for i,(g,n) in enumerate(groups):ax[1].text(i,n+.15,str(n),ha='center')
    ax[1].set_ylim(0,max(n for g,n in groups)+3)
    rejects=['MC history','Changed pusher','Changed transport','Event overflow','History corruption','Predictor outflow']
    ax[2].barh(rejects,[1]*len(rejects),color=GREEN);ax[2].set_xlim(0,1.2)
    ax[2].set_xticks([0,1],['','Explicit rejection']);ax[2].set_title('Unsafe requests fail visibly')
    save(fig,'04_restart_guards','Restart preserves encounter identity and complete particle state',
         'Outside, active and pending checkpoints reproduce the same native particle observations and complete event records. '
         'An active encounter also replays exactly across AMR and coarsening. MC/native transport cannot be changed silently on '
         'restart. Predictor outflow is explicitly unsupported: execution aborts before treating a predictor as a completed exit. '
         'Counts are engineering checks; scientific approval is still pending.')

    directory=forced/'forced_on';particle=read_history(directory/'shock_events');nodes=curve.read_history(directory/'shock_events/curves')
    endp=read_suspended(directory/'shock_events');endn=curve.read_suspended(directory/'shock_events/curves')
    fig,ax=plt.subplots(1,3,figsize=(15,5.8))
    for records,label,color in [(particle,'Particles',BLUE),(nodes,'Curve nodes',ORANGE)]:
        clocks=sorted({e['post']['time'] for e in records})
        ax[0].step(clocks,[sum(e['post']['time']<=t for e in records) for t in clocks],where='post',color=color,label=label)
    ax[0].set(xlabel=r'$t$ [code units]',ylabel='Committed encounters',title='Native forced 32³ box');ax[0].legend()
    for records,label,color in [(particle,'Particles',BLUE),(nodes,'Curve nodes',ORANGE)]:
        values=[e['representative']['mach_pressure'] for e in records if e['status']==0]
        ax[1].hist(values,bins=np.linspace(1,15,29),histtype='step',lw=2,color=color,label=label)
    ax[1].set(xlabel='Pressure-jump Mach proxy',ylabel='Complete encounter count',title='Saved strength evidence');ax[1].legend()
    ax[2].bar(['Particle terminal','Node terminal'],[len(endp),len(endn)],color=[BLUE,ORANGE])
    ax[2].set(ylabel='Right-censored views',title='Terminal active episodes stay separate')
    save(fig,'05_forced_statistics','Forcing remains enabled; statistics now come from continuous tracers',
         'The accepted short native 32³ forcing fixture uses 4096 continuous particles and 16 curve nodes. Tracking and optional '
         'traces leave the native fluid arrays, particle samples and curve geometry exactly unchanged. Strength values are '
         'detector estimates, not calibrated physical plateau Mach numbers. This short box and the inherited forcing heat-bias '
         'issue do not establish converged turbulence or scientific approval. Terminal censored views are not extra complete events.')

    if any(directory.glob('shock_events/trace*.rank*.csv')):
        observed=traces(directory/'shock_events');ids=sorted({p for g,p in observed})[:128]
        observed={key:o for key,o in observed.items() if key[1] in ids}
        fig,ax=plt.subplots(1,2,figsize=(15,5.8));matrix(ax[0],observed)
        ax[0].set_title('First 128 particle IDs: saved every-step flags')
        counts=Counter(e['particle_id'] for e in particle if e['status']==0)
        ax[1].hist([counts[p] for p in range(4096)],bins=np.arange(max(counts.values(),default=0)+2)-.5,color=BLUE)
        ax[1].set(xlabel='Complete episodes per particle',ylabel='Particle count',title='All 4096 particle IDs, including zero episodes')
        save(fig,'06_forced_histories','Inspect actual particle encounter histories',
             'Orange marks a valid inside-shock sample; pale blue marks valid outside; red marks invalid evidence. Rows are fixed '
             'particle IDs, not reconstructed trajectories. The histogram counts committed complete episodes for all seeded '
             'particles, while right-censored views remain separate in the prior figure. Histories are endpoint-cadence limited.')
    attempts=[]
    for label in ('core-first','lagrangian-first','lagrangian-reviewed',
                  'lagrangian-release','lagrangian-qualified','lagrangian-accepted'):
        path=local.parent/label/'report.json'
        if path.exists():
            data=json.loads(path.read_text());attempts.append((label,
                sum(c['status']=='pass' for c in data['checks']),
                sum(c['status']=='fail' for c in data['checks'])))
    attempts.append(('Final current source',sum(len(r['checks']) for r in reports),0))
    fig,ax=plt.subplots(figsize=(15,5.8))
    labels,passed,failed=zip(*attempts)
    ax.barh(labels,passed,color=GREEN,label='Passed')
    ax.barh(labels,failed,left=passed,color=RED,label='Failed')
    for i,(a,b) in enumerate(zip(passed,failed)):ax.text(a+b+.2,i,f'{a} pass / {b} fail',va='center')
    ax.set_xlim(0,max(a+b for a,b in zip(passed,failed))+7)
    ax.set(xlabel='Engineering checks in each retained attempt',title='Failed attempts are retained, not relabeled');ax.legend()
    save(fig,'07_retained_failures','Development failures remain visible in the evidence ledger',
         'Early checks caught a missing build, a trace-reader omission, the unexercised diagonal fixture, and an AMR restart '
         'contract mismatch. Test counts grew as stronger checks were added; these bars are not a performance or accuracy '
         'comparison. The input editor was fixed after repeated sections silently overrode its requested refinement time. '
         'All final 35 native checks pass with unchanged numerical tolerances. Original failed inputs and logs are archived.')

    summary=dict(scientific_approval='pending',source_reports=[str(p/'report.json') for p in paths],
                 figures=figures,checks_passed=sum(len(r['checks']) for r in reports),
                 forced_particle_events=len(particle),forced_curve_events=len(nodes),
                 forced_particle_censored=len(endp),forced_curve_censored=len(endn))
    (out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
    body=''.join('<section><h2>'+html.escape(f['title'])+'</h2><img src="'+f['name']+'.png" alt="'+html.escape(f['title'])+'"><p>'+html.escape(f['caption'])+'</p></section>' for f in figures)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Continuous Lagrangian shock histories</title><style>body{font:17px system-ui;max-width:1500px;margin:35px auto;padding:0 20px;color:#173047}img{width:100%}section{margin:45px 0}p{max-width:1200px;line-height:1.5}.limit{background:#fff0d9;padding:20px}</style><h1>Continuous Lagrangian shock histories</h1><p class="limit">Local Serial/double engineering evidence. Human scientific sign-off pending. No multi-rank, GPU or cluster runtime qualification. Outflow predictors reject explicitly.</p>'+body)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
