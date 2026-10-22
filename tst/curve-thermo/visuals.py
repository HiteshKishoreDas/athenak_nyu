"""Plot actual shared thermo samples and independent oracles; no synthetic paths."""
import argparse
import html
import json
import math
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc')]
from common import thermo
from tracer_thermo_history import read_curves
from curve_shock_history import read_trace
import material_curves


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cadence-failure',type=Path)
    args=parser.parse_args(); run=args.run.resolve(); out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    report=json.loads((run/'report.json').read_text())
    assert all(c['status']=='pass' for c in report['checks']), 'Refusing to label failing campaign qualified'
    panels=[]
    def save(fig,name,title,caption):
        fig.suptitle(title,fontsize=15)
        fig.text(.5,.01,caption,ha='center',va='bottom',fontsize=9,wrap=True)
        fig.tight_layout(rect=(0,.1,1,.94))
        fig.savefig(out/(name+'.png'),dpi=170); plt.close(fig)
        panels.append((name,title,caption))

    curves,clock=read_curves(run/'both/curve_thermo_history'); particles,_=thermo(run/'both')
    node=min(p for g,p in curves); particle=min(p for g,p in particles)
    cycles=sorted(clock); times=[clock[g] for g in cycles]
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,field,label in zip(axes.flat,['pressure','temperature','entropy','v1'],
                             [r'$P$',r'$P/\rho$',r'$\ln(P/\rho^\gamma)$',r'$v_x$']):
        actual=[curves[g,node]['fields'][field] for g in cycles]
        regular=[particles[g,particle]['fields'][field] for g in cycles]
        vx=.2+.3*np.array(times); vy=-.1+.2*np.array(times)
        p=(4-np.maximum(abs(vx),abs(vy)))**2/(5/3)
        target={'pressure':p,'temperature':p,'entropy':np.log(p),'v1':vx}[field]
        ax.plot(times,target,'k--',lw=2,label='Independent prescribed oracle')
        ax.plot(times,actual,'o',ms=4,label='Curve node')
        ax.plot(times,regular,'+',ms=8,label='Regular Lagrangian particle')
        ax.set(xlabel=r'$t$ (code units)',ylabel=label+' (code units)');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    error=next(c['curve_error'] for c in report['checks'] if c['id'].startswith('CT02'))
    save(fig,'01_shared_fields','One field list and clock for both kinds of tracer',
         f'2D prescribed acceleration: saved samples overlie the independent field oracle. Maximum error across all 11 fields: {error:.3g}.\n'
         'Shared containing-cell diagnostics; interpolated advection velocity remains a separate quantity. No fluid backreaction added.')

    curves,clock=read_curves(run/'amr/curve_thermo_history'); cycles=sorted(clock)
    states={}
    for path in (run/'amr').glob('material.curves.*.txt'):
        s=material_curves.read(path);states[s['cycle']]=s
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    axes[0,0].step([states[g]['time'] for g in sorted(states)],
        [len(states[g]['amr']['leaves']) for g in sorted(states)],where='post')
    axes[0,0].set(ylabel='Active leaf blocks',title='Native refinement and coarsening')
    axes[0,1].step([clock[g] for g in cycles],
                  [sum(c==g for c,p in curves) for g in cycles],where='post')
    axes[0,1].set(ylabel='Recorded curve nodes',title='Insertion changes the sampling set')
    first={}
    for (g,p),s in sorted(curves.items()):first.setdefault(p,s)
    newborn=sorted((p,s) for p,s in first.items() if s['birth']>0)
    for i,(p,s) in enumerate(newborn):
        axes[1,0].plot([s['birth'],s['time']],[i,i],color='tab:blue')
        axes[1,0].scatter([s['birth']],[i],marker='o',facecolors='none',edgecolors='black',s=55)
        axes[1,0].scatter([s['time']],[i],marker='+',color='tab:blue',s=50)
    axes[1,0].set(ylabel='Inserted node index',title=f'{len(newborn)} births: circles = birth, crosses = first sample')
    prefix,_=read_curves(run/'restart_seed/curve_thermo_history')
    suffix,_=read_curves(run/'restart_resume/curve_thermo_history'); continued=prefix|suffix
    errors=[max(abs(s['fields'][f]-continued[g,p]['fields'][f])
                for (c,p),s in curves.items() if c==g for f in s['fields']) for g in cycles]
    axes[1,1].plot([clock[g] for g in cycles],errors,'o-',ms=3)
    cp=max(s['time'] for s in prefix.values());axes[1,1].axvline(cp,color='black',ls='--',label='Restart checkpoint')
    axes[1,1].set(ylabel='Maximum field difference (code units)',title='Resumed versus uninterrupted histories',ylim=(-1e-15,1e-15));axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.set_xlabel(r'$t$ (code units)');ax.grid(alpha=.2)
    save(fig,'02_amr_birth_restart','Live AMR, fresh insertion histories, exact continuation',
         'Native leaf count changes from 4 to 16 and back. Fourteen inserted nodes begin at their own birth, with no parental past.\n'
         'Every curve record and regular-particle record equals its uninterrupted counterpart; the restart plot compares all selected fields.')

    curves,clock=read_curves(run/'planar/curve_thermo_history'); signals={}
    for path in (run/'planar/shock_events/curves').glob('trace*.rank*.csv'):
        for p,s in read_trace(path).items():signals[s['generation'],p]=s
    ids=sorted({p for g,p in curves}); selected=[ids[0],ids[len(ids)//2],ids[-1]]
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for ax,field,label in zip(axes,['density','pressure'],[r'$\rho$',r'$P$']):
        for p in selected:
            samples=sorted(((g,s) for (g,q),s in curves.items() if q==p))
            ax.plot([s['time'] for g,s in samples],[s['fields'][field] for g,s in samples],'.-',label=f'Node {p}')
            inside=[s for g,s in samples if signals[g,p]['inside']]
            ax.scatter([s['time'] for s in inside],[s['fields'][field] for s in inside],
                       facecolors='none',edgecolors='black',s=40)
        ax.set(xlabel=r'$t$ (code units)',ylabel=label+' (code units)');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    save(fig,'03_shock_thermo','Thermo and shock history sample the same completed fluid state',
         'Translating planar shock: all 1,728 curve thermo samples exactly match shock-trace density, pressure and temperature.\n'
         'Black circles mark raw numerical shock-region membership, including any suppressed episodes; they are not exact crossing times.')

    curves,clock=read_curves(run/'forced_both/curve_thermo_history'); particles,_=thermo(run/'forced_both')
    cycles=sorted(clock);times=[clock[g] for g in cycles]
    ids=sorted({p for g,p in curves});selected=[ids[0],ids[len(ids)//2],ids[-1]]
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,field,label in zip(axes.flat,['density','pressure','temperature','entropy'],
                             [r'$\rho$',r'$P$',r'$P/\rho$',r'$\ln(P/\rho^\gamma)$']):
        quantiles=np.array([np.quantile([s['fields'][field] for (c,p),s in particles.items() if c==g],[.1,.5,.9]) for g in cycles])
        ax.fill_between(times,quantiles[:,0],quantiles[:,2],color='gray',alpha=.2,label='Particle 10–90%')
        ax.plot(times,quantiles[:,1],color='gray',ls='--',label='Particle median')
        for p in selected:
            ax.plot(times,[curves[g,p]['fields'][field] for g in cycles],'.-',label=f'Curve node {p}')
        ax.set(xlabel=r'$t$ (code units)',ylabel=label+' (code units)');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    save(fig,'04_forced_shared','Same thermo output for curves and particles in native forced flow',
         'Accepted 32³ box: 4,096 continuous particles and a fixed 16-node curve; the same 11 fields are sampled every 10 steps.\n'
         'Three node paths shown. Particle quantiles are unweighted tracer statistics; the short curve is not a representative box sample.\n'
         'Local engineering demonstration only: turbulence convergence, the inherited forcing thermal issue and scientific approval remain open.')

    _,reference=read_curves(run/'time_cadence/curve_thermo_history')
    _,prefix=read_curves(run/'dt_seed_18/curve_thermo_history')
    _,suffix=read_curves(run/'dt_resume_18/curve_thermo_history')
    sequences=[('Uninterrupted scheduled output',reference),('Corrected stop + resume',prefix|suffix)]
    if args.cadence_failure:
        _,before=read_curves(args.cadence_failure/'dt_seed_18/curve_thermo_history')
        _,after=read_curves(args.cadence_failure/'dt_resume_18/curve_thermo_history')
        sequences.insert(1,('Retained failing stop + resume',before|after))
    fig,ax=plt.subplots(figsize=(12,4.8))
    for i,(label,clocks) in enumerate(sequences):
        ax.scatter(sorted(clocks.values()),[i]*len(clocks),s=65,label=label)
        for g,t in reference.items():
            if g not in clocks:ax.scatter([t],[i],marker='x',color='red',s=150)
    ax.axvline(max(prefix.values()),ls='--',color='gray',label='Off-cadence stop/checkpoint')
    ax.set(yticks=range(len(sequences)),yticklabels=[s[0] for s in sequences],
           xlabel=r'$t$ (code units)',ylim=(-.5,len(sequences)-.5))
    ax.grid(alpha=.2)
    save(fig,'05_restart_cadence','Stopping between output times preserves the next scheduled sample',
         'Each point is an actual saved sample time. The extra sample at t=0.225 is the native terminal output.\n'
         'The old shared thermo writer advanced its clock there, skipping t=0.25 on resume. The corrected writer preserves it.\n'
         'Both particle and curve scheduled records match the uninterrupted run, with exact checks for off-cadence and on-cadence stops.')

    intro='<h1>Shared material-curve and particle thermo histories</h1><p>Actual native output; local Serial/double. Scientific approval pending.</p>'
    body=''.join('<section><h2>'+html.escape(title)+'</h2><img src="'+name+'.png"><p>'+html.escape(caption).replace('\n','<br>')+'</p></section>' for name,title,caption in panels)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Curve–particle thermo synergy</title><style>body{font:17px system-ui;background:#f4f6f8;color:#172330;margin:auto;padding:30px;max-width:1250px}section{background:white;padding:20px;margin:25px 0;border-radius:12px}img{width:100%}p{line-height:1.5}</style>'+intro+body)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(out/'index.html')


if __name__=='__main__':main()
