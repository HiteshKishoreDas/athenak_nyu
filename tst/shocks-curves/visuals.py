"""Plot actual interaction fixtures using the production readers; no new simulations.

Usage: python visuals.py --core CORE --forced FORCED --remap-report REPORT --output NEW
The intentionally blocked coexisting-AMR gate is plotted as a failure.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'vis/python'), str(ROOT/'tst/shocks-mc')]
import acceptance as a
import curve_shock_history as curve
from bin_convert import read_binary
from common import thermo
from shock_history import read_history as mc_history, read_suspended as mc_suspended

BLUE, ORANGE, RED, GREEN = '#2166ac', '#e08214', '#b2182b', '#23734a'
plt.rcParams.update({'font.size':11, 'axes.spines.top':False,
                     'axes.spines.right':False, 'axes.titleweight':'bold',
                     'savefig.facecolor':'white'})


def hashfile(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def rows(directory):
    return a.trace_rows(Path(directory)/'shock_events/curves', curve)


def clocks(observed):
    result = {}
    for (g, _), o in observed.items():
        if g in result:assert result[g] == o['time']
        result[g] = o['time']
    return result


def membership(ax, observed, node_ids=None):
    ids = sorted({p for g,p in observed}) if node_ids is None else node_ids
    times = clocks(observed); generations = sorted(times)
    matrix = np.full((len(ids),len(generations)),np.nan)
    pi = {p:i for i,p in enumerate(ids)}; gi = {g:i for i,g in enumerate(generations)}
    for (g,p), o in observed.items():
        matrix[pi[p],gi[g]] = int(o['inside']) if o['valid'] else 2
    ax.pcolormesh([times[g] for g in generations],np.arange(1,len(ids)+1),matrix,
                  shading='nearest',cmap=ListedColormap(['#d9e3ed',ORANGE,RED]),vmin=0,vmax=2)
    ax.set(xlabel=r'$t$ [code units]',ylabel='Node order (IDs retained in data)')
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('core','forced','remap-report','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args(); core=args.core.resolve(); forced=args.forced.resolve()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((core/'report.json').read_text())
    force_report=json.loads((forced/'report.json').read_text())
    checks={c['id']:c for c in report['checks']}
    assert sum(c['status']=='pass' for c in checks.values())==19
    assert checks['CSH04_CSH07_coexisting_live_amr_births']['status']=='fail'
    assert force_report['status']=='pass'
    for directory,test_report in ((core,report),(forced,force_report)):
        for run in test_report['runs']:
            assert hashfile(directory/run['name']/'run.log')==run['log_sha256']
            assert hashfile(directory/run['name']/'input.athinput')==run['input_sha256']
            binary=Path(run['command'][0]);assert hashfile(binary)==run['executable_sha256']
    captions=[];metrics={'scientific_approval':'pending','fresh_run_checks_passed':19,
                        'fresh_run_checks_total':20,'figures':[]}
    def save(fig,name,title,caption):
        fig.suptitle(title,fontsize=16,x=.05,ha='left',y=.98)
        fig.text(.05,.025,caption,fontsize=10,va='bottom',wrap=True)
        fig.tight_layout(rect=(.035,.12,.995,.92),w_pad=2.5)
        fig.savefig(out/(name+'.png'),dpi=175);plt.close(fig)
        captions.append((name,title,caption));metrics['figures'].append(name+'.png')

    # Independently reconstructed fluid density front versus recorded node samples.
    planar=core/'planar'; observed=rows(planar); events=curve.read_history(planar/'shock_events/curves')
    ids=sorted({p for g,p in observed});times=clocks(observed)
    assert len(events)==len(ids)==16 and all(e['status']==0 for e in events)
    fig,ax=plt.subplots(1,3,figsize=(15,5.3))
    for p in ids:
        samples=sorted((o for (g,q),o in observed.items() if p==q),key=lambda o:o['time'])
        ax[0].plot([o['time'] for o in samples],[o['x'][0] for o in samples],color=BLUE,lw=.7,alpha=.5)
        inside=[o for o in samples if o['inside']]
        ax[0].scatter([o['time'] for o in inside],[o['x'][0] for o in inside],s=3,color=ORANGE)
    fronts=a.density_fronts(planar)
    ax[0].plot([times[g] for g,t,x in fronts],[x for g,t,x in fronts],color='black',lw=2,label='Native density midpoint')
    ax[0].legend(fontsize=9);ax[0].set(xlabel=r'$t$ [code units]',ylabel=r'$x_1$ [code units]',title='Recorded node paths and fluid front')
    membership(ax[1],observed);ax[1].set_title('Raw shock membership: all 16 nodes')
    event_by_id={e['node_id']:e for e in events}
    for i,p in enumerate(ids,1):
        e=event_by_id[p];ax[2].plot(i,e['representative']['mach_pressure'],'o',color=BLUE)
        ax[2].plot(i,e['representative']['mach_temperature'],'x',color=ORANGE)
    ax[2].axhline(2,color='black',ls='--',label=r'Input RH Mach $M=2$')
    ax[2].plot([],[],'o',color=BLUE,label='Pressure estimate');ax[2].plot([],[],'x',color=ORANGE,label='Temperature estimate')
    ax[2].legend(fontsize=9);ax[2].set(xlabel='Initial node order',ylabel='Representative Mach estimate',title='Strength evidence, no tuned tolerance')
    metrics['planar']={'complete_events':len(events),'nodes':len(ids),
                       'density_brackets_covered':checks['CSH02_CSH03_planar']['density_brackets_covered_by_event_envelope']}
    save(fig,'01_planar','Planar shock: fluid motion, node sampling and encounter evidence',
         'CSH02–03. Blue lines join stored node positions; orange samples lie inside the detector region. Pale blue means outside.\n'
         'Each node records one completed passage. All 16 independent density-front brackets lie within encounter envelopes; no exact crossing time is inferred.')

    # Birth provenance: before a node exists there must be no history.
    amr=core/'amr_curveonly_live_shock'; ar=rows(amr); br=a.trace_rows(amr/'shock_events/curves',curve,True)
    states=a.snapshots(amr);ordered=sorted(states);earliest={}
    for (g,p),o in list(ar.items())+list(br.items()):
        if p not in earliest or o['time']<earliest[p]['time']:earliest[p]=o
    newborn={p:o for p,o in earliest.items() if o['birth']>0}
    born_inside={p:o for p,o in newborn.items() if o['inside']}
    assert len(newborn)==39 and len(born_inside)==17
    assert all(o['time']==o['birth'] and o['sequence']==0 for o in newborn.values())
    assert all(o['flags']&1 for o in born_inside.values())
    fig,ax=plt.subplots(1,3,figsize=(15,5.3))
    ax[0].step([states[g]['time'] for g in ordered],[len(states[g]['amr']['leaves']) for g in ordered],where='post',color=BLUE)
    ax[0].set(xlabel=r'$t$ [code units]',ylabel='Native AMR leaf blocks',title='Actual refinement and coarsening',yticks=[2,8])
    by_birth=sorted(earliest,key=lambda p:(earliest[p]['birth'],p));membership(ax[1],ar,by_birth)
    ax[1].set_title('Blank = no sample; orange = inside')
    for inside,color,label,marker in [(False,BLUE,'22 born outside','o'),(True,ORANGE,'17 born inside','x')]:
        samples=[o for o in newborn.values() if o['inside']==inside]
        ax[2].scatter([o['birth'] for o in samples],[o['time']-o['birth'] for o in samples],s=55,color=color,label=label,marker=marker)
    ax[2].axhline(0,color='black',lw=.6);ax[2].legend(fontsize=9)
    ax[2].set(xlabel='Insertion time [code units]',ylabel=r'First observation $t-t_{birth}$',title='Fresh at birth: zero invented past',ylim=(-.001,.001))
    metrics['amr_births']={'inside':17,'outside':22,'maximum_birth_delay':0.0,
                          'all_sequences_start_zero':True,'inside_left_censored':True}
    save(fig,'02_amr_births','Curve-only AMR: inserted nodes start independent histories',
         'CSH04/07. 39 positive-time births: 17 inside and 22 outside. All first observations occur at birth with sequence zero.\n'
         'Inside births carry unknown-entry censoring. Coincident zero-delay points overlap. This tests storage/timing; AMR detector calibration remains pending.')

    # Native 3D forced box projected onto x-y; never invent MC trajectories.
    native=forced/'forced_on';fr=rows(native); ft=clocks(fr); fe=curve.read_history(native/'shock_events/curves')
    me=mc_history(native/'shock_events');ms=mc_suspended(native/'shock_events')
    sp=a.snapshots(native);final=sp[max(sp)]
    binary=sorted((native/'bin').glob('*.bin'))[-1];data=read_binary(binary)
    assert data['n_mbs']==1 and data['Nx1']==data['Nx2']==data['Nx3']==32
    density=np.asarray(data['mb_data']['dens'])[0]
    column=density.mean(axis=0)*(data['x3max']-data['x3min'])
    history,ht=thermo(native);last_generation=max(ht)
    particles=[o['x'] for (g,p),o in history.items() if g==last_generation]
    positions=np.asarray(particles);assert len(positions)==4096
    fig=plt.figure(figsize=(14,6.5));grid=fig.add_gridspec(2,2)
    ax=[fig.add_subplot(grid[:,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,1])]
    image=ax[0].imshow(np.log10(column),origin='lower',extent=[0,1,0,1],cmap='magma',aspect='equal')
    fig.colorbar(image,ax=ax[0],label=r'$\log_{10}\Sigma$ [code units]')
    ax[0].scatter(positions[:,0],positions[:,1],s=2,c='cyan',alpha=.2,label='4096 MC positions (projection)')
    for c in final['curves']:
        points=np.asarray([n['x'] for n in c['nodes']]);wrapped=np.mod(points,1)
        ax[0].scatter(wrapped[:,0],wrapped[:,1],s=28,edgecolor='white',c=GREEN,label='Curve nodes (projection)')
    ax[0].legend(fontsize=8,loc='upper left');ax[0].set(xlabel=r'$x_1$ [code units]',ylabel=r'$x_2$ [code units]',title=r'3D box projected along $x_3$ at $t=0.2$')
    for target,label,events_,color in [(ax[1],'16 curve nodes',fe,GREEN),
                                       (ax[2],'4096 MC particles',me,BLUE)]:
        counts=Counter(e['emission_generation'] for e in events_)
        generations=sorted(ft);totals=np.cumsum([counts[g] for g in generations])
        target.step([ft[g] for g in generations],totals,where='post',color=color)
        target.set(xlabel=r'$t$ [code units]',ylabel='Completed records',
                   title=f'{label}: {len(events_)} completed records')
        target.grid(alpha=.15)
    assert all(e['status']==0 for e in me)
    metrics['forcing']={'curve_events':len(fe),'mc_events':len(me),'mc_population':4096,
                        'curve_nodes':16,'mc_terminal_censored':len(ms),'short_engineering_run':True}
    save(fig,'03_forced_coexistence','Native forced box: MC particles and material curves coexist',
         'CSH09. Actual 32³ native forcing run, unchanged physical input. Density is a line-of-sight integral; symbols are projected stored positions.\n'
         f'Right panels use separate count scales and populations. {len(ms)} MC terminal censored views excluded. Short-run engineering evidence; no converged turbulence claim.')
    del history,positions

    fig,ax=plt.subplots(1,2,figsize=(14,5.5));fi=sorted({p for g,p in fr})
    membership(ax[0],fr,fi);ax[0].set_title('Every-step flags, all 16 curve nodes')
    counts=Counter(e['node_id'] for e in fe)
    ax[1].bar(np.arange(1,17),[counts[p] for p in fi],color=GREEN)
    ax[1].set(xlabel='Node order (curve 7, low ID bits 1–16)',ylabel='Completed encounters per node',title='28 complete episodes; repeated passages retained',xticks=np.arange(1,17),yticks=range(5))
    assert len(fe)==28 and all(e['status']==0 for e in fe)
    save(fig,'04_forced_histories','Forced-box node histories: flags and completed passages',
         'Pale blue = outside; orange = unreduced shock-region membership. An episode can contain multiple inside samples and brief reentries.\n'
         'Bars count completed encounter records, not orange samples. Existing rearming rules distinguish passages; whole-curve statistics are derived afterward.')

    # Independent reanalysis of tracking passivity and actual stored checkpoint replay.
    fig,ax=plt.subplots(1,3,figsize=(15,5.3));fluid_errors=[];fluid_times=[]
    for path in sorted((native/'bin').glob('*.bin')):
        before=read_binary(path);after=read_binary(forced/'forced_off/bin'/path.name)
        maximum=0.0
        for key,value in before['mb_data'].items():
            x=np.asarray(value);y=np.asarray(after['mb_data'][key]);assert np.array_equal(x,y)
            scale=float(np.max(np.abs(x)));difference=float(np.max(np.abs(x-y)))
            maximum=max(maximum,difference/scale if scale else difference)
        fluid_errors.append(maximum);fluid_times.append(before['time'])
    ax[0].plot(fluid_times,fluid_errors,'o-',color=GREEN)
    ax[0].set(xlabel=r'$t$ [code units]',ylabel='Maximum normalized fluid difference',title='Tracking on/off: exact fluid arrays',ylim=(-1e-15,1e-15))
    expected=curve.read_history(planar/'shock_events/curves');comparison=[]
    for phase in ('outside','active','pending'):
        seed=core/('seed_'+phase)
        actual=curve.read_history(seed/'shock_events/curves')
        assert a.canonical(actual)==a.canonical(expected)
        assert a.canonical(curve.read_suspended(seed/'shock_events/curves'))==a.canonical(curve.read_suspended(planar/'shock_events/curves'))
        comparison.append(len(actual))
    ax[1].bar(np.arange(3)-.16,[len(expected)]*3,width=.32,label='Uninterrupted',color=BLUE)
    ax[1].bar(np.arange(3)+.16,comparison,width=.32,label='Restarted',color=GREEN)
    ax[1].set(xticks=range(3),xticklabels=['Outside','Active','Pending exit'],ylabel='Completed planar records',title='Restart from actual encounter phases',ylim=(0,20))
    ax[1].legend(fontsize=9)
    source=amr/'shock_events/curves';resumed=core/'amr_birth_seed/shock_events/curves'
    before=rows(amr);after=rows(core/'amr_birth_seed');assert before.keys()==after.keys()
    deltas={}
    for key,o in before.items():
        r=after[key];assert json.dumps(o,sort_keys=True)==json.dumps(r,sort_keys=True)
        g,p=key;deltas[g]=max(deltas.get(g,0),max(abs(x-y) for x,y in zip(o['x'],r['x'])))
    assert a.canonical(curve.read_history(source))==a.canonical(curve.read_history(resumed))
    gtime=clocks(before);ax[2].plot([gtime[g] for g in sorted(deltas)],[deltas[g] for g in sorted(deltas)],color=GREEN)
    ax[2].set(xlabel=r'$t$ [code units]',ylabel='Maximum node-position replay difference',title='AMR active-newborn replay: exact',ylim=(-1e-15,1e-15))
    metrics['reanalysis']={'fluid_on_off_exact':True,'phase_restarts_event_evidence_exact':True,
                           'amr_newborn_full_traces_exact':True,'zero_errors_real_not_clipped':True}
    save(fig,'05_passivity_restart','Interaction checks: tracking preserves evolution and replay preserves history',
         'Left and right values are exactly zero; the small symmetric axis range only makes zero visible. No positive error is clipped or replaced.\n'
         'Middle bars summarize counts, while assertions also compare full event evidence and censored records. AMR replay compares every node trace, not only final positions.')

    # Physical MC regrid failure, displayed beside every individual acceptance gate.
    remap=json.loads(args.remap_report.read_text());record=next(r for r in remap['runs'] if r['mode']=='occupancy')
    occ=record['occupancy'];assert occ['children']==[0,0,0,4096]
    for name,digest in remap['source_sha256'].items():assert hashfile(ROOT/name)==digest
    fig,ax=plt.subplots(1,2,figsize=(15,9),gridspec_kw={'width_ratios':[1,1.7]})
    x=np.arange(4);ax[0].bar(x-.18,occ['children'],.36,color=RED,label='Observed native remap')
    ax[0].bar(x+.18,[occ['count']*occ['uniform_child_mass_fraction']]*4,.36,color=BLUE,label='Equal-mass expectation')
    ax[0].set(xticks=x,xticklabels=['Child 0','Child 1','Child 2','Child 3'],ylabel='MC tracer count',title='Inherited MC AMR allocation fails')
    ax[0].legend(fontsize=9);ax[0].text(.04,.8,'4096 tracers all select one child.\nHistory-state preservation does not\nmake this a mass-following remap.',transform=ax[0].transAxes,fontsize=10)
    labels=['Planar + MC','Planar curve-only','Initial inside censoring','Exact passivity','Restart phases','Relative directory restart','Sidecar / prefix integrity','Bounded overflow rejection','Old checkpoint rejection','Consumer disable rejection','Native hydro contract','Finite forcing inputs','Smooth zero','Smooth translation','Smooth shear','Smooth refine/coarsen','Refinement callback contract','Live AMR births (curve-only)','AMR newborn restart','Live AMR + MC + curves']
    checks_list=report['checks'];assert len(labels)==len(checks_list)==20
    values=[c['status']=='pass' for c in checks_list]
    ax[1].barh(np.arange(20),[1]*20,color=[GREEN if v else RED for v in values],height=.72)
    ax[1].set(yticks=range(20),yticklabels=labels,xticks=[],xlim=(0,1.15),title='Fresh core suite: 19 pass; 1 blocked')
    ax[1].invert_yaxis()
    for i,v in enumerate(values):ax[1].text(.96,i,'PASS' if v else 'BLOCKED',ha='right',va='center',color='white',weight='bold',fontsize=9)
    ax[1].spines[['left','bottom']].set_visible(False)
    metrics['mc_amr']={'observed_children':occ['children'],'expected_child_count':1024,
                       'qualification':'failed; startup guard retained'}
    save(fig,'06_amr_blocker','Full integration is incomplete: the combined AMR gate stays red',
         'Left: focused native remap harness with explicit constructor/count stand-ins, not a full fluid transaction. Source hashes are checked against this branch.\n'
         'Right: fresh native acceptance results. MC+curve AMR deliberately rejects at startup. RUSTY/multi-rank/GPU execution and human scientific approval remain pending.')

    inputs=[core/'report.json',forced/'report.json',args.remap_report,Path(__file__)]
    metrics['input_sha256']={str(p):hashfile(p) for p in inputs}
    metrics['figure_sha256']={p.name:hashfile(p) for p in out.glob('*.png')}
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    text=['# Curve shock interaction plots — 2026-10-19','',
          'Native test runs repeated with unchanged inputs. Full integration is incomplete: 19/20 core gates pass; forcing suite passes.',
          'Scientific approval and RUSTY/multi-rank/GPU validation remain pending.','']
    for name,title,caption in captions:text.extend(['## '+title,'',caption,'','!['+title+']('+name+'.png)','',''])
    (out/'README.md').write_text('\n'.join(text))
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
