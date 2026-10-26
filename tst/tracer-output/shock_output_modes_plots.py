"""Actual saved-output plots for sampled shocks and binary encounter modes."""
import argparse
import json
from pathlib import Path
import re
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python')]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tracer_history import read_run
from plot_curves import select_curve,plot_curve
import shock_binary


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
 p.add_argument('--reference-runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 report=json.loads((a.runs/'report.json').read_text());assert len(report['checks'])==21 and all(c['status']=='pass' for c in report['checks'])
 data={name:read_run(a.runs/name) for name in ('sampled','complete','coarse','amr_complete','forced','oracle_dt')}
 figures=[]
 def save(fig,name,title,caption):
  fig.savefig(out/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
  figures.append(dict(path=name+'.png',title=title,caption=caption))
 def counts(name):
  log=(a.runs/name/'run.log').read_text()
  return [int(re.search(pattern+r'=(\d+)',log)[1]) for pattern in ('particle detector cache calls this process','curve stencil batches this process')]
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');names=['no_shock_fields','sampled','complete']
 values=np.array([counts(name) for name in names])
 for i,ax in enumerate(axes):
  bars=ax.bar(['No shock fields','Sampled only','Complete + snapshots'],values[:,i],color=['gray','tab:blue','tab:orange'])
  ax.bar_label(bars);ax.set(ylabel='Actual calls in this process',title=['Particle full-mesh cache builds','Curve node-stencil batches'][i]);ax.tick_params(axis='x',labelrotation=12)
 save(fig,'detection-work','No per-step finding in snapshot-only mode',
      'Same accepted 107-timestep planar run. Snapshot-only detection occurs at exactly 23 saved clocks. Complete mode evaluates 108 clocks including initial state; curve 131 includes 108 history batches plus 23 saved-field batches. Calls measure work frequency, not wall-time or equivalent CPU cost across the two kernels.')
 fig,axes=plt.subplots(2,1,figsize=(12,7),layout='constrained');full=data['complete'];coarse=data['coarse']
 node=coarse['curve_events'][0]['node_id'];samples=sorted((o for (g,n),o in coarse['curve_nodes'].items() if n==node),key=lambda o:o['time'])
 axes[0].plot([o['time'] for o in samples],[o['fields']['shock_region'] for o in samples],'o',ms=9,label='Only saved CRV shock flags')
 e=next(e for e in coarse['curve_events'] if e['node_id']==node)
 axes[0].axvspan(e['first']['time'],e['last']['time'],alpha=.25,label='Recorded inside evidence interval')
 axes[0].set(ylabel='Saved shock-region flag [0/1]',title=f'Node {node}: sparse snapshots miss the entire passage',ylim=(-.1,1.2));axes[0].legend()
 for i,e in enumerate(coarse['curve_events']):axes[1].plot([e['first']['time'],e['last']['time']],[i,i],lw=4)
 axes[1].set(ylabel='Curve-node encounter index',title='Separate complete-history block retains all 16 passages')
 for ax in axes:ax.set_xlabel(r'$t$ (code units)');ax.grid(alpha=.2)
 save(fig,'between-snapshots','Sampled flags and complete histories are distinct products',
      'Only cycles 0 and 107 are saved in ordinary CRV; all curve flags there are zero. The explicit complete-history block still records 16 curve passages. Shaded/line intervals use first/last numerical inside evidence, not inferred exact physical crossing times.')
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for carrier,ax in enumerate(axes):
  folder=a.runs/'complete'/'shock_events';folder=folder/'curves' if carrier else folder
  packets=[shock_binary.read(p) for p in shock_binary.binary_files(folder,carrier)];packets=[p for p in packets if p['kind']==0]
  for packet in packets:
   ax.axvline(packet['generation'],color='gray',alpha=.35)
   ax.scatter([packet['generation']]*packet['count'],[e['emission_generation'] for e in packet['events']],s=40)
  ax.set(xlabel='Binary flush cycle',ylabel='Original completion cycle',title=['Particle batches:11-cycle cadence','Curve batches:17-cycle cadence'][carrier]);ax.grid(alpha=.2)
 save(fig,'completion-versus-flush','Writing later preserves the original completion clock',
      'Dots represent actual complete records in typed event-only THP/CRV batches; vertical lines include empty flushes. Completion generation remains the timestep when the encounter engine finalized the episode. Normal snapshot fields are not duplicated.')
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 packets=[shock_binary.read(p) for p in (a.runs/'pending-packets').glob('*')]
 axes[0].bar(['Particles','Curve nodes'],[next(p['count'] for p in packets if p['carrier']==c) for c in (0,1)])
 axes[0].set(ylabel='Unflushed completed records',title='Actual checkpoint at cycle 80, before flush100')
 old=read_run(a.reference_runs/'planar');d=read_run(a.runs/'pending_resume',shock_directory=str((a.runs/'pending-journal').resolve()))
 for carrier,key in enumerate(('particle_events','curve_events')):
  axes[1].bar(carrier-.18,len(old[key]),.36,color='gray',label='Uninterrupted' if carrier==0 else None)
  axes[1].bar(carrier+.18,len(d[key]),.36,color='tab:blue',label='Resumed pending buffer' if carrier==0 else None)
 axes[1].set(xticks=[0,1],xticklabels=['Particles','Curve nodes'],ylabel='Completed records',title='All identities and evidence match exactly');axes[1].legend(fontsize=8)
 save(fig,'pending-restart','Checkpointed pending records are not lost or retimed',
      'The paired cycle 80 checkpoint contains 15 completed curve records that have not reached their configured flush. After restart the complete histories recover exactly 15 particle and 16 curve events, retaining every original emission clock and evidence field. Active/zero-step and clean-prefix restart tests also pass.')
 metrics=[]
 for carrier in (0,1):
  oldfolder=a.reference_runs/'planar/shock_events';newfolder=a.runs/'complete/shock_events'
  if carrier:oldfolder/='curves';newfolder/='curves'
  oldfiles=list(oldfolder.glob('cycle*.csv'))+list(oldfolder.glob('suspended*.csv'))+list(oldfolder.glob('manifest.rank*'))
  newfiles=shock_binary.binary_files(newfolder,carrier)
  metrics.append(dict(carrier=['particles','curves'][carrier],old_files=len(oldfiles),new_files=len(newfiles),old_bytes=sum(p.stat().st_size for p in oldfiles),new_bytes=sum(p.stat().st_size for p in newfiles)))
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');x=np.arange(2)
 for ax,oldkey,newkey,label in [(axes[0],'old_files','new_files','Actual encounter files'),(axes[1],'old_bytes','new_bytes','Actual encounter bytes')]:
  ax.bar(x-.18,[m[oldkey] for m in metrics],.36,label='Previous every-step CSV')
  ax.bar(x+.18,[m[newkey] for m in metrics],.36,label='Binary at requested cadence')
  ax.set(xticks=x,xticklabels=['Particles','Curve nodes'],ylabel=label);ax.legend(fontsize=8)
 save(fig,'encounter-storage','Same complete evidence, fewer encounter files',
      'Actual planar files include completed batches, censored views and old manifests; optional traces excluded. Both runs preserve exactly the same complete/partial episodes. The comparison intentionally includes both format and flush-cadence changes; it is not a throughput or memory benchmark.')
 (out/'storage.json').write_text(json.dumps(metrics,indent=2)+'\n')
 for case,title,caption in [('sampled','Sampled-only planar CRV: geometry + thermo + flags','No complete encounter block is enabled here. Geometry and sampled thermo/shock properties remain synchronized in one CRV; detections only occur at saved clocks.'),('amr_complete','Live AMR with explicit complete tracking','Negative signed IDs, new-node birth times, sampled fields and complete binary histories match the previous accepted native output. New nodes inherit no prior encounters.'),('forced','Short forced-box coexistence','Accepted 32³ engineering fixture. Binary output recovers 6904 particle and 28 curve events exactly; this short run does not establish converged turbulence statistics.')]:
  name=case+'-combined.png';plot_curve(select_curve(data[case],field='pressure'),out/name)
  figures.append(dict(path=name,title=title,caption=caption))
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');d=data['oracle_dt'];node=next(iter(d['curve_nodes']))[1]
 rows=sorted((o for (g,n),o in d['curve_nodes'].items() if n==node),key=lambda o:o['time']);t=np.array([o['time'] for o in rows])
 expected=(4-np.maximum(abs(.2+.3*t),abs(-.1+.2*t)))**2/(5/3);actual=np.array([o['fields']['pressure'] for o in rows])
 axes[0].plot(t,expected,label='Independent prescribed oracle');axes[0].plot(t,actual,'o',label='Saved CRV');axes[0].set(ylabel='Pressure (code units)');axes[0].legend()
 axes[1].plot(t,actual-expected,'o-');axes[1].set(ylabel='Saved minus oracle pressure (code units)',ylim=(-2e-12,2e-12))
 for ax in axes:ax.set_xlabel(r'$t$ (code units)');ax.grid(alpha=.2)
 save(fig,'independent-thermo','Independent thermo oracle at requested output times',
      'Accepted accelerating-flow fixture, dt 0.05 output. Nine saved clocks produce nine curve detector batches. Maximum error over all selected thermo fields is below the unchanged 2e-12 CT budget; no fitted tolerance.')
 (out/'figures.json').write_text(json.dumps(figures,indent=2)+'\n');print(json.dumps(dict(figures=len(figures),output=str(out))))


if __name__=='__main__':main()
