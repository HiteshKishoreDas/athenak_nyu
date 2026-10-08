"""Independent analytic-front bins on actual MC cell-parent observations."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tst/shocks-mc'))
sys.path.insert(0,str(ROOT/'vis/python'))
from common import thermo
from shock_history import read_history
from qualification import fixtures


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--qualification',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    source=json.loads(args.qualification.read_text())
    cases=dict(fixtures());measurements=[]
    for item in source['runs']:
        label=item['label'];changes=cases[label]
        directory=Path(item['command'][item['command'].index('-i')+1]).parent
        observations,times=thermo(directory)
        reflected=label=='reflected'
        xfront=changes.get('problem/xshock',.3)
        # Forward RH constructor supplies u_upstream = speed +/- Mach sqrt(gamma).
        gamma=changes.get('hydro/gamma',1.4)
        if reflected:
            mach=((changes['problem/pr']*(gamma+1)+(gamma-1))/(2*gamma))**.5
            speed=changes['problem/ul']-mach*gamma**.5
        else:
            mach=((changes['problem/pl']*(gamma+1)+(gamma-1))/(2*gamma))**.5
            speed=changes['problem/ur']+mach*gamma**.5
        walks=defaultdict(list)
        for (generation,pid),o in observations.items():
            signed=(o['x'][0]-xfront-speed*o['time'])*(-1 if reflected else 1)
            walks[pid].append((generation,o['time'],signed))
        bins=[]
        for pid,walk in walks.items():
            walk.sort()
            for previous,current in zip(walk,walk[1:]):
                if previous[2]>0 and current[2]<=0:
                    bins.append({'particle_id':pid,'times':[previous[1],current[1]]})
        events=read_history(directory/'shock_events')
        matched=[]
        for interval in bins:
            candidates=[]
            for e in events:
                if e['particle_id']!=interval['particle_id']:continue
                lo=e['pre']['time'] if e['pre']['valid'] else e['first']['time']
                hi=e['post']['time'] if e['post']['valid'] else e['last']['time']
                if lo<=interval['times'][0] and hi>=interval['times'][1]:
                    candidates.append(e['sequence'])
            matched.append(dict(interval,event_sequences=candidates))
        measurements.append({'label':label,'parent_front_change_bins':len(bins),
            'committed_events':len(events),'bins_covered_by_event_envelope':
            sum(bool(x['event_sequences']) for x in matched),'matching':matched})
    report={'source_fingerprint':source.get('compiled_source_fingerprint'),
        'definition':'positive-to-nonpositive analytic-front distance at consecutive actual MC parent samples',
        'interpretation':'observable cell-parent transfer bins only; stochastic returns can create multiple bins; neither deterministic material paths nor exact crossing times',
        'tolerance':'inclusive recorded timestep endpoints, no fitted time tolerance',
        'physical_acceptance':'pending','measurements':measurements}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
