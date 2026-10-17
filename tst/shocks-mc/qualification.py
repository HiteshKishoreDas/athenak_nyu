"""Held-out MC planar recovery evidence, independent forward RH state construction."""
import argparse
import json
import math
from pathlib import Path
import shutil
import time

from common import ROOT, provenance, run, sha
from integration import recovery


def rh_states(mach=2,gamma=1.4,speed=.5,reflected=False):
    compression=(gamma+1)*mach*mach/((gamma-1)*mach*mach+2)
    pressure=(2*gamma*mach*mach-(gamma-1))/(gamma+1)
    relative=mach*gamma**.5*(1 if reflected else -1)
    upstream={'d':1,'p':1,'u':speed+relative}
    downstream={'d':compression,'p':pressure,'u':speed+relative/compression}
    left,right=(upstream,downstream) if reflected else (downstream,upstream)
    result={'hydro/gamma':gamma,'problem/xshock':.7 if reflected else .3}
    for side,state in [('l',left),('r',right)]:
        result.update({'problem/'+key+side:value for key,value in state.items()})
    # Independent mass/momentum/energy conservation in the front frame.
    u1=upstream['u']-speed;u2=downstream['u']-speed
    for a,b in [(u1,compression*u2),(1+u1*u1,pressure+compression*u2*u2),
                (gamma/(gamma-1)+u1*u1/2,gamma*pressure/((gamma-1)*compression)+u2*u2/2)]:
        assert abs(a-b)<=1e-12*max(abs(a),abs(b),1)
    return result


def fixtures():
    result=[]
    for seed in (11,29,47,83,127,211,307,401):
        result.append(('seed'+str(seed),dict(rh_states(),**{'particles/random_seed':seed})))
    for mach in (1.2,1.5,3,5):
        changes=rh_states(mach);changes['time/tlim']=.45 if mach==1.2 else .3 if mach==1.5 else .25
        result.append(('mach'+str(mach),changes))
    result.extend([('gamma53',rh_states(gamma=5/3)),('stationary',rh_states(speed=0)),
                   ('boost2',rh_states(speed=2.5)),
                   ('reflected',dict(rh_states(speed=-.5,reflected=True),
                       **{'tracer_seed_initial/slab_min':.25,'tracer_seed_initial/slab_max':.55}))])
    for nx in (64,256):result.append(('nx'+str(nx),dict(rh_states(),**{'mesh/nx1':nx,'meshblock/nx1':nx})))
    for cfl in (.15,.45):result.append(('cfl'+str(cfl),dict(rh_states(),**{'time/cfl_number':cfl})))
    result.append(('planar3d',dict(rh_states(),**{'mesh/nx3':4,'meshblock/nx3':4,
        'mesh/ix3_bc':'periodic','mesh/ox3_bc':'periodic'})))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if not args.prepare_only and not args.executable:parser.error('--executable required to run')
    work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    report=args.report.resolve();report.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();results={'scope':'MC full-region recovery measured; no exact physical path/crossing-time estimator',
        'physical_acceptance':'owner stochastic recovery/strength/timing thresholds pending',
        'experimental_defaults':{'width':4,'min_mach':1.1,'integrator':'rk2'},'runs':[],'measurements':[]}
    if args.executable:results.update(provenance(args.executable))
    def save():
        results['elapsed_seconds']=time.monotonic()-start
        report.write_text(json.dumps(results,indent=2,allow_nan=True)+'\n')
    for label,changes in fixtures():
        if time.monotonic()-start>300:raise TimeoutError('qualification300s budget exhausted')
        changes['shock_tracking/directory']=str(work/label/'shock_events')
        if args.prepare_only:
            from common import edit_input
            directory=work/label;directory.mkdir()
            inp=directory/'input.athinput'
            inp.write_text(edit_input((ROOT/'inputs/tests/shock_mc.athinput').read_text(),changes))
            results['runs'].append({'label':label,'input_sha256':sha(inp),'status':'prepared_unrun'})
        else:
            item=run(args.executable,work/label,changes);item['label']=label;results['runs'].append(item)
            try:
                measurement=recovery(work/label) if item['exit']==0 else {'error':'process exit '+str(item['exit'])}
            except Exception as exc:measurement={'error':repr(exc)}
            measurement['label']=label
            expected_mach=float(label[4:]) if label.startswith('mach') else 2
            measurement['expected_input_mach']=expected_mach
            for field in ('pressure','temperature'):
                errors=[e[field+'_mach']-expected_mach for e in measurement.get('event_evidence',[])
                        if math.isfinite(e[field+'_mach'])]
                measurement[field+'_mach_error_range']=[min(errors),max(errors)] if errors else None
            results['measurements'].append(measurement)
            target=report.parent/'runs'/label;target.mkdir(parents=True,exist_ok=True)
            for name in ('input.athinput','run.log'):shutil.copy2(work/label/name,target/name)
            (target/'manifest.json').write_text(json.dumps(item,indent=2)+'\n')
        save()
    results['status']='prepared_unrun' if args.prepare_only else 'measured_scientific_acceptance_pending'
    save();print(json.dumps({'status':results['status'],'report':str(report)}))


if __name__=='__main__':main()
