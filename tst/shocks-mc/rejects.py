"""Explicit safety/config rejection fixtures, retaining mutated checkpoint copies."""
import argparse
import ctypes
import json
from pathlib import Path
import shutil

from common import ROOT, provenance, run, sha, uniform_changes


class SidecarPrefix(ctypes.Structure):
    _fields_=[('magic',ctypes.c_char*8),('version',ctypes.c_uint32),('endian',ctypes.c_uint32),
              ('real_bytes',ctypes.c_uint32),('state_bytes',ctypes.c_uint32),('int_bytes',ctypes.c_uint32),
              ('rank',ctypes.c_int),('ranks',ctypes.c_int),('gids',ctypes.c_int),('nmb',ctypes.c_int),
              ('nrdata',ctypes.c_int),('nidata',ctypes.c_int),('generation',ctypes.c_uint64),
              ('count',ctypes.c_uint64),('contract_bytes',ctypes.c_uint64)]


def fnv(data):
    digest=1469598103934665603
    for byte in data:digest=((digest^byte)*1099511628211)&((1<<64)-1)
    return digest


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    report=args.report.resolve();report.parent.mkdir(parents=True,exist_ok=True)
    result=dict(provenance(args.executable),runs=[])
    def execute(label,changes=None,restart=None,overrides=(),reason=None):
        item=run(args.executable,work/label,changes,restart,overrides)
        item['label']=label;item['expected_diagnostic']=reason
        log=(work/label/'run.log').read_text()
        item['passed']=item['exit']==0 if reason is None else isinstance(item['exit'],int) and item['exit']!=0 and reason in log
        result['runs'].append(item);report.write_text(json.dumps(result,indent=2)+'\n')
        target=report.parent/'reject_runs'/label;target.mkdir(parents=True,exist_ok=True)
        for name in ('input.athinput','run.log'):
            if (work/label/name).exists():shutil.copy2(work/label/name,target/name)
        (target/'manifest.json').write_text(json.dumps(item,indent=2)+'\n')
        return work/label
    execute('duplicate_seed_ids',{'tracer_seed_duplicate/id':1,
        'tracer_seed_duplicate/count_per_event':1},reason='unique tracer seed schedule IDs')
    execute('pernode_restart',{'output1/single_file_per_node':'true','output1/dt':1},reason='per-node')
    execute('shearing_box',{'shearing_box/omega0':1,'shearing_box/qshear':1.5},reason='requires uniform')
    execute('fofc',{'hydro/fofc':'true'},reason='requires uniform')
    execute('tag_overflow',{'particles/next_tracer_tag':2147483647,
        'tracer_seed_initial/count_per_event':2},reason='identity/schedule range exhausted')
    changes=uniform_changes(count=2);changes.update({'shock_tracking/enabled':'true',
        'output1/dt':1,'time/tlim':.02})
    seed=execute('bound_seed',changes)
    if result['runs'][-1]['passed']:
        cp=sorted((seed/'rst').glob('*.rst'))[-1]
        execute('changed_floor',restart=cp,
            overrides=('output1/dt=-1','hydro/dfloor=.1'),reason='contract changed')
        execute('changed_rng',restart=cp,
            overrides=('output1/dt=-1','particles/random_seed=29'),reason='contract changed')
        for label,action,reason in [('future_sidecar','future','incompatible shock sidecar'),
                                    ('truncated_sidecar','truncated','sidecar'),
                                    ('corrupt_fluid','fluid','checksum')]:
            artifact=work/(label+'_artifacts');shutil.copytree(seed,artifact)
            checkpoint=artifact/'rst'/cp.name;side=Path(str(checkpoint)+'.shock.rank0')
            if action=='future':
                raw=bytearray(side.read_bytes());raw[8:12]=(99).to_bytes(4,'little')
                raw[-8:]=fnv(raw[:-8]).to_bytes(8,'little');side.write_bytes(raw)
            elif action=='truncated':side.write_bytes(side.read_bytes()[:-1])
            else:
                raw=bytearray(checkpoint.read_bytes());raw[-1]^=1;checkpoint.write_bytes(raw)
            execute(label,restart=checkpoint,overrides=('output1/dt=-1',
                'shock_tracking/directory='+str(artifact/'shock_events')),reason=reason)
        copy=work/'bound_artifacts';shutil.copytree(seed,copy)
        cp=copy/'rst'/cp.name;sidecar=Path(str(cp)+'.shock.rank0')
        raw=bytearray(sidecar.read_bytes());offset=SidecarPrefix.contract_bytes.offset
        raw[offset:offset+8]=(1<<20).to_bytes(8,'little')
        raw[-8:]=fnv(raw[:-8]).to_bytes(8,'little');sidecar.write_bytes(raw)
        result['bound_contract_mutation']={'native_header_offset':offset,'contract_bytes':1<<20,
            'sidecar_sha256':sha(sidecar),'footer_recomputed':True,
            'purpose':'reject impossible contract length before string allocation despite valid whole-sidecar checksum'}
        execute('bound_contract',restart=cp,
            overrides=('output1/dt=-1','shock_tracking/directory='+str(copy/'shock_events')),
            reason='incompatible shock sidecar')
    result['status']='passed_explicit_reject_checks' if all(x['passed'] for x in result['runs']) else 'failed'
    report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'report':str(report)}))


if __name__=='__main__':main()
