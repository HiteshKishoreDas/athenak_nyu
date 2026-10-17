"""Five-step same-partition replay; immutable original outputs are read-only."""
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu')
ORIGINAL=Path('/private/tmp/shock-highres128-long-20261014')
VIEW=ORIGINAL/'checkpoint-views/generation-1443'
OUT=Path(__file__).resolve().parent
BINARY=ROOT/'build-release/src/athena'
CHECKPOINT=ORIGINAL/'rst/shock_highres128_long.00001.rst'
START,END=1443,1448
sys.path.insert(0,str(ROOT/'tst/shocks-mc'))
sys.path.insert(0,str(ROOT/'vis/python'))
from common import provenance,sha
from shock_history import manifests,read_chunk,read_suspended,parse_observation
report={'status':'in_progress','scope':'Five-step same-partition Serial CPU/double replay of the existing128^3 native-forced run; no HPC/scientific qualification.',
        'source_run':str(ORIGINAL),'checkpoint_generation':START,'end_generation':END,
        'command':[], 'checks':[], 'scientific_signoff':'pending'}

def save():
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

try:
    report.update(provenance(BINARY))
    recovery=json.loads((VIEW/'recovery.json').read_text())
    shutil.copy2(VIEW/'recovery.json',OUT/'quarter-recovery.json')
    shutil.copy2(ORIGINAL/'input.athinput',OUT/'source-input.athinput')
    report['source_input_sha256']=sha(ORIGINAL/'input.athinput')
    report['checkpoint']=str(CHECKPOINT)
    report['checkpoint_sha256']=sha(CHECKPOINT)
    sidecar=Path(str(CHECKPOINT)+'.shock.rank0')
    report['sidecar_sha256']=sha(sidecar)
    for name in ('executable_sha256','compiled_source_fingerprint','checkpoint_sha256','sidecar_sha256'):
        assert report[name]==recovery[name],name+' changed since quarter recovery'
    cmake=(ROOT/'CMakeLists.txt').read_text()
    source_files=re.search(r'set\(shock_source_files (.*?)\)\n',cmake,re.S)[1].split()
    report['source_sha256']={name:sha(ROOT/name) for name in source_files}
    computed=hashlib.sha256(''.join(name+':'+report['source_sha256'][name]+';' for name in source_files).encode()).hexdigest()
    assert computed==report['compiled_source_fingerprint'],'current source differs from compiled fingerprint'
    report['helper_sha256']={name:sha(ROOT/name) for name in ('tst/shocks-mc/common.py','vis/python/shock_history.py')}
    report['script_sha256']=sha(__file__)
    clone=OUT/'shock_events';clone.mkdir()
    copied={}
    source=VIEW/'shock_events'
    for name in ['manifest.rank0']+[f'cycle{g}.rank0.csv' for g in range(START+1)]:
        src=source/name
        shutil.copy2(src,clone/name)
        copied[name]=sha(src)
        assert (clone/name).read_bytes()==src.read_bytes()
    (OUT/'copied-prefix-sha256.json').write_text(json.dumps(copied,indent=2)+'\n')
    report['copied_prefix_files']=len(copied)
    quarter=read_suspended(source,generation=START)
    shutil.copy2(source/f'suspended{START}.rank0.csv',OUT/f'quarter-suspended{START}.rank0.csv')
    report['quarter_censored_view']={'generation':START,'records':len(quarter),'identity_count':len({e['identity'] for e in quarter}),
                                    'separate_from_canonical_events':True}
    command=[str(BINARY),'-r',str(CHECKPOINT),'time/nlim='+str(END),
             'shock_tracking/directory='+str(clone),'shock_tracking/trace=true',
             'output4/dt=-1','job/basename=replay128']
    report['command']=command;save()
    started=time.monotonic()
    with (OUT/'run.log').open('w') as log:
        report['exit']=subprocess.run(command,cwd=OUT,stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
    report['elapsed_seconds']=time.monotonic()-started
    report['log_sha256']=sha(OUT/'run.log')
    save();assert report['exit']==0,'AthenaK failed'
    assert sha(CHECKPOINT)==report['checkpoint_sha256'] and sha(sidecar)==report['sidecar_sha256'],'original checkpoint pair changed'
    metadata=manifests(clone);assert list(metadata)==[0]
    lineage,contract=metadata[0]
    identities=set();event_count=0
    chunk_hashes={}
    for g in range(END+1):
        target=clone/f'cycle{g}.rank0.csv';original=ORIGINAL/'shock_events'/target.name
        assert target.read_bytes()==original.read_bytes(),target.name+' differs'
        for event in read_chunk(target,lineage,contract,g):
            assert event['identity'] not in identities,'duplicate canonical event identity'
            identities.add(event['identity']);event_count+=1
        if g>=START:chunk_hashes[target.name]=sha(target)
    report['event_chunks']={'all_prefix_and_continuation_byte_identical':True,'generations':[0,END],
                            'canonical_event_count':event_count,'unique_identities':len(identities),'continuation_sha256':chunk_hashes}
    traces=sorted(clone.glob('trace*.rank0.csv'),key=lambda p:int(re.fullmatch(r'trace(\d+)\.rank0\.csv',p.name)[1]))
    trace_generations=[int(re.fullmatch(r'trace(\d+)\.rank0\.csv',p.name)[1]) for p in traces]
    assert trace_generations==list(range(START+1,END+1)),trace_generations
    trace_hashes={};trace_rows={};terminal={}
    for path,g in zip(traces,trace_generations):
        original=ORIGINAL/'shock_events'/path.name
        assert path.read_bytes()==original.read_bytes(),path.name+' differs'
        seen=set()
        for row in csv.reader(path.read_text().splitlines()):
            assert len(row)==28,'trace schema'
            pid,sequence,phase,flags=map(int,row[:4]);assert pid not in seen,'duplicate particle observation'
            seen.add(pid);observation=parse_observation(row[4:]);assert observation['generation']==g
            if g==END and phase!=0:terminal[lineage,pid,sequence]=flags
        assert seen==set(range(8192)),'particle IDs lost or added'
        trace_rows[str(g)]=len(seen);trace_hashes[path.name]=sha(path)
    report['particle_traces']={'byte_identical':True,'generation_range':[START+1,END],'rows_per_generation':trace_rows,
                              'sha256':trace_hashes,'checkpoint_generation_reobserved':False,
                              'note':'ObserveShocks skips restored generation1443; comparisons cover exactly the five evolved generations1444..1448.'}
    censored=read_suspended(clone,generation=END)
    assert {e['identity'] for e in censored}==set(terminal),'censored identities do not equal terminal active/pending states'
    for event in censored:
        assert event['identity'] not in identities,'suspended view duplicates a completed/partial committed identity'
        assert event['status']==1 and event['flags']==(terminal[event['identity']]|8|64)
    report['terminal_censored_view']={'generation':END,'records':len(censored),'unique_identities':len(terminal),
        'exact_identity_flags_match_terminal_trace':True,'separate_from_canonical_events':True,
        'compared_to_uninterrupted_terminal_view':False,
        'note':'The uninterrupted run continues and has no terminal view at1448; validate resumed censoring independently against its byte-identical final trace, without counting views as completed events.'}
    report['checks']=['Original native checkpoint/sidecar/executable/source fingerprint match the isolated quarter recovery metadata.',
                      'All1449immutable canonical generations and rank manifest strictly parsed; events unique.',
                      'All canonical chunks0..1448equal uninterrupted bytes.',
                      'All five evolved trace generations equal uninterrupted bytes; exactly8192particle IDs each.',
                      'Terminal censored identities/status/flags agree with active/pending final state trace and remain separate.']
    report['status']='passed_same_partition_five_step_replay'
except Exception as exc:
    report['status']='failed';report['error']=repr(exc);save();raise
save()
print(json.dumps({'status':report['status'],'elapsed_seconds':report['elapsed_seconds'],
                  'events':report['event_chunks']['canonical_event_count'],'trace_rows':report['particle_traces']['rows_per_generation'],
                  'censored':report['terminal_censored_view']['records']}))
