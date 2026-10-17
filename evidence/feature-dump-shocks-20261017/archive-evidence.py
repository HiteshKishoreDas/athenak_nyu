from pathlib import Path
import hashlib,json,shutil,subprocess,tarfile
n=Path('/private/tmp/athenak-feature-dump-shocks-20261017')
w=Path('/private/tmp/athenak-feature-dump-shocks-tests')
a=Path('/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/feature-dump-merge-20261017');a.mkdir(exist_ok=False)
e=n/'evidence/feature-dump-shocks-20261017';e.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Preserve actual runs, failures, checkpoints and histories durably; original paths stay valid.
shutil.move(str(w),str(a/'raw'));w.symlink_to(a/'raw',target_is_directory=True)
for f in Path('/private/tmp').glob('athenak-merge-*.log'):shutil.copy2(f,e/f.name)
for f in Path('/private/tmp').glob('athenak-feature-dump-*-*.log'):shutil.copy2(f,e/f.name)
for folder in ('parity','integration','replay','qualification','forcing','curve-core',
               'curve-regression','curve-unsupported','old-contract','pre-fix-silent-curve-request',
               'interactions','interactions-fixed','interactions-final','interactions-reviewed'):
 src=w/folder
 if not src.exists():continue
 for p in src.rglob('*'):
  if p.is_file() and (p.suffix in ('.json','.log','.athinput') or p.name.endswith('.curves.final.txt')):
   dest=e/folder/p.relative_to(src);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
shutil.copy2(n/'tst/shocks-mc/feature_dump.py',e/'executed-interaction-check.py')
shutil.copy2('/private/tmp/archive-feature-dump-merge.py',e/'archive-evidence.py')
shutil.copy2('/private/tmp/athenak-merge-glue.diff',e/'merge-glue.diff')
for kind,path in (('built-in','/private/tmp/athenak-feature-dump-shocks-build'),('curves','/private/tmp/athenak-feature-dump-curves-build')):
 build=Path(path);runtime=a/'runtime'/kind;runtime.mkdir(parents=True)
 for filename in ('CMakeCache.txt','config.hpp','shock_version.hpp','src/athena','shock_detector_test','shock_encounter_test'):
  p=build/filename
  if p.is_file():
   target=runtime/p.name;shutil.copy2(p,target)
   if filename!='src/athena' and p.suffix in ('.txt','.hpp'):shutil.copy2(p,e/(kind+'-'+p.name))
with tarfile.open(a/'runtime/merged-source.tar.gz','w:gz') as archive:
 for filename in ('src','CMakeLists.txt','config.hpp.in','tst/shocks-mc/feature_dump.py'):
  archive.add(n/filename,arcname=filename)
raw_hashes={str(p.relative_to(w)):sha(p) for p in sorted(w.rglob('*')) if p.is_file()}
(e/'raw-sha256.json').write_text(json.dumps(raw_hashes,indent=2)+'\n')
reviews={
 'integration':{'finding':'Enabled shock request without Particles was silently ignored. Reproduced and fixed by seven-line physics-assembly fatal guard.','final_review':'No production blocker;11native interactions (8success/3expected rejects).','test_changes_applied':['Require512actual samples before comparing traces.','Restart consumed-input hash distinguished from unused fixture template.']},
 'numerical':{'result':'No timing/dataflow regression. Curve RK2 previous-stage primitives preserved; MC shock observations remain step-complete. Curve regrid/adaptation precedes curve outputs. Curves still reject particles/forcing.','scope':'Read-only source review, no runtime or scientific approval.'},
 'persistence':{'result':'Both independent loaders retained. Shock preflight precedes output mutation/open; checked metadata/close and durable shock sidecar retained; curve sidecar Save follows native close.','reviewed_main_sha256':sha(n/'src/main.cpp'),'reviewed_restart_sha256':sha(n/'src/outputs/restart.cpp'),'scope':'Read-only merge review, no combined curve/MC histories or hardware approval.'},
 'ponytail':'Minimal guard and existing fixtures; no new dependency, abstraction or physics algorithm.',
 'scientific_signoff':'pending'}
(e/'independent-reviews.json').write_text(json.dumps(reviews,indent=2)+'\n')
source=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu')
old=json.loads((source/'evidence/shocks-mc/forcing-enabled-20261013/regression/qualification-report.json').read_text())
new=json.loads((w/'qualification/report.json').read_text());assert old['measurements']==new['measurements']
record={
 'target_before':'cea1dafcf6830bfecba3e2aea5622b5c943a761c',
 'source_unchanged':subprocess.check_output(['git','rev-parse','codex/shock-multiphase'],cwd=n,text=True).strip(),
 'compiled_source_fingerprint':'48c53dc5c167604269d121412d4316215f09be8224e053d98af2eae9908a2b7b',
 'main_binary_sha256':sha(a/'runtime/built-in/athena'),
 'curve_binary_sha256':sha(a/'runtime/curves/athena'),
 'kokkos_commit':'08ceff92bcf3a828844480bc1e6137eb74028517',
 'source_sha256':{str(p.relative_to(n)):sha(p) for p in sorted((n/'src').rglob('*')) if p.is_file()},
 'new_check_sha256':sha(n/'tst/shocks-mc/feature_dump.py'),
 'raw_directory':str(a/'raw'),'raw_file_count':len(raw_hashes),
 'prior_21_heldout_measurement_dicts_exact':True,
 'checks':{'detector_parity':53,'mc_integration_processes':20,'mc_replay_processes':4,
           'heldout_measured_runs':21,'forcing_processes':30,'curve_core_records':31,
           'curve_regression_records':7,'curve_unsupported_processes':18,'interaction_processes':11,
           'curve_reader_tests':10,'curve_cpp_serial_io':'passed','encounter_state':'passed',
           'old_shock_contract_rejection':'passed_without_mutating_original_files'},
 'limitations':['Local Serial CPU/double only; no new MPI/GPU/AMR qualification.',
                'Curve nodes have no shock encounter state; curve+particle/forcing still rejected.',
                'Inherited native RK2 source-work physical oracle fails;21heldouts measured, scientific acceptance pending.',
                'Old shock checkpoints require matching old binary; new merged source contract rejects them.'],
 'scientific_signoff':'pending','release_status':'local engineering merge, no push'}
(e/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'packet':str(e),'raw_files':len(raw_hashes),'main_binary_sha256':record['main_binary_sha256'],'source_branch':record['source_unchanged']}))
