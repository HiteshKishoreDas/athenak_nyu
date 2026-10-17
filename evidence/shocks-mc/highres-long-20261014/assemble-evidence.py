"""Freeze completed run evidence, keeping large raw histories outside Git."""
from pathlib import Path
import csv,hashlib,json,shutil,subprocess,sys,zipfile
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');run=Path('/private/tmp/shock-highres128-long-20261014');analysis=Path('/private/tmp/shock-highres128-analysis-20261014');trace=Path('/private/tmp/shock-highres128-trace-validation-20261014')
artifact=Path('/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/highres128-long-20261014');packet=Path('/private/tmp/shock-highres128-evidence-20261014')
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
r=json.loads((analysis/'report.json').read_text());assert r['status']=='completed_scientific_review_pending' and r['end_time']==1
execution=json.loads((run/'execution.json').read_text());assert execution['status']=='completed' and execution['exit']==0
validation=json.loads((trace/'progress.json').read_text());assert validation['status']=='completed' and validation['last_time']==1 and validation['validated_observations']==r['counters']['sampled']
assert subprocess.run(['git','diff','--quiet','ad4225edea78d255b67df97fe8f593c0c268fa1e','--','src','CMakeLists.txt'],cwd=root).returncode==0
packet.mkdir(exist_ok=False);artifact.mkdir(parents=True,exist_ok=False)
runtime=artifact/'runtime';runtime.mkdir()
source_paths=subprocess.check_output(['git','ls-files','src','CMakeLists.txt','.gitmodules','cmake','tst/shocks','tst/shocks-mc','vis/python/bin_convert.py','vis/python/shock_history.py','docs/shocks-mc','inputs/tests/shock_mc_highres_long.athinput','inputs/tests/shock_mc_showcase.athinput','evidence/shocks/reference'],cwd=root,text=True).splitlines()
subprocess.run(['git','archive','--format=tar.gz','--output='+str(runtime/'athenak-source.tar.gz'),'HEAD',*source_paths],cwd=root,check=True)
assert subprocess.check_output(['git','status','--porcelain'],cwd=root/'kokkos',text=True)==''
subprocess.run(['git','archive','--format=tar.gz','--output='+str(runtime/'kokkos-source.tar.gz'),'HEAD'],cwd=root/'kokkos',check=True)
for source,name in [(root/'build-release/src/athena','athena'),(root/'build-release/shock_detector_test','shock_detector_test'),(root/'build-release/CMakeCache.txt','CMakeCache.txt')]:
 shutil.copy2(source,runtime/name);assert sha(source)==sha(runtime/name)
# Same-volume renames keep every full history and checkpoint without duplicating tens of GiB.
run.rename(artifact/'raw');run.symlink_to(artifact/'raw',target_is_directory=True)
analysis.rename(artifact/'plots');analysis.symlink_to(artifact/'plots',target_is_directory=True)
trace.rename(artifact/'trace-validation');trace.symlink_to(artifact/'trace-validation',target_is_directory=True)
for path in (artifact/'plots').iterdir():
 if path.is_file() and (path.suffix in ('.png','.svg','.pdf','.npz','.json') or path.name in ('episodes.csv','evolution.csv','selected-traces.csv','frames.jsonl','baseline64-starts.csv','index.html')):
  shutil.copy2(path,packet/path.name)
for name in ('input.athinput','execution.json','run.log'):shutil.copy2(run/name,packet/name)
for name in ('progress.json','trace-sha256.json'):shutil.copy2(trace/name,packet/('trace-'+name))
files={
 'executed-launcher.py':'/private/tmp/start-highres-long.py',
 'executed-analysis.py':str(root/'tst/shocks-mc/highres_long_analysis.py'),
 'executed-pdf.py':str(root/'tst/shocks-mc/highres_long_pdf.py'),
 'executed-trace-validator.py':'/private/tmp/collect-highres-traces.py',
 'executed-checkpoint-recovery.py':'/private/tmp/recover-highres-checkpoint.py',
 'executed-checkpoint-monitor.py':'/private/tmp/highres-checkpoint-monitor.py',
 'assemble-evidence.py':__file__,
 'analysis.log':'/private/tmp/highres-analysis-final.log',
 'pdf-generation.log':'/private/tmp/highres-pdf-final.log',
 'trace-validation.log':'/private/tmp/shock-highres128-trace-validation-20261014.log',
 'previous-manifest.json':str(root/'evidence/shocks-mc/final-manifest.json'),
 'corrected-portable-input.athinput':str(root/'inputs/tests/shock_mc_highres_long.athinput')}
for name,source in files.items():shutil.copy2(source,packet/name)
for folder,source in [('width8-parity','/private/tmp/shock-highres128-width8-20261014'),('quarter-replay','/private/tmp/shock-highres128-quarter-replay-20261014'),('quarter-replay80','/private/tmp/shock-highres128-quarter-replay80-20261014'),('cadence-smoke','/private/tmp/shock-highres128-cadence-smoke-fixed-20261014'),('cadence-smoke-portable','/private/tmp/shock-highres128-cadence-smoke-portable-20261014'),('cadence-smoke-failed','/private/tmp/shock-highres128-cadence-smoke-20261014')]:
 target=packet/folder;target.mkdir()
 for path in Path(source).iterdir():
  if path.is_file() and path.suffix in ('.py','.json','.log','.athinput','.txt'):
   shutil.copy2(path,target/path.name)
 # Width8 tiny native CSV fixtures are part of the retained contract evidence.
 if folder=='width8-parity':
  for path in Path(source).glob('*.csv'):shutil.copy2(path,target/path.name)
shutil.copy2('/private/tmp/check-highres-width8.py',packet/'width8-parity/executed-check.py')
shutil.copy2('/private/tmp/highres-cadence-smoke.py',packet/'cadence-smoke/executed-check.py')
shutil.copy2('/private/tmp/highres-cadence-smoke-portable.py',packet/'cadence-smoke-portable/executed-check.py')
shutil.copy2('/private/tmp/highres-prefix-equivalence.py',packet/'cadence-smoke/executed-prefix-check.py')
views=packet/'checkpoint-recovery';views.mkdir()
for path in (run/'checkpoint-views').glob('*/recovery.json'):
 target=views/path.parent.name;target.mkdir();shutil.copy2(path,target/path.name);shutil.copy2(path.parent/'run.log',target/'run.log')
# Trace hashes were computed while validating every immutable original file.
rawhashes={'shock_events/'+name:digest for name,digest in json.loads((trace/'trace-sha256.json').read_text()).items()}
for path in sorted(run.rglob('*')):
 if path.is_file() and not path.relative_to(run).parts[0]=='checkpoint-views' and not path.name.startswith('trace'):
  rawhashes[str(path.relative_to(run))]=sha(path)
(packet/'raw-sha256.json').write_text(json.dumps(rawhashes,indent=2)+'\n')
sources=['src/shocks/detector.hpp','src/shocks/encounter.hpp','src/particles/shock_particles.cpp','src/srcterms/turb_driver.cpp','tst/shocks/detector_test.cpp','tst/shocks-mc/highres_long_analysis.py','tst/shocks-mc/highres_long_pdf.py','inputs/tests/shock_mc_highres_long.athinput']
source_record={'analysis_tools_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'runtime_start_commit':execution['source_commit'],'production_core_commit':'ad4225edea78d255b67df97fe8f593c0c268fa1e','production_src_and_cmake_unchanged':True,'compiled_source_fingerprint':execution['compiled_source_fingerprint'],'source_sha256':{name:sha(root/name) for name in sources},'scientific_signoff':'pending','large_raw_current_path':str(artifact/'raw'),'original_raw_alias':str(run),'current_plot_path':str(artifact/'plots'),'trace_validation_path':str(artifact/'trace-validation'),'preserved_runtime':{str(path):sha(path) for path in runtime.iterdir()},'kokkos_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root/'kokkos',text=True).strip(),'kokkos_clean':True}
(packet/'source-verification.json').write_text(json.dumps(source_record,indent=2)+'\n')
base=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/evidence/shocks-mc/supersonic-showcase-20261013')
# Independent compact CSV recomputation of final counts and common-window starts.
identities=set();counts={kind:0 for kind in ('complete','partial','open')};starts=0;last_sequence={}
with (packet/'episodes.csv').open() as f:
 for row in csv.DictReader(f):
  pid,seq=int(row['particle_id']),int(row['sequence']);key=(pid,seq);assert key not in identities and seq==last_sequence.get(pid,-1)+1
  identities.add(key);last_sequence[pid]=seq;counts[row['outcome']]+=1
  first,last=float(row['first_inside']),float(row['last_inside']);assert 0<=first<=last<=1
  if row['pre_time']:assert float(row['pre_time'])<=first
  if row['post_time']:assert last<=float(row['post_time'])<=1
  starts+=first<=.3
assert all(counts[k]==r['outcomes'].get(k,0) for k in counts)
assert starts==r['common_window_comparison']['highres128_starts']
with (base/'episodes.csv').open() as f:old=sum(x['case']=='baseline' and float(x['first_inside'])<=.3 for x in csv.DictReader(f))
assert old==r['common_window_comparison']['baseline64_starts']
(record:=packet/'independent-csv-check.json').write_text(json.dumps({'unique_identities':len(identities),'outcomes':counts,'common_t0p3_starts128':starts,'common_t0p3_starts64':old,'status':'passed'},indent=2)+'\n')
wall=execution['wall_seconds'];last=r['field_evolution'][-1]
readme=f'''# Larger longer forced-box evidence - 2026-10-14\n\nActual local Serial CPU/double run:128^3,8192MC particles, initial Mach5,\ncontinuous native edot50, gamma1.4/RK2/PLM/HLLC/CFL.30, t=1.0.\nCompared with64^3/t.3:8times cells,3.333times code-time span.\nWidth8/nghost8 keeps physical detector half-span1/16box (earlier width4/64).\nNo production source changed. Every-step raw shock flags and phases saved.\n\nRun exit0; generation{r['through_generation']}; wall{wall/60:.3f}minutes.\nComplete{counts['complete']:,}, partial{counts['partial']:,}, open{counts['open']:,};\ncomplete/particle{r['complete_per_particle']:.6f}. Native sampled\n{r['counters']['sampled']:,}, invalid0; all8192IDs every frame.\nFinal volume rmsMach{last['mach_volume']:.6f}, supersonic volume\n{100*last['supersonic_volume_fraction']:.3f}percent. Adiabatic heating changes\nsound speed; this is an evolving forced realization, not stationary Mach5.\n\nExact width8 NUM08 parity/NUM10 independent planarMach3 oracle retain1e-12\nstrength tolerance and exact masks. Strict bounded history parsing verifies\ncanonical identities, sequences, counters, brackets and censored records.\nFive-step native same-partition quarter restart matches40,960trace rows and\nall1449canonical chunks byte for byte; terminal records are checked separately.\nCheckpoint recovery evolves zero steps. Full8192thermo samples agree with\noriginal traces; native plotted z64 slices match actual live particles' masks,\nstrengths, normals and endpoint evidence at1e-12. Crop support applies only\nto the central plane, not whole-volume mask reconstruction.\n\nOriginal dcycle=-1 key presence suppressed intermediate volume output; this\nfrozen failing cadence configuration remains in input.athinput. Same-binary\npaired checkpoints supply actual quarter-time maps. Corrected portable input\nomits those keys; smoke produces0/1/2, with an exact original two-step sample\nprefix check. Initial CLI nlim attempt failed before evolution; successful\nvariant explicitly declares nlim. SystemPython lackedNumPy on two helper\nattempt; the installed scientific environment executes the actual checks.\n\nScientific approval is pending. The inherited RK2 source thermal oracle\nremains a failure; plateau/physical strength truth, stationary ensembles,\nphysical convergence tolerance, MPI/GPU/RUSTY/Frontier budgets remain open.\nNumerical elapsed spans retain entry/exit bounds; no exact crossing time or\ncontinuous MC trajectory is invented. All flags and open identities remain.\n\n## Inspect and reproduce\n\nThe PDF has explanations adjacent to each of7vector figure families. PNG/SVG/\nPDF plots, episodes.csv, evolution.csv, selected-traces.csv and frames.jsonl\nare included. Complete per-step all-particle histories and paired checkpoints\nare preserved at `{artifact/'raw'}` (outsideGit). Exact executable, detector fixture binary and build cache are retained in `{runtime}`. The original `/private/tmp`\npath is a symlink to this durable copy, keeping executed metadata valid.\n`raw-sha256.json` binds every primary raw file; checkpoint-recoveryJSON binds\nthe input pair for every recovered map. The packet intentionally omits large\nraw cubes/traces and native sliceCSV fixtures; these remain in the raw/plot\nartifact directories. Source-verificationJSON gives the current absolute paths.\n\nFresh run, in an empty directory:\n\n```sh\nATHENA={root}/build-release/src/athena\nINPUT={root}/inputs/tests/shock_mc_highres_long.athinput\n"$ATHENA" -i "$INPUT" shock_tracking/directory=/absolute/new/events\n```\n\nUse a new directory and a distinct lineage when combining separate realizations.\nAllow tens ofGiB for fulltrace=true. All native restart files must remain paired\nwith their .shock.rankN sidecars and canonical event prefix. Samebinary and\nqualified samepartition only; the feature restart contract rejects incompatible\nstate. Five-step replay_check.py records a concrete restart command.\n\nAnalysis (existingNumPy/Matplotlib scientific environment):\n\n```sh\n/Users/hitesh/hitesh/git/shockit/.venv/bin/python {root}/tst/shocks-mc/highres_long_analysis.py \\\n  --root {root} --run {artifact/'raw'} \\\n  --trace {artifact/'trace-validation'} --baseline {base} \\\n  --output /absolute/new/plots --pdf\n```\n\nThe frozen executed-trace-validator.py provides bounded fulltrace collection\nfor a fresh equivalent run; it takes --run and --output. Original input/config\ncommands, all logs/scripts, independent reviews/checks and frozen source hashes\nare retained. BundledReportLab/pypdf Python runs highres_long_pdf.py with the\nplot packet and --root{root}; PDF pages were parsed/rendered/visually checked.\nNew source/evidence commits are local, not pushed.\n'''
(packet/'README.md').write_text(readme)
shutil.copy2(packet/'README.md',artifact/'README.md')
# No self-hash; final archive checks its own members and contains the hash ledger.
hashes={str(p.relative_to(packet)):sha(p) for p in sorted(packet.rglob('*')) if p.is_file()}
(packet/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
with zipfile.ZipFile(packet/'figures-and-data.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for path in sorted(packet.rglob('*')):
  if path.is_file() and path.name!='figures-and-data.zip':z.write(path,str(path.relative_to(packet)))
with zipfile.ZipFile(packet/'figures-and-data.zip') as z:assert z.testzip() is None
shutil.copy2(packet/'figures-and-data.zip',artifact/'figures-and-data.zip')
print(json.dumps({'status':'evidence_completed','packet':str(packet),'artifact':str(artifact),'raw_files':len(rawhashes),'archive_sha256':sha(packet/'figures-and-data.zip'),'counts':counts}))
