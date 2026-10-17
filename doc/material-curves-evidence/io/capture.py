"""Capture target-only I/O evidence without changing repository sources."""
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import runpy
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path('/Users/hitesh/.codex/worktrees/multiphase-material-curves/athenak_nyu')
BUILD = Path('/private/tmp/athenak-multiphase-curves-debug')
OUT = Path(__file__).resolve().parent
LIB = BUILD/'kokkos/core/src/libkokkoscore.a'
LAUNCHER = '/opt/homebrew/bin/prterun'
RUNTIME_COMMIT = '313977256c95e0aee17433afed27379f0b5579bc'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            h.update(block)
    return h.hexdigest()

def git(path, *args):
    return subprocess.check_output(['git', '-C', str(path), *args], text=True).strip()

source_files = ['src/material_curves/material_curves_io.cpp',
                'src/material_curves/material_curves.hpp', 'src/outputs/restart.cpp',
                'src/athena.hpp', 'src/globals.hpp', 'src/mesh/mesh.hpp',
                'src/mesh/mesh_refinement.hpp', 'src/mesh/meshblock_pack.hpp',
                'src/mesh/logical_location.hpp', 'vis/python/material_curves.py',
                'tst/scripts/material_curves_io/test_reader.py',
                'tst/scripts/material_curves_io/test_cpp.py',
                'tst/scripts/material_curves_io/roundtrip.cpp']
source_files = [s for s in source_files if (ROOT/s).is_file()]
for name in source_files:
    target = OUT/'source'/name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ROOT/name, target)
configured_files = ['CMakeCache.txt', 'compile_commands.json', 'config.hpp',
                    'git-state.txt', 'kokkos/KokkosCore_config.h',
                    'src/CMakeFiles/athena.dir/flags.make',
                    'src/CMakeFiles/athena.dir/link.txt']
for name in configured_files:
    if (BUILD/name).is_file():
        target = OUT/'build'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(BUILD/name, target)
manifest = {'recorded_utc': datetime.now(timezone.utc).isoformat(),
            'runtime_commit': RUNTIME_COMMIT, 'source_head': git(ROOT, 'rev-parse', 'HEAD'),
            'source_status': git(ROOT, 'status', '--short'),
            'build_git_state': (BUILD/'git-state.txt').read_text().strip(),
            'build_git_state_note': 'Configure-time value; source HEAD and exact source hashes below record the tested port.',
            'source_root': str(ROOT), 'build': str(BUILD),
            'kokkos_commit': git(ROOT/'kokkos', 'rev-parse', 'HEAD'),
            'kokkos_status': git(ROOT/'kokkos', 'status', '--short'),
            'python': sys.version, 'python_executable': sys.executable,
            'launcher': LAUNCHER,
            'hashes': {str(ROOT/s): digest(ROOT/s) for s in source_files}}
for p in [LIB, BUILD/'src/athena', Path(LAUNCHER), BUILD/'config.hpp',
          BUILD/'compile_commands.json']:
    manifest['hashes'][str(p)] = digest(p)
manifest['commands'] = [
    [sys.executable, str(ROOT/'tst/scripts/material_curves_io/test_reader.py')],
    [sys.executable, str(ROOT/'tst/scripts/material_curves_io/test_cpp.py'),
     '--build', str(BUILD), '--kokkos-library', str(LIB), '--launcher', LAUNCHER]]
(OUT/'provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
(OUT/'commands.txt').write_text('cwd='+str(ROOT)+'\nPYTHONDONTWRITEBYTECODE=1\n'+
    '\n'.join(shlex.join(c) for c in manifest['commands'])+'\n')
(OUT/'restart-port.diff').write_text(git(ROOT, 'diff', '526a187f', '--',
    'src/outputs/restart.cpp')+'\n')
original_run = subprocess.run
calls = OUT/'native-subprocesses.jsonl'
calls.write_text('')

def recording_run(command, *args, **kwargs):
    start = time.monotonic()
    try:
        result = original_run(command, *args, **kwargs)
    except BaseException as error:
        with calls.open('a') as f:
            f.write(json.dumps({'command': command, 'cwd': str(kwargs.get('cwd', ROOT)),
                'elapsed_seconds': time.monotonic()-start, 'exception': repr(error)})+'\n')
        raise
    with calls.open('a') as f:
        f.write(json.dumps({'command': command, 'cwd': str(kwargs.get('cwd', ROOT)),
            'elapsed_seconds': time.monotonic()-start, 'returncode': result.returncode,
            'stdout': result.stdout, 'stderr': result.stderr})+'\n')
    return result

results = []
for label, command in zip(['reader', 'native'], manifest['commands']):
    start = time.monotonic()
    with (OUT/(label+'.stdout')).open('w') as stdout, (OUT/(label+'.stderr')).open('w') as stderr:
        if label == 'reader':
            result = original_run(command, cwd=ROOT, stdout=stdout, stderr=stderr)
            code = result.returncode
        else:
            sys.argv = command[1:]
            subprocess.run = recording_run
            try:
                with redirect_stdout(stdout), redirect_stderr(stderr):
                    runpy.run_path(command[1], run_name='__main__')
                code = 0
            except BaseException:
                import traceback
                traceback.print_exc(file=stderr)
                code = 1
            finally:
                subprocess.run = original_run
    results.append({'check': label, 'command': command, 'returncode': code,
                    'elapsed_seconds': time.monotonic()-start,
                    'assertions': 'passed' if code == 0 else 'failed',
                    'stdout': label+'.stdout', 'stderr': label+'.stderr'})
(OUT/'results.json').write_text(json.dumps(results, indent=2)+'\n')
print(json.dumps(results, indent=2))
if any(r['returncode'] for r in results):
    raise SystemExit(1)
