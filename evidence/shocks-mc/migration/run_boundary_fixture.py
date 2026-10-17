"""Compile bounded actual-method fixtures; no full Athena build, no physics approval."""
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[3]
# Run from the AthenaK checkout. Kokkos libraries come from the baseline CPU build.
flags = (root / 'build-baseline/src/CMakeFiles/athena.dir/flags.make').read_text()
args = []
for key in ('CXX_DEFINES', 'CXX_INCLUDES', 'CXX_FLAGS'):
    args.extend(shlex.split(re.search(r'^' + key + r' = (.*)$', flags, re.M)[1]))
libs = [root / 'build-baseline/kokkos' / path for path in (
    'containers/src/libkokkoscontainers.a', 'core/src/libkokkoscore.a',
    'simd/src/libkokkossimd.a')]
results = {'scope': 'actual boundary methods, stand-in constructors/count/exit callback; '
                    'MPI singleton is not inter-rank migration qualification',
           'source_sha256': {}, 'builds': [], 'cases': []}
for path in ('src/bvals/bvals.hpp', 'src/bvals/bvals_part.cpp',
             'src/shocks/encounter.hpp', 'src/particles/particles.hpp',
             'evidence/shocks-mc/migration/boundary_fixture.cpp'):
    results['source_sha256'][path] = hashlib.sha256((root / path).read_bytes()).hexdigest()
with tempfile.TemporaryDirectory(prefix='shock-mc-boundary-') as temp:
    temp = Path(temp)
    for mpi in (False, True):
        exe = temp / ('boundary-mpi' if mpi else 'boundary-serial')
        extra = []
        if mpi:
            config = (root / 'build-baseline/config.hpp').read_text()
            config = config.replace('#define MPI_PARALLEL_ENABLED 0',
                                    '#define MPI_PARALLEL_ENABLED 1')
            (temp / 'config.hpp').write_text(config)
            extra = ['-I' + str(temp), '-include', 'mpi.h']
        cmd = [('/opt/homebrew/bin/mpicxx' if mpi else '/usr/bin/c++'),
               *extra, *args, str(Path(__file__).with_name('boundary_fixture.cpp')),
               str(root / 'src/bvals/bvals_part.cpp'), str(root / 'src/globals.cpp'),
               *(str(p) for p in libs), '-o', str(exe)]
        built = subprocess.run(cmd, capture_output=True, text=True)
        results['builds'].append({'mpi': mpi, 'command': cmd,
                                 'exit': built.returncode, 'output': built.stdout + built.stderr})
        if built.returncode:
            continue
        modes = ['periodic2d', 'periodic3d', 'outflow', 'image_max', 'far', 'extent', 'disabled', 'duplicate']
        if mpi:
            modes.extend(['mixed', 'byte_limit'])
        for mode in modes:
            run = subprocess.run([str(exe), mode], capture_output=True, text=True, timeout=30)
            output = run.stdout + run.stderr
            expected = {'image_max': 'image counter limit',
                        'far': 'beyond adjacent MeshBlock', 'extent': 'extent mismatch',
                        'duplicate': 'duplicate/invalid shock particle removal',
                        'byte_limit': 'MPI byte count exceeds INT_MAX'}
            passed = (run.returncode != 0 and expected[mode] in output) if mode in expected \
                else run.returncode == 0 and 'PASS ' in output
            results['cases'].append({'mpi': mpi, 'mode': mode, 'exit': run.returncode,
                                     'output': output, 'expected_result_pass': passed})
results['all_passed'] = len(results['builds']) == 2 and \
    all(row['exit'] == 0 for row in results['builds']) and \
    len(results['cases']) == 18 and all(row['expected_result_pass'] for row in results['cases'])
Path(__file__).with_name('boundary_results.json').write_text(json.dumps(results, indent=2) + '\n')
for row in results['builds']:
    print('compile', 'MPI' if row['mpi'] else 'Serial', row['exit'])
for row in results['cases']:
    print('MPI' if row['mpi'] else 'Serial', row['mode'], row['expected_result_pass'],
          row['exit'], row['output'].splitlines()[0] if row['output'] else '')
raise SystemExit(0 if results['all_passed'] else 1)
