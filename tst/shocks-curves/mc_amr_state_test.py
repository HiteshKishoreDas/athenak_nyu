"""Bounded actual-native remap checks; constructor/count stand-ins are explicit."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kokkos-build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mpi', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    build = args.kokkos_build.resolve()
    flags = (build/'src/CMakeFiles/athena.dir/flags.make').read_text()
    compiler_args = []
    for key in ('CXX_DEFINES', 'CXX_INCLUDES', 'CXX_FLAGS'):
        compiler_args.extend(shlex.split(re.search(r'^'+key+r' = (.*)$',flags,re.M)[1]))
    overlay = output/'config'
    overlay.mkdir(exist_ok=True)
    config = (build/'config.hpp').read_text()
    config = re.sub(r'#define MPI_PARALLEL_ENABLED [01]',
                    '#define MPI_PARALLEL_ENABLED '+str(int(args.mpi)),config)
    (overlay/'config.hpp').write_text(config)
    cpp = Path(__file__).with_suffix('.cpp')
    command = [('/opt/homebrew/bin/mpicxx' if args.mpi else '/usr/bin/c++'),
               '-I'+str(overlay), '-I'+str(root/'src'), *compiler_args,
               '-ffunction-sections', '-fdata-sections', '-Wl,-dead_strip']
    if args.mpi:
        command.extend(['-include','mpi.h'])
    command.extend([str(cpp),str(root/'src/particles/particles_lagrangian_mc.cpp'),
                    str(root/'src/bvals/bvals_part.cpp'),str(root/'src/globals.cpp')])
    command.extend(str(build/'kokkos'/p) for p in (
        'containers/src/libkokkoscontainers.a','core/src/libkokkoscore.a',
        'simd/src/libkokkossimd.a'))
    exe = output/'mc_amr_state_test'
    command.extend(['-o',str(exe)])
    report = {'scope':'actual remap method with constructor/count stand-ins; '
                       'not full fluid AMR transaction or physical oracle',
              'mpi':args.mpi,'command':command,'source_sha256':{}}
    for p in (cpp,root/'src/particles/particles_lagrangian_mc.cpp',
              root/'src/bvals/bvals_part.cpp',root/'src/shocks/encounter.hpp'):
        report['source_sha256'][str(p.relative_to(root))] = hashlib.sha256(p.read_bytes()).hexdigest()
    result = subprocess.run(command,capture_output=True,text=True,timeout=120)
    report.update(compile_exit=result.returncode,compile_output=result.stdout+result.stderr)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('compile_exit',result.returncode,flush=True)
    if result.returncode:
        print(report['compile_output'])
        raise SystemExit(1)
    report['runs'] = []
    for mode in ('refine_coarsen','outside','extent','occupancy'):
        result = subprocess.run([str(exe),mode],capture_output=True,text=True,timeout=30)
        text = result.stdout+result.stderr
        passed = result.returncode == 0 if mode in ('refine_coarsen','occupancy') else \
            result.returncode != 0 and ('outside domain' if mode == 'outside' else 'extent') in text
        record = {'mode':mode,'exit':result.returncode,'output':text,'passed':passed}
        if mode == 'occupancy' and result.returncode == 0:
            line = next(line for line in text.splitlines() if line.startswith('AMR_OCCUPANCY '))
            record['occupancy'] = json.loads(line.removeprefix('AMR_OCCUPANCY '))
            record['scientific_mass_following_qualification'] = 'not passed: native remap collapses to one child'
        report['runs'].append(record)
        print(mode,result.returncode,passed,text[:200],flush=True)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert all(r['passed'] for r in report['runs'])


if __name__ == '__main__':
    main()
