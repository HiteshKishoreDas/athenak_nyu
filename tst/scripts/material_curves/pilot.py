"""Bounded S01/S02 local CPU resource pilot, not production scaling acceptance.

Use --executable, --artifacts and optional --source/--revision. Runs the same
32x32 prescribed-translation fixture for 128 fluid steps: curves off and on with
128/512/4096 fixed samples, three interleaved repeats. Native /usr/bin/time -l
captures complete-process wall time and peak RSS (macOS bytes). MPI singleton
startup, setup, host sampling/geometry and initial/final output are all included;
these totals do not isolate host compute from IO. No performance threshold is
imposed. Final full-grid float32 binary payload identity is compared to baseline;
this output-level identity does not establish bitwise double internal-state equality.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import re
import statistics
import shlex
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
COUNTS = (0, 128, 512, 4096)
REPEATS = 3
STEPS = 128


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def fluid_identity(path):
    """Version1.1 native output payload: mesh metadata and float32 conserved cells."""
    with path.open('rb') as stream:
        require(stream.readline().strip() == b'Athena binary output version=1.1',
                'unsupported native binary output version')
        for _ in range(7):
            stream.readline()
        offset_line = stream.readline().decode('ascii')
        match = re.fullmatch(r'\s*header offset=(\d+)\s*', offset_line)
        require(match is not None, 'bad native output header offset')
        offset = int(match[1])
        require(offset < path.stat().st_size, 'bad native parameter length')
        stream.seek(offset, 1)
        payload = stream.read()
    require(payload, 'empty fluid payload')
    return {'sha256': hashlib.sha256(payload).hexdigest(), 'bytes': len(payload),
            'scope': 'native float32 conserved-grid payload plus mesh metadata',
            'precision': 'float32 output; internal fluid double equality unverified'}


def timed_metrics(path):
    text = path.read_text()
    wall = re.search(r'([0-9.]+)\s+real\s+([0-9.]+)\s+user\s+([0-9.]+)\s+sys', text)
    rss = re.search(r'(\d+)\s+maximum resident set size', text)
    require(wall and rss, 'missing /usr/bin/time -l wall/RSS measurement')
    return {'wall_seconds_native_time': float(wall[1]), 'user_seconds': float(wall[2]),
            'system_seconds': float(wall[3]), 'peak_rss_bytes': int(rss[1])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--artifacts', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--revision')
    args = parser.parse_args()
    executable, artifacts = args.executable.resolve(), args.artifacts.resolve()
    require(platform.system() == 'Darwin', '/usr/bin/time -l pilot requires macOS')
    require(executable.is_file() and not artifacts.exists(),
            'missing executable or existing artifact directory (preserve old results)')
    build = executable.parent.parent
    cache = (build/'CMakeCache.txt').read_text()
    source_match = re.search(r'^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$', cache, re.MULTILINE)
    source = args.source.resolve() if args.source else Path(source_match[1])
    fixture = (source/'inputs/tests/material_curves/prescribed.athinput').read_text()
    config = (build/'config.hpp').read_text()
    require('"material_curve_test"' in config, 'requires prescribed validation fixture')
    revision = args.revision or subprocess.check_output(
        ['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    spec = importlib.util.spec_from_file_location('curve_reader',
                                                source/'vis/python/material_curves.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    artifacts.mkdir(parents=True)
    original_executable = executable
    executable = artifacts/"athena-frozen"
    shutil.copy2(original_executable, executable)
    metadata = {'status': 'running', 'revision': revision, 'executable': str(executable),
                'original_executable': str(original_executable),
                'executable_sha256': sha(executable), 'hardware': platform.platform(),
                'config': config, 'ranks': 1, 'backend': 'CPU Serial',
                'mpi_singleton_startup_included': '#define MPI_PARALLEL_ENABLED 1' in config,
                'compiler_flags': (build/'src/CMakeFiles/athena.dir/flags.make').read_text(),
                'controls': {'grid': [32, 32, 1], 'steps': STEPS, 'repeats': REPEATS,
                             'node_counts': list(COUNTS), 'adaptive': False,
                             'interior_outputs': False}, 'records': [],
                'limitations': ['Whole-process wall/RSS includes startup/setup and IO.',
                                'IO bytes reported; compute-only cost not isolated.',
                                'Small CPU singleton pilot, not production scaling.',
                                'No speed, memory or scaling acceptance threshold.',
                                'Fluid identity applies to native float32 output only.']}
    command_entry = json.loads((build/'compile_commands.json').read_text())[0]
    size_command = shlex.split(command_entry['command'])
    output_index = size_command.index('-o')
    del size_command[output_index:output_index+2]
    size_command.remove('-c')
    size_command.pop()
    size_source = artifacts/'struct_sizes.cpp'
    size_source.write_text('#include <iostream>\n'
        '#include "material_curves/material_curves.hpp"\n'
        'int main() { std::cout << sizeof(material_curves::Node) << " " '
        '<< sizeof(material_curves::Curve) << "\\n"; }\n')
    size_binary = artifacts/'struct_sizes'
    size_command += [str(size_source), '-o', str(size_binary)]
    subprocess.run(size_command, cwd=command_entry['directory'], check=True,
                   capture_output=True, text=True)
    sizes = subprocess.check_output([str(size_binary)], text=True).split()
    metadata['struct_sizes_bytes'] = dict(zip(('Node', 'Curve'), map(int, sizes)))
    metadata['push_workspace_model'] = {
        'device_bytes_per_node': 72, 'mirror_bytes_per_node_if_distinct': 72,
        'total_bytes_per_node_if_distinct': 144,
        'assumptions': 'Three Nx3 double device vectors plus three distinct host mirrors; '
                       'CPU create_mirror_view may alias. This model is not measured RSS.'}
    write(artifacts/'environment.json', metadata)
    records = metadata['records']
    try:
        for repeat in range(REPEATS):
            for count in COUNTS:
                directory = artifacts/f'r{repeat}_nodes{count}'
                directory.mkdir()
                text = fixture
                overrides = {'job/basename': 'pilot', 'time/tlim': 1000,
                             'time/nlim': STEPS, 'mesh/nx1': 32, 'mesh/nx2': 32,
                             'material_curves/enabled': 'true' if count else 'false',
                             'material_curves/adaptive': 'false',
                             'material_curves/output_dcycle': 0,
                             'curve1/count': max(count, 2)}
                for key, value in overrides.items():
                    block, name = key.split('/')
                    text += f'\n<{block}>\n{name} = {value}\n'
                text += ('\n<output1>\nfile_type = bin\nvariable = hydro_u\n'
                         'dt = 1000\n')
                input_path = directory/'input.athinput'
                input_path.write_text(text)
                command = ['/usr/bin/time', '-l', '-o', str(directory/'time.txt'),
                           str(executable), '-i', str(input_path)]
                require(sha(executable) == metadata['executable_sha256'],
                        'executable changed during pilot')
                started = time.monotonic()
                process = subprocess.run(command, cwd=directory, capture_output=True,
                                         text=True, timeout=120)
                (directory/'stdout.txt').write_text(process.stdout)
                (directory/'stderr.txt').write_text(process.stderr)
                require(process.returncode == 0, 'pilot run failed: '+str(directory))
                record = {'repeat': repeat, 'nodes': count, 'ranks': 1,
                          'command': command, 'input_sha256': sha(input_path),
                          'revision': revision, 'artifacts': str(directory),
                          'wall_seconds': time.monotonic()-started,
                          **timed_metrics(directory/'time.txt')}
                step_match = re.search(r'time=([^ ]+) cycle=(\d+)', process.stdout)
                require(step_match and int(step_match[2]) == STEPS, 'wrong pilot step count')
                record.update(fluid_steps=int(step_match[2]), final_time=float(step_match[1]))
                cpu = re.search(r'cpu time used\s*=\s*([0-9.eE+-]+)', process.stdout)
                record['reported_simulation_cpu_seconds'] = float(cpu[1]) if cpu else None
                record['fluid_identity'] = fluid_identity(
                    sorted((directory/'bin').glob('*.hydro_u.*.bin'))[-1])
                if count:
                    first = reader.read(directory/'pilot.curves.initial.txt')
                    last = reader.read(directory/'pilot.curves.final.txt')
                    curve = last['curves'][0]
                    require(last['cycle'] == STEPS and len(curve['nodes']) == count,
                            'curve activity/count mismatch')
                    record.update(node_samples=curve['samples']-first['curves'][0]['samples'],
                                  inserted=curve['inserted'], migrations=curve['migrations'],
                                  curve_reason=curve['reason'])
                    require(record['node_samples'] > 0 and record['inserted'] == 0,
                            'inactive sampling or unexpected adaptation')
                else:
                    require(not list(directory.glob('pilot.curves.*')),
                            'feature-off run wrote curve data')
                    record.update(node_samples=0, inserted=0, migrations=0)
                record['output_bytes'] = sum(path.stat().st_size for path in
                    list((directory/'bin').glob('*.bin'))+list(directory.glob('pilot.curves.*')))
                records.append(record)
                write(directory/'measurement.json', record)
                write(artifacts/'results.json', metadata)
        identities = {r['fluid_identity']['sha256'] for r in records}
        require(len(identities) == 1, 'feature on/off changed native fluid output payload')
        metadata['baseline_fluid_identity'] = records[0]['fluid_identity']
        baseline = statistics.median(r['wall_seconds'] for r in records if not r['nodes'])
        summary = []
        for count in COUNTS:
            group = [r for r in records if r['nodes'] == count]
            walls = [r['wall_seconds'] for r in group]
            rss = [r['peak_rss_bytes'] for r in group]
            native_walls = [r['wall_seconds_native_time'] for r in group]
            summary.append({'nodes': count, 'wall_median_seconds': statistics.median(walls),
                            'wall_min_seconds': min(walls), 'wall_max_seconds': max(walls),
                            'wall_median_ratio_to_off': statistics.median(walls)/baseline,
                            'native_wall_median_seconds': statistics.median(native_walls),
                            'native_wall_min_seconds': min(native_walls),
                            'native_wall_max_seconds': max(native_walls),
                            'rss_median_bytes': statistics.median(rss),
                            'rss_min_bytes': min(rss), 'rss_max_bytes': max(rss),
                            'output_bytes': group[0]['output_bytes']})
        metadata.update(status='pass', summary=summary)
    except Exception as error:
        metadata.update(status='fail', failure=str(error))
        write(artifacts/'results.json', metadata)
        raise
    write(artifacts/'results.json', metadata)
    print(json.dumps({'status': 'pass', 'summary': summary, 'artifacts': str(artifacts)}))


if __name__ == '__main__':
    main()
