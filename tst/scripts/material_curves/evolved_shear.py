"""Independent stationary Euler shear check through unforced HD and curve sampling.

Compile PROBLEM=material_curve_test. This qualifies this smooth stationary shear,
not a general evolved-flow or turbulence accuracy envelope.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

if __package__:
    from . import analytic
    from .production import INPUT, ROOT, load_reader, require
    from ..utils import athena
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from scripts.material_curves import analytic
    from scripts.material_curves.production import INPUT, ROOT, load_reader, require
    from scripts.utils import athena

# Frozen before unforced candidate output reads. Unit characteristic length.
POSITION_BUDGET = 1e-4
MINIMUM_ORDER = 1.8
FLUID_BUDGET = 2e-12
_results = []


def exercise(executable, output_dir, allow_unqualified_accelerator=False):
    require(executable.is_file(), 'missing AthenaK executable')
    reader = load_reader(ROOT/'vis/python/material_curves.py')
    records, errors = [], []
    amplitude, wavenumber, time = 0.2, math.pi/2, 0.4
    for nx in (16, 32, 64, 128):
        directory = output_dir/str(nx)
        require(not directory.exists(), 'existing output directory '+str(directory))
        directory.mkdir(parents=True)
        options = {'job/basename': 'shear', 'mesh/nx1': nx, 'mesh/nx2': nx,
            'problem/user_srcs': 'false', 'problem/flow': 'reversing_shear',
            'problem/amplitude': amplitude, 'problem/wavenumber': wavenumber,
            'problem/frequency': 0, 'problem/v1': 0, 'problem/v2': 0,
            'problem/constant_pressure': 1, 'hydro/rsolver': 'hllc',
            'curve1/x1': 0.137, 'curve1/end1': 0.137, 'curve1/x2': -1.137,
            'curve1/end2': 1.137, 'curve1/count': 257, 'time/tlim': time,
            'output1/file_type': 'tab', 'output1/variable': 'hydro_w',
            'output1/slice_x1': 0.137, 'output1/data_format': '%24.17e',
            'output1/dt': time}
        if allow_unqualified_accelerator:
            options['material_curves/allow_unqualified_accelerator'] = 'true'
        text = INPUT.read_text()
        for key, value in options.items():
            block, field = key.split('/')
            text += '\n<{}>\n{} = {}\n'.format(block, field, value)
        actual_input = directory/'input.athinput'
        actual_input.write_text(text)
        command = [str(executable), '-i', str(actual_input)]
        completed = subprocess.run(command, cwd=directory, capture_output=True,
                                   text=True, timeout=120)
        (directory/'stdout.txt').write_text(completed.stdout)
        (directory/'stderr.txt').write_text(completed.stderr)
        require(completed.returncode == 0, 'unforced shear failed; see '+str(directory))
        initial = reader.read(directory/'shear.curves.initial.txt')
        final = reader.read(directory/'shear.curves.final.txt')
        require(final['time'] == time and final['cycle'] > 0, 'missing actual HD evolution')
        start, end = initial['curves'][0], final['curves'][0]
        require(end['samples'] == 2*final['cycle']*len(start['nodes']),
                'unexpected actual mesh sampler stage count')
        require(end['reason'] == 0 and end['inserted'] == 0, 'invalid fixed curve state')
        require([n['id'] for n in start['nodes']] == [n['id'] for n in end['nodes']],
                'unforced run changed material identity')
        exact = [analytic.reversing_shear(n['x'], time, amplitude, wavenumber, 0)
                 for n in start['nodes']]
        error = analytic.position_error([n['x'] for n in end['nodes']], exact)
        errors.append(error)
        same_labels = analytic.polyline_length(exact)
        reference_length = analytic.integrate(lambda y: math.hypot(
            amplitude*wavenumber*time*math.cos(wavenumber*y), 1), -1.137, 1.137)
        polygon_length_error = analytic.length_error(end['length'], same_labels)
        smooth_length_error = analytic.length_error(end['length'], reference_length)
        # Full sampled y profile at each output. Stored primitives use internal
        # energy, so derive pressure by P=(gamma-1)*eint independently.
        tables = sorted((directory/'tab').glob('*.tab'))
        require(len(tables) >= 2, 'missing initial/final fluid profiles')
        fluid_error = 0.0
        for table in tables:
            rows = [line.split() for line in table.read_text().splitlines()
                    if line and not line.startswith('#')]
            require(len(rows) == nx, 'fluid profile did not cover all y cells')
            for row in rows:
                y = float(row[-6])
                density, vx, vy, vz, eint = map(float, row[-5:])
                actual = (density, vx, vy, vz, (5/3-1)*eint)
                expected = (1, amplitude*math.sin(wavenumber*y), 0, 0, 1)
                differences = [abs(a-b) for a, b in zip(actual, expected)]
                require(all(math.isfinite(e) for e in differences), 'nonfinite fluid state')
                fluid_error = max(fluid_error, max(differences))
        require(fluid_error <= FLUID_BUDGET,
                'stationary unforced Euler shear drifted: '+str(fluid_error))
        record = {'case': 'F01', 'name': 'unforced_stationary_shear_'+str(nx),
            'nx': nx, 'time': time, 'cycle': final['cycle'], 'samples': end['samples'],
            'nodes': len(end['nodes']), 'position_error': error,
            'fluid_profile_error': fluid_error, 'fluid_budget': FLUID_BUDGET,
            'polygon_length_error': polygon_length_error,
            'smooth_length_error': smooth_length_error,
            'source_hook_enabled': False, 'rsolver': 'hllc', 'command': command,
            'input_sha256': hashlib.sha256(actual_input.read_bytes()).hexdigest(),
            'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
            'artifacts': str(directory), 'status': 'pass'}
        records.append(record)
        (directory/'result.json').write_text(json.dumps(record, indent=2)+'\n')
    slopes = analytic.check_order(errors, budget=POSITION_BUDGET, minimum=MINIMUM_ORDER)
    records.append({'case': 'E03', 'name': 'unforced_stationary_shear_order',
        'errors': errors, 'slopes': slopes, 'minimum_order': MINIMUM_ORDER,
        'position_budget': POSITION_BUDGET, 'status': 'pass'})
    return records


def run(**kwargs):
    global _results
    _results = exercise(Path.cwd()/'build/src/athena',
                        Path.cwd()/'material-curve-evolved-shear')


def analyze():
    return len(_results) == 5 and all(row['status'] == 'pass' for row in _results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--allow-unqualified-accelerator', action='store_true')
    args = parser.parse_args()
    result = {'schema': 'athenak-material-curves-evolved-shear-v1',
              'implementation_execution': True}
    try:
        result['tests'] = exercise(args.executable.resolve(), args.output_dir.resolve(),
                                   args.allow_unqualified_accelerator)
        result['status'] = 'pass'
    except (ValueError, OSError, subprocess.TimeoutExpired) as error:
        result.update(status='fail', error=str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir/'results.json'
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'results': str(path)}))
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
