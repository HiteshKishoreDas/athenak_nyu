"""Independent unforced HD controls for material curves (selected F02/F03/F05/F06).

Run the default built-in AthenaK executable with --executable and --artifacts.
The exact Sod oracle and acceptance budgets are saved before launching any case.
This is a finite-grid shock/rarefaction/contact check, not a shock-finder or
curvature qualification. Every failure retains its inputs, logs and measurements.
Only Python's standard library is required; simulations use production HD RK2,
PLM/HLLC and fluid primitive sampling, without a prescribed-field/source hook.
"""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

GAMMA = 1.4
FINAL_TIME = 0.1
GRIDS = (64, 128, 256, 512)
ANCHORS = 64  # even count avoids an initial anchor exactly on the discontinuity
SOD_BUDGET = 2e-3
TRANSLATION_BUDGET = 2e-12
SHOCK_EXCLUSION_CELLS = 3
ROOT = Path(__file__).resolve().parents[3]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def pressure_function(p, rho, pressure, velocity):
    """Exact ideal-gas Riemann wave function, independent of AthenaK solvers."""
    if p > pressure:
        a = 2 / ((GAMMA+1)*rho)
        b = (GAMMA-1)*pressure / (GAMMA+1)
        return (p-pressure)*math.sqrt(a/(p+b))
    sound = math.sqrt(GAMMA*pressure/rho)
    return 2*sound/(GAMMA-1)*((p/pressure)**((GAMMA-1)/(2*GAMMA))-1)


def sod_parameters():
    left, right = (1.0, 1.0, 0.0), (0.125, 0.1, 0.0)

    def residual(p):
        return pressure_function(p, *left) + pressure_function(p, *right)

    lo, hi = 0.0, 1.0
    require(residual(lo) < 0 < residual(hi), 'Sod pressure root not bracketed')
    for _ in range(100):
        mid = (lo+hi)/2
        if residual(mid) > 0:
            hi = mid
        else:
            lo = mid
    pstar = (lo+hi)/2
    ustar = (pressure_function(pstar, *right)-pressure_function(pstar, *left))/2
    a_left = math.sqrt(GAMMA)
    shock = math.sqrt(GAMMA*right[1]/right[0])*math.sqrt(
        (GAMMA+1)/(2*GAMMA)*(pstar/right[1])+(GAMMA-1)/(2*GAMMA))
    tail = ustar-a_left*pstar**((GAMMA-1)/(2*GAMMA))
    require(abs(residual(pstar)) < 2e-15, 'Sod bisection residual')
    require(abs(pstar-0.3031301780506468) < 2e-15, 'Sod pressure sanity check')
    require(abs(ustar-0.9274526200489499) < 2e-15, 'Sod velocity sanity check')
    return {'pstar': pstar, 'ustar': ustar, 'a_left': a_left,
            'shock_speed': shock, 'fan_head_speed': -a_left, 'fan_tail_speed': tail}


def sod_trajectory(x0, t, parameters):
    """Exact material position; shock arrival and fan transit use initial label."""
    require(x0 != 0 and t >= 0, 'ambiguous interface label or negative time')
    ustar = parameters['ustar']
    if x0 > 0:
        entry = x0/parameters['shock_speed']
        return x0+ustar*max(0, t-entry), ('right_star' if t > entry
                                         else 'right_unaffected')
    entry = -x0/parameters['a_left']
    if t <= entry:
        return x0, 'left_unaffected'
    k = 2/(GAMMA+1)
    q = k*parameters['a_left']/(1-k)
    coefficient = (x0-q*entry)/entry**k
    tail = parameters['fan_tail_speed']
    exit_time = (coefficient/(tail-q))**(1/(1-k))
    if t <= exit_time:
        return coefficient*t**k+q*t, 'fan'
    return tail*exit_time+ustar*(t-exit_time), 'left_star'


def initial_anchors():
    return [[-0.25+0.5*i/(ANCHORS-1), -0.2+0.4*i/(ANCHORS-1), 0.0]
            for i in range(ANCHORS)]


def polygon_length(points):
    return math.fsum(math.dist(a, b) for a, b in zip(points, points[1:]))


def frozen_oracle():
    parameters = sod_parameters()
    anchors = initial_anchors()
    expected, phases = [], []
    for point in anchors:
        x, phase = sod_trajectory(point[0], FINAL_TIME, parameters)
        expected.append([x, point[1], point[2]])
        phases.append(phase)
        # Analytic continuity and branch timing checks precede any candidate data.
        if point[0] < 0:
            entry = -point[0]/parameters['a_left']
            require(abs(sod_trajectory(point[0], entry, parameters)[0]-point[0]) < 1e-15,
                    'fan entry position mismatch')
        else:
            entry = point[0]/parameters['shock_speed']
            require(sod_trajectory(point[0], entry, parameters)[0] == point[0],
                    'shock entry position mismatch')
    require({'fan', 'left_star', 'right_star', 'left_unaffected',
             'right_unaffected'} <= set(phases), 'oracle does not cover all wave regions')
    return {'gamma': GAMMA, 'time': FINAL_TIME, 'parameters': parameters,
            'initial_anchors': anchors, 'sod_final_anchors': expected, 'phases': phases,
            'contact_uniform_final_anchors': [[x[0]+0.4*FINAL_TIME, x[1], x[2]]
                                             for x in anchors],
            'reference_polygon_length': polygon_length(expected),
            'grids': list(GRIDS), 'ny': 8, 'anchor_count': ANCHORS,
            'acceptance': {'sod_final_max_position_error': SOD_BUDGET,
                           'sod_error_must_decrease_each_refinement': True,
                           'contact_uniform_translation_error': TRANSLATION_BUDGET,
                           'contact_pressure_momentum_ratio_error': TRANSLATION_BUDGET,
                           'shock_exclusion_cells': SHOCK_EXCLUSION_CELLS},
            'limits': ['No second-order assertion at discontinuities.',
                       'Shock-near exclusions are separately reported, not removed '
                       'from the final maximum-error acceptance test.',
                       'Length reference is a polygon through exact initial anchors '
                       'advected analytically, not continuous-curve quadrature.',
                       'Finite-grid shock width limits local stretch and curvature; '
                       'physical curvature claims are excluded.',
                       'Periodic second Riemann jump remains beyond tested anchors '
                       'at t=0.1; this is not a general periodic Sod solution.',
                       'Fixed sampling (adaptive=false) is not an adaptation '
                       'or geometric-tolerance qualification.']}


def input_text(nx, kind):
    contact = kind == 'contact'
    uniform = kind == 'uniform'
    dl, dr = 1.0, (1.0 if uniform else 0.125)
    pl, pr = 1.0, (1.0 if contact or uniform else 0.1)
    velocity = 0.4 if contact or uniform else 0.0
    return f'''<job>
basename = physical
<mesh>
nghost = 2
nx1 = {nx}
x1min = -0.5
x1max = 0.5
ix1_bc = periodic
ox1_bc = periodic
nx2 = 8
x2min = -0.5
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic
nx3 = 1
x3min = -0.5
x3max = 0.5
ix3_bc = periodic
ox3_bc = periodic
<meshblock>
nx1 = 32
nx2 = 8
nx3 = 1
<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.4
tlim = {FINAL_TIME}
nlim = -1
ndiag = 100000
<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = {GAMMA}
<problem>
pgen_name = shock_tube
shock_dir = 1
xshock = 0
dl = {dl}
dr = {dr}
pl = {pl}
pr = {pr}
ul = {velocity}
ur = {velocity}
vl = 0
vr = 0
wl = 0
wr = 0
<material_curves>
enabled = true
ncurves = 1
max_nodes = 10000
spacing_cells = 1000000
max_spacing = 1000
max_angle = 3
chord_error = 1000
adaptive = false
output_dcycle = 0
<curve1>
id = 1
type = line
closed = false
count = {ANCHORS}
x1 = -0.25
x2 = -0.2
x3 = 0
end1 = 0.25
end2 = 0.2
end3 = 0
<output1>
file_type = tab
variable = hydro_w
slice_x2 = 0
slice_x3 = 0
data_format = %24.17e
dt = {FINAL_TIME}
<output2>
file_type = tab
variable = hydro_u
slice_x2 = 0
slice_x3 = 0
data_format = %24.17e
dt = {FINAL_TIME}
'''


def digest(path):
    hasher = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            hasher.update(block)
    return hasher.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def table(path, nx):
    lines = path.read_text().splitlines()
    require(len(lines) > 2 and lines[0].startswith('# Athena++'), 'bad fluid table')
    time_match = re.search(r'time=([^ ]+)', lines[0])
    require(time_match and abs(float(time_match[1])-FINAL_TIME) < 1e-13,
            'fluid table time mismatch')
    columns = lines[1].removeprefix('#').split()
    rows = []
    for line in lines[2:]:
        if not line.strip() or line.startswith('#'):
            continue
        values = [float(token) for token in line.split()]
        require(len(values) == len(columns) and all(math.isfinite(v) for v in values),
                'fluid table column/nonfinite error')
        rows.append(dict(zip(columns, values)))
    require(len(rows) == nx, 'fluid table does not contain every x cell')
    rows.sort(key=lambda row: row['x1v'])
    return rows


def build_metadata(executable, declared_revision):
    build = executable.parent.parent
    metadata = {'executable': str(executable), 'executable_sha256': digest(executable),
                'hardware': platform.platform(), 'python': sys.version,
                'ranks': 1, 'threads': 1, 'revision': declared_revision,
                'numerics': {'integrator': 'rk2', 'reconstruct': 'plm',
                             'rsolver': 'hllc', 'cfl': 0.4, 'gamma': GAMMA,
                             'sources': 'none', 'adaptive_curves': False}}
    config = build/'config.hpp'
    require(config.is_file(), 'missing generated executable build config')
    config_text = config.read_text()
    require('#define USER_PROBLEM_ENABLED 0' in config_text,
            'physical tests require the built-in unforced shock_tube pgen')
    require('#define SINGLE_PRECISION_ENABLED 0' in config_text,
            'frozen budgets require double precision')
    metadata['build_config'] = config_text
    flags = build/'src/CMakeFiles/athena.dir/flags.make'
    if flags.is_file():
        metadata['compiler_flags'] = flags.read_text()
    cache = build/'CMakeCache.txt'
    if cache.is_file():
        text = cache.read_text()
        metadata['build_settings'] = [line for line in text.splitlines()
                                      if re.match(r'(CMAKE_BUILD_TYPE|Athena_ENABLE|'
                                                  r'Kokkos_ENABLE|CMAKE_CXX_COMPILER):', line)]
        match = re.search(r'^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$', text, re.MULTILINE)
        if match:
            source = Path(match[1])
            result = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'],
                                    capture_output=True, text=True)
            metadata['source_revision_observed'] = result.stdout.strip()
            if metadata['revision'] is None:
                metadata['revision'] = metadata['source_revision_observed']
    return metadata


def measure(initial, final, expected, phases=None, shock_position=None, dx=None):
    require(initial['time'] == 0 and abs(final['time']-FINAL_TIME) <= 1e-14,
            'initial/final material time mismatch')
    require(initial['curve_count'] == final['curve_count'] == 1, 'wrong curve count')
    first, last = initial['curves'][0], final['curves'][0]
    require(len(first['nodes']) == len(last['nodes']) == ANCHORS,
            'fixed-anchor fixture changed node count')
    require(first['inserted'] == last['inserted'] == 0, 'fixed sampling inserted nodes')
    require(last['samples'] > first['samples'] and final['cycle'] > 0,
            'material sampling/fluid evolution path was inactive')
    require(last['reason'] == 0, 'structurally invalid fixed curve')
    current = {node['id']: node for node in last['nodes']}
    actual, errors = [], []
    for i, anchor in enumerate(first['nodes']):
        require(anchor['anchor'] == 1 and anchor['birth'] == 0,
                'fixture contains a born sample')
        require(math.dist(anchor['x'], initial_anchors()[i]) < 2e-15,
                'candidate initial geometry differs from frozen oracle')
        node = current[anchor['id']]
        require(node['a'] == anchor['a'] and node['birth'] == anchor['birth'] and
                node['anchor'] == 1, 'changed anchor label/history')
        actual.append(node['x'])
        errors.append(math.dist(node['x'], expected[i]))
    groups = {}
    if phases is not None:
        for phase in sorted(set(phases)):
            selected = [errors[i] for i, p in enumerate(phases) if p == phase]
            groups[phase] = {'count': len(selected), 'max_error': max(selected)}
    near = ([i for i, x in enumerate(expected)
             if abs(x[0]-shock_position) <= SHOCK_EXCLUSION_CELLS*dx]
            if shock_position is not None else [])
    far = [errors[i] for i in range(ANCHORS) if i not in near]
    require(far, 'shock exclusion removed all anchors')
    reference_length, measured_length = polygon_length(expected), polygon_length(actual)
    return {'max_position_error': max(errors),
            'rms_position_error': math.sqrt(math.fsum(e*e for e in errors)/ANCHORS),
            'y_displacement_error': max(abs(x[1]-initial_anchors()[i][1])
                                        for i, x in enumerate(actual)),
            'shock_near_count': len(near), 'shock_near_indices': near,
            'shock_near_exclusion_cells': SHOCK_EXCLUSION_CELLS,
            'max_position_error_excluding_shock_near': max(far), 'regions': groups,
            'reference_polygon_length': reference_length,
            'measured_polygon_length': measured_length,
            'polygon_length_absolute_error': abs(measured_length-reference_length),
            'polygon_length_relative_error': abs(measured_length-reference_length)
                                               / reference_length,
            'stored_length_reconstruction_error': abs(last['length']-measured_length),
            'fluid_steps': final['cycle'], 'node_samples': last['samples']-first['samples'],
            'nodes': len(last['nodes']), 'curve_reason': last['reason'],
            'anchor_errors': errors, 'expected_anchors': expected, 'actual_anchors': actual}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--artifacts', type=Path, required=True)
    parser.add_argument('--reader', type=Path,
                        default=ROOT/'vis/python/material_curves.py')
    parser.add_argument('--revision')
    parser.add_argument('--allow-unqualified-accelerator', action='store_true')
    args = parser.parse_args()
    executable, artifacts = args.executable.resolve(), args.artifacts.resolve()
    require(executable.is_file(), 'missing executable')
    require(not artifacts.exists(), 'artifact directory already exists; preserve old runs')
    artifacts.mkdir(parents=True)
    # This freeze happens before any executable invocation or candidate snapshot read.
    oracle = frozen_oracle()
    write_json(artifacts/'oracle.json', oracle)
    metadata = build_metadata(executable, args.revision)
    write_json(artifacts/'environment.json', metadata)
    spec = importlib.util.spec_from_file_location('material_curves_reader', args.reader)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    records = []
    status = {'status': 'running', 'records': records, 'artifacts': str(artifacts),
              'allow_unqualified_accelerator': args.allow_unqualified_accelerator,
              'cases': ['F02 selected contact', 'F03 selected planar Sod shock',
                        'F05 selected Sod rarefaction', 'F06 oblique string crossing']}
    try:
        for kind, nx in [('sod', nx) for nx in GRIDS] + [('contact', 128), ('uniform', 128)]:
            directory = artifacts/f'{kind}_{nx}'
            directory.mkdir()
            text = input_text(nx, kind)
            if args.allow_unqualified_accelerator:
                text += ('\n<material_curves>\n'
                         'allow_unqualified_accelerator = true\n')
            input_path = directory/'input.athinput'
            input_path.write_text(text)
            command = [str(executable), '-i', str(input_path)]
            record = {'name': directory.name, 'nx': nx, 'ny': 8, 'dx': 1/nx,
                      'command': command, 'input_sha256': digest(input_path),
                      'revision': metadata['revision'], 'ranks': 1, 'threads': 1,
                      'artifacts': str(directory)}
            records.append(record)
            write_json(artifacts/'results.json', status)
            require(digest(executable) == metadata['executable_sha256'],
                    'executable changed during qualification')
            started = time.monotonic()
            process = subprocess.run(command, cwd=directory, capture_output=True,
                                     text=True, timeout=120)
            record.update(returncode=process.returncode, wall_seconds=time.monotonic()-started)
            (directory/'stdout.txt').write_text(process.stdout)
            (directory/'stderr.txt').write_text(process.stderr)
            require(process.returncode == 0, 'AthenaK failed: '+str(directory))
            initial = reader.read(directory/'physical.curves.initial.txt')
            final = reader.read(directory/'physical.curves.final.txt')
            expected = (oracle['sod_final_anchors'] if kind == 'sod' else
                        oracle['contact_uniform_final_anchors'])
            record.update(measure(initial, final, expected,
                                  oracle['phases'] if kind == 'sod' else None,
                                  oracle['parameters']['shock_speed']*FINAL_TIME
                                  if kind == 'sod' else None, 1/nx))
            if kind != 'sod':
                primitive = table(sorted((directory/'tab').glob('*.hydro_w.*.tab'))[-1], nx)
                conserved = table(sorted((directory/'tab').glob('*.hydro_u.*.tab'))[-1], nx)
                pressure_error = ratio_error = conserved_pressure_error = 0.0
                for prim, cons in zip(primitive, conserved):
                    require(prim['x1v'] == cons['x1v'] and cons['dens'] > 0,
                            'fluid primitive/conserved row mismatch')
                    pressure_error = max(pressure_error, abs((GAMMA-1)*prim['eint']-1))
                    ratio_error = max(ratio_error, abs(cons['mom1']/cons['dens']-0.4),
                                      abs(prim['velx']-0.4))
                    kinetic = math.fsum(cons[m]*cons[m] for m in ('mom1', 'mom2', 'mom3'))
                    pressure = (GAMMA-1)*(cons['ener']-kinetic/(2*cons['dens']))
                    conserved_pressure_error = max(conserved_pressure_error, abs(pressure-1))
                record.update(primitive_pressure_error=pressure_error,
                              conserved_pressure_error=conserved_pressure_error,
                              momentum_density_ratio_error=ratio_error,
                              density_range=[min(p['dens'] for p in primitive),
                                             max(p['dens'] for p in primitive)])
            write_json(directory/'measurement.json', record)
            write_json(artifacts/'results.json', status)
            require(record['y_displacement_error'] <= TRANSLATION_BUDGET,
                    'unexpected transverse material displacement')
            if kind != 'sod':
                require(record['max_position_error'] <= TRANSLATION_BUDGET,
                        kind+' translation exceeds frozen budget')
                require(max(record['primitive_pressure_error'],
                            record['conserved_pressure_error'],
                            record['momentum_density_ratio_error']) <= TRANSLATION_BUDGET,
                        kind+' pressure/momentum-ratio control failed')
        errors = [r['max_position_error'] for r in records if r['name'].startswith('sod_')]
        status.update(sod_max_errors=errors,
                      sod_refinement_ratios=[a/b for a, b in zip(errors, errors[1:])])
        require(all(a > b for a, b in zip(errors, errors[1:])),
                'Sod material-position error did not decrease at every refinement')
        require(errors[-1] <= SOD_BUDGET, 'finest Sod anchor error exceeds frozen 2e-3 budget')
        status['status'] = 'pass'
    except Exception as error:
        status.update(status='fail', failure=str(error), failure_type=type(error).__name__)
        write_json(artifacts/'results.json', status)
        raise
    write_json(artifacts/'results.json', status)
    print(json.dumps({'status': 'pass', 'sod_errors': status['sod_max_errors'],
                      'contact_error': records[-2]['max_position_error'],
                      'uniform_error': records[-1]['max_position_error'],
                      'artifacts': str(artifacts)}))


if __name__ == '__main__':
    main()
