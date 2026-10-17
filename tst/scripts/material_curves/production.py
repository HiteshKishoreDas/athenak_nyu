"""Bounded production-path curve tests using the prescribed mesh-field pgen.

Compile with -DPROBLEM=material_curve_test -DAthena_ENABLE_MPI=ON. Direct invocation
requires --executable and --output-dir. No output directory is deleted by this test.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

if __package__:
    from . import analytic
    from ..utils import athena
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from scripts.material_curves import analytic
    from scripts.utils import athena

ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT/'inputs/tests/material_curves/prescribed.athinput'
_results = []


def load_reader(path):
    spec = importlib.util.spec_from_file_location('curve_output_reader', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exercise(executable, output_dir, reader, suite='core', launcher=None, replay=False,
             allow_unqualified_accelerator=False):
    records = []
    output_dir.mkdir(parents=True, exist_ok=True)
    require(executable.is_file(), 'missing AthenaK executable '+str(executable))

    def case(name, overrides=None, ranks=1, restart=None, expect_curves=True, reject=None):
        directory = output_dir/name
        if replay:
            require(directory.is_dir(), 'missing replay artifacts: '+str(directory))
        else:
            require(not directory.exists(), 'result directory already exists: '+str(directory))
            directory.mkdir()
        arguments = {"job/basename": "material"}
        arguments.update(overrides or {})
        if allow_unqualified_accelerator:
            arguments['material_curves/allow_unqualified_accelerator'] = 'true'
        # AthenaK command-line overrides cannot add parameters absent from the input.
        # Repeated blocks are explicitly supported by its input reader and overwrite
        # existing fields, so append each bounded fixture override to the copied input.
        text = INPUT.read_text()
        for key, value in arguments.items():
            block, field = key.split('/')
            text += '\n<{}>\n{} = {}\n'.format(block, field, value)
        actual_input = directory/'input.athinput'
        if replay:
            require(actual_input.read_text() == text, 'replay input differs from frozen fixture')
        else:
            actual_input.write_text(text)
        command = [str(executable), '-i', str(actual_input)]
        if restart is not None:
            command += ['-r', str(restart)]
        if ranks > 1:
            command = (shlex.split(launcher or os.environ.get('CURVES_MPI_LAUNCHER',
                                                             'mpiexec'))
                       + ['-n', str(ranks)] + command)
        if replay:
            saved = json.loads((directory/'command.json').read_text())
            require(saved['command'] == command, 'replay command differs from recorded run')
            stdout = (directory/'stdout.txt').read_text()
            stderr = (directory/'stderr.txt').read_text()
            returncode = saved['returncode']
        else:
            process = subprocess.run(command, cwd=directory, capture_output=True,
                                     text=True, timeout=120)
            stdout, stderr, returncode = process.stdout, process.stderr, process.returncode
            (directory/'stdout.txt').write_text(stdout)
            (directory/'stderr.txt').write_text(stderr)
        metadata = {'name': name, 'command': command, 'ranks': ranks,
                    'input_sha256': hashlib.sha256(actual_input.read_bytes()).hexdigest(),
                    'artifacts': str(directory), 'returncode': returncode, 'replay': replay}
        if not replay:
            metadata['executable_sha256'] = hashlib.sha256(executable.read_bytes()).hexdigest()
            (directory/'command.json').write_text(json.dumps(metadata, indent=2)+'\n')
        if reject is not None:
            require(returncode != 0 and reject in stdout+stderr,
                    name+' did not explicitly reject '+reject)
            metadata.update(expected_rejection=reject, case='R03')
            records.append(metadata)
            return None, None, stdout
        require(returncode == 0, name+' executable failed; see '+str(directory))
        initial_path = directory/'material.curves.initial.txt'
        final_path = directory/'material.curves.final.txt'
        if not expect_curves:
            require(not initial_path.exists() and not final_path.exists(),
                    'disabled feature unexpectedly wrote curve snapshots')
            records.append(metadata)
            return None, None, stdout
        require((initial_path.is_file() or restart is not None) and final_path.is_file(),
                name+' missing initial/final curve snapshots')
        initial = (reader.read(initial_path) if restart is None else
                   reader.read(str(restart)+'.curves', fluid_path=restart))
        final = reader.read(final_path)
        require(final['time'] > initial['time'], name+' did not evolve')
        require(final['curve_count'] == 1, name+' unexpected curve count')
        curve = final['curves'][0]
        require(curve['samples'] > initial['curves'][0]['samples'],
                name+' did not sample fluid primitives during this run')
        metadata.update(time=final['time'], cycle=final['cycle'], nodes=len(curve['nodes']),
                        samples=curve['samples'], inserted=curve['inserted'],
                        migrations=curve['migrations'], reason=curve['reason'])
        records.append(metadata)
        return initial, final, stdout

    def positions(initial, final, mapping, tolerance=analytic.POSITION_BUDGET):
        start = initial['curves'][0]
        end = final['curves'][0]
        by_id = {n['id']: n for n in end['nodes']}
        actual, expected = [], []
        for node in start['nodes']:
            if not node['anchor']:
                continue
            require(node['id'] in by_id, 'lost initial material anchor')
            current = by_id[node['id']]
            require(current['anchor'] and current['a'] == node['a']
                    and current['birth'] == node['birth'], 'changed anchor identity/history')
            actual.append(current['x'])
            expected.append(mapping(node['x'], final['time']))
        error = analytic.check_positions(actual, expected, tolerance)
        records[-1]['position_error'] = error
        records[-1]['position_budget'] = tolerance
        return error

    def diagnostic_match(final):
        curve = final['curves'][0]
        geometric = curve['reconstructed']
        require(abs(curve['length']-geometric['length']) <= 2e-11,
                'stored and independently reconstructed lengths disagree')
        require(curve['reason'] == 0, 'unexpected unresolved curve')

    if suite in ('core', 'all'):
        initial, final, _ = case('zero', {'problem/flow': 'zero', 'problem/v1': 0,
                                         'problem/v2': 0})
        positions(initial, final, lambda x, t: x, analytic.ROUNDING_TOL)
        diagnostic_match(final)
        initial, final, _ = case('translation_wrap', {'problem/v1': 2,
                                                     'problem/v2': -1, 'time/tlim': 2.2})
        positions(initial, final, lambda x, t: analytic.translation(x, t, (2, -1, 0)),
                  analytic.ROUNDING_TOL)
        require(any(abs(i) > 0 for n in final['curves'][0]['nodes'] for i in n['image']),
                'translation test did not cross periodic boundaries')
        diagnostic_match(final)
        initial, final, _ = case('acceleration', {'problem/flow': 'acceleration',
                                                 'problem/a1': 0.3, 'time/tlim': 0.8})
        positions(initial, final, lambda x, t: analytic.acceleration(
            x, t, (0.2, -0.1, 0), (0.3, 0, 0)), analytic.ROUNDING_TOL)
        initial, final, _ = case('physical_acceleration', {'problem/user_srcs': 'false',
            'problem/flow': 'translation', 'hydro_srcterms/const_accel': 'true',
            'hydro_srcterms/const_accel_val': 0.3, 'hydro_srcterms/const_accel_dir': 1, 'time/tlim': 0.8})
        positions(initial, final, lambda x, t: analytic.acceleration(
            x, t, (0.2, -0.1, 0), (0.3, 0, 0)), analytic.ROUNDING_TOL)
        for flow in ('sinusoidal', 'rotation', 'strain'):
            errors = []
            for level, cfl in enumerate((0.4, 0.2, 0.1, 0.05)):
                arguments = {'problem/flow': flow, 'problem/v1': 0, 'problem/v2': 0,
                             'time/tlim': 0.8, 'time/cfl_number': cfl,
                             'problem/a1': 1.0, 'problem/amplitude': 0.5,
                             'problem/s1': 0.2, 'problem/s2': -0.2}
                initial, final, _ = case(flow+'_time_'+str(level), arguments)
                maps = {'sinusoidal': lambda x, t: analytic.sinusoidal(x, t, (1, 0, 0), 1),
                        'rotation': lambda x, t: analytic.rotation(x, t, 0.5),
                        'strain': lambda x, t: analytic.strain(x, t, (0.2, -0.2, 0))}
                errors.append(positions(initial, final, maps[flow]))
            records.append({'case': 'E02', 'name': flow+'_temporal_order',
                            'errors': errors, 'slopes': analytic.check_order(errors),
                            'minimum_order': analytic.MIN_RK2_ORDER,
                            'budget': analytic.POSITION_BUDGET})
        initial, final, _ = case('shear', {'problem/flow': 'shear', 'problem/v1': 0,
                                          'problem/v2': 0, 'problem/amplitude': 0.4})
        positions(initial, final, lambda x, t: analytic.shear(x, t, 0.4),
                  analytic.ROUNDING_TOL)
        initial, final, _ = case('nonuniform_density', {'problem/density_amplitude': 0.5})
        positions(initial, final, lambda x, t: analytic.translation(x, t, (0.2, -0.1, 0)),
                  analytic.ROUNDING_TOL)
        length_errors = []
        for count in (8, 16, 32, 64):
            initial, final, _ = case('circle_'+str(count), {'curve1/type': 'circle',
                'curve1/count': count, 'curve1/x1': 0, 'curve1/x2': 0,
                'curve1/radius': 0.7, 'problem/v1': 0, 'problem/v2': 0})
            curve = final['curves'][0]
            exact = analytic.circle_geometry(0.7, count)
            require(abs(curve['length']-exact['polygon_length']) <= 2e-11,
                    'circle polygon length is wrong')
            require(all(n['curvature_mask'] == 0 and
                        abs(n['curvature']-exact['curvature']) <= 2e-10
                        for n in curve['nodes']), 'circle curvature is wrong')
            length_errors.append(analytic.length_error(curve['length'], exact['length']))
            diagnostic_match(final)
        records.append({'case': 'E04', 'name': 'circle_geometry_order',
                        'errors': length_errors,
                        'slopes': analytic.check_order(length_errors, budget=5e-4),
                        'minimum_order': 1.8, 'budget': 5e-4})
        spatial_errors = []
        for nx in (16, 32, 64, 128):
            initial, final, _ = case('nonlinear_spatial_'+str(nx), {
                'mesh/nx1': nx, 'mesh/nx2': nx, 'problem/flow': 'reversing_shear',
                'problem/v1': 0, 'problem/v2': 0, 'problem/amplitude': 0.2,
                'problem/wavenumber': math.pi/2, 'problem/frequency': 0,
                'curve1/x1': 0.137, 'curve1/end1': 0.137, 'curve1/count': 257,
                'curve1/x2': -1.137, 'curve1/end2': 1.137, 'time/tlim': 0.4})
            spatial_errors.append(positions(initial, final, lambda x, t:
                analytic.reversing_shear(x, t, 0.2, math.pi/2, 0), tolerance=0.1))
        records.append({'case': 'E03', 'name': 'nonlinear_spatial_order',
                        'errors': spatial_errors,
                        'slopes': analytic.check_order(spatial_errors, budget=1e-4),
                        'minimum_order': 1.8, 'budget': 1e-4})
    if suite in ('mpi', 'all'):
        finals = []
        for ranks in (1, 2, 4):
            initial, final, _ = case('rank_crossing_'+str(ranks), {
                'problem/v1': 1.0, 'problem/v2': 1.0, 'time/tlim': 0.8,
                'curve1/x1': -0.003, 'curve1/end1': -0.003,
                'curve1/x2': -0.003, 'curve1/end2': 0.797}, ranks)
            positions(initial, final, lambda x, t: analytic.translation(x, t, (1, 1, 0)),
                      analytic.ROUNDING_TOL)
            if ranks > 1:
                require(final['curves'][0]['migrations'] > 0,
                        'rank crossing test did not migrate nodes')
            finals.append(final)
        for final in finals[1:]:
            first = finals[0]['curves'][0]['nodes']
            second = final['curves'][0]['nodes']
            require([n['id'] for n in first] == [n['id'] for n in second],
                    'rank count changed fixed anchor IDs')
            analytic.check_positions([n['x'] for n in first], [n['x'] for n in second],
                                     analytic.ROUNDING_TOL)
        records.append({'case': 'P01', 'name': 'rank_parity', 'ranks': [1, 2, 4]})
    if suite in ('extended', 'all'):
        errors, bending_errors = [], []
        for count in (8, 16, 32, 64):
            points = [analytic.helix(-math.pi+2*math.pi*i/count, 0.7, 0.2)
                      for i in range(count+1)]
            coordinates = ';'.join(' '.join(format(x, '.17g') for x in p) for p in points)
            initial, final, _ = case('helix_'+str(count), {
                'curve1/type': 'polyline', 'curve1/closed': 'false',
                'curve1/coordinates': 'lifted', 'curve1/points': coordinates,
                'mesh/nx3': 16, 'meshblock/nx3': 16, 'time/nlim': 1,
                'problem/v1': 0, 'problem/v2': 0})
            exact = analytic.helix_geometry(0.7, 0.2, 2*math.pi)
            curve = final['curves'][0]
            errors.append(analytic.length_error(curve['length'], exact['length']))
            bending_errors.append(max(abs(n['curvature']-exact['curvature'])
                                      for n in curve['nodes'] if n['curvature_mask'] == 0))
            positions(initial, final, lambda x, t: x, analytic.ROUNDING_TOL)
            diagnostic_match(final)
        for name, values in (('helix_length', errors), ('helix_curvature', bending_errors)):
            records.append({'case': 'G03', 'name': name+'_order', 'errors': values,
                            'slopes': analytic.check_order(values, budget=5e-4),
                            'minimum_order': 1.8, 'budget': 5e-4})
        initial, final, _ = case('oblique_rotation', {'mesh/nx3': 16, 'meshblock/nx3': 16,
            'curve1/type': 'circle', 'curve1/x1': 0, 'curve1/x2': 0, 'curve1/radius': 0.7,
            'curve1/normal1': 1, 'curve1/normal2': 2, 'curve1/normal3': 3,
            'problem/flow': 'rotation', 'problem/axis1': 1, 'problem/axis2': 2,
            'problem/axis3': 3, 'problem/amplitude': 0.2, 'problem/v1': 0,
            'problem/v2': 0})
        positions(initial, final, lambda x, t: analytic.rotation(x, t, 0.2, (1, 2, 3)))
        diagnostic_match(final)
        initial, final, _ = case('winding', {'curve1/type': 'line', 'curve1/closed': 'true',
            'curve1/winding1': 1, 'curve1/x1': -2, 'curve1/end1': 2,
            'curve1/x2': 0.2, 'curve1/end2': 0.2, 'problem/v1': 1,
            'problem/v2': 0, 'time/tlim': 0.8})
        positions(initial, final, lambda x, t: analytic.translation(x, t, (1, 0, 0)),
                  analytic.ROUNDING_TOL)
        require(final['curves'][0]['winding'] == [1, 0, 0], 'lost winding')
        require(abs(final['curves'][0]['length']-4) < 2e-11, 'wrong winding closure length')
        diagnostic_match(final)
        adaptive = {'material_curves/adaptive': 'true', 'material_curves/max_spacing': 0.18,
                    'material_curves/max_angle': 0.2, 'material_curves/chord_error': 0.001,
                    'curve1/count': 8, 'problem/flow': 'strain', 'problem/s1': 1,
                    'problem/v1': 0, 'problem/v2': 0, 'time/tlim': 0.8}
        initial, final, _ = case('adaptive_birth', adaptive)
        positions(initial, final, lambda x, t: analytic.strain(x, t, (1, 0, 0)))
        curve = final['curves'][0]
        require(curve['inserted'] > 0 and any(n['birth'] > 0 and not n['anchor']
                                             for n in curve['nodes']),
                'adaptation fixture did not insert after initial time')
        require(sum(n['anchor'] for n in curve['nodes']) == 8, 'lost history anchors')
        require(max(n['edge_length'] for n in curve['nodes']) <= 0.18+2e-12,
                'spacing contract not restored')
        diagnostic_match(final)
        cap = dict(adaptive, **{'material_curves/max_nodes': 8})
        initial, final, _ = case('adaptive_cap', cap)
        positions(initial, final, lambda x, t: analytic.strain(x, t, (1, 0, 0)))
        curve = final['curves'][0]
        require(len(curve['nodes']) == 8 and curve['inserted'] == 0
                and curve['reason'] == 1 and curve['unresolved_since'] > 0,
                'cap did not retain coherent flagged state')
        reversing = {'problem/flow': 'reversing_shear', 'problem/v1': 0, 'problem/v2': 0,
                     'problem/amplitude': 0.2, 'problem/wavenumber': math.pi/2,
                     'problem/frequency': math.pi/2, 'time/tlim': 4}
        initial, final, _ = case('reversing_return', reversing)
        positions(initial, final, lambda x, t: analytic.reversing_shear(
            x, t, 0.2, math.pi/2, math.pi/2), analytic.ROUNDING_TOL)
        reversing['time/tlim'] = 1
        initial, final, _ = case('reversing_max_fold', reversing)
        positions(initial, final, lambda x, t: analytic.reversing_shear(
            x, t, 0.2, math.pi/2, math.pi/2), tolerance=1e-3)
        # Full paired fluid/curve checkpoints, including post-birth insertion state.
        uninterrupted = dict(adaptive)
        _, full, _ = case('restart_reference', uninterrupted, ranks=2)
        checkpoint = dict(adaptive, **{'time/tlim': 0.4, 'output1/file_type': 'rst',
                                     'output1/dt': 0.4})
        _, _, _ = case('restart_source', checkpoint, ranks=2)
        files = sorted((output_dir/'restart_source/rst').glob('*.rst'))
        require(files, 'no shared fluid checkpoint written')
        fluid = files[-1]
        resume = dict(adaptive, **{'output1/file_type': 'rst', 'output1/dt': 0.4})
        _, resumed, _ = case('restart_four_ranks', resume, ranks=4, restart=fluid)
        first, second = full['curves'][0], resumed['curves'][0]
        require([n['id'] for n in first['nodes']] == [n['id'] for n in second['nodes']],
                'restart changed curve allocator/inserted IDs')
        for p, q in zip(first['nodes'], second['nodes']):
            require(all(p[key] == q[key] for key in
                        ('id', 'prev', 'next', 'parent_left', 'parent_right',
                         'a', 'anchor', 'placement')), 'restart changed exact metadata')
            require(abs(p['birth']-q['birth']) <= 2e-12, 'restart changed birth time')
        error = analytic.check_positions([n['x'] for n in first['nodes']],
                                         [n['x'] for n in second['nodes']], 2e-11)
        records.append({'case': 'O02', 'name': 'restart_rank_change_parity',
                        'source_ranks': 2, 'resume_ranks': 4,
                        'position_error': error, 'budget': 2e-11})
    if suite in ('amr', 'all'):
        amr = {'mesh_refinement/refinement':
            'adaptive', 'mesh_refinement/num_levels': 2,
            'mesh_refinement/refinement_interval': 1,
            'mesh_refinement/max_nmb_per_rank': 128,
            'problem/refine_time': 0.1, 'problem/coarsen_time': 0.25,
            'amr_criterion0/method': 'user'}
        initial, final, stdout = case('refine_coarsen', amr)
        positions(initial, final, lambda x, t: analytic.translation(x, t, (0.2, -0.1, 0)),
                  analytic.ROUNDING_TOL)
        counts = re.findall(r'(\d+) MeshBlocks created, (\d+) deleted by AMR', stdout)
        require(counts and int(counts[-1][0]) > 0 and int(counts[-1][1]) > 0,
                'AMR fixture did not refine and coarsen')
        records[-1]['amr_created'], records[-1]['amr_deleted'] = map(int, counts[-1])
        adaptive_amr = dict(amr, **{'material_curves/adaptive': 'true',
            'material_curves/max_spacing': 0, 'material_curves/spacing_cells': 1,
            'material_curves/max_angle': 0.2, 'material_curves/chord_error': 0,
            'curve1/count': 8})
        initial, final, stdout = case('refine_coarsen_with_insertion', adaptive_amr)
        positions(initial, final, lambda x, t: analytic.translation(x, t, (0.2, -0.1, 0)),
                  analytic.ROUNDING_TOL)
        curve = final['curves'][0]
        require(curve['inserted'] > initial['curves'][0]['inserted']
                and any(n['birth'] >= 0.1 for n in curve['nodes']),
                'mesh refinement did not trigger actual curve insertion')
        require(sum(n['anchor'] for n in curve['nodes']) == 8, 'AMR changed history anchors')
        diagnostic_match(final)
    if suite in ('adaptation', 'all'):
        arguments = {'curve1/type': 'circle', 'curve1/count': 32, 'curve1/x1': 0,
                     'curve1/x2': 0, 'curve1/radius': 0.7,
                     'material_curves/adaptive': 'true', 'material_curves/max_spacing': 0.12,
                     'material_curves/max_angle': 0.2, 'material_curves/chord_error': 0.001,
                     'problem/flow': 'strain', 'problem/s1': 1, 'problem/v1': 0,
                     'problem/v2': 0, 'time/tlim': 0.8}
        initial, final, _ = case('strained_circle_placement', arguments)
        curve = final['curves'][0]
        start = initial['curves'][0]
        phase = math.atan2(start['nodes'][0]['x'][1], start['nodes'][0]['x'][0])
        p, q = start['nodes'][0]['x'], start['nodes'][1]['x']
        orientation = 1 if p[0]*q[1]-p[1]*q[0] > 0 else -1
        expected = [analytic.strain(analytic.circle(phase + orientation*2*math.pi
                    *n['a']/curve['period'], 0.7), final['time'], (1, 0, 0))
                    for n in curve['nodes']]
        error = analytic.check_positions([n['x'] for n in curve['nodes']], expected)
        require(curve['reason'] == 0, 'adaptive ellipse recorded unresolved geometry')
        by_id = {n['id']: n for n in curve['nodes']}
        actual_error, linear_error = [], []
        for node, exact in zip(curve['nodes'], expected):
            if node['birth'] <= 0 or node['placement'] != 2:
                continue
            left, right = by_id[node['parent_left']], by_id[node['parent_right']]
            right_label = right['a'] + (curve['period'] if right['a'] <= left['a'] else 0)
            fraction = (node['a']-left['a'])/(right_label-left['a'])
            linear = [a+fraction*(b-a) for a, b in zip(left['x'], right['x'])]
            actual_error.append(analytic.position_error([node['x']], [exact]))
            linear_error.append(analytic.position_error([linear], [exact]))
        require(actual_error and max(actual_error) < max(linear_error),
                'cubic placement did not improve on linear baseline at common labels')
        dense = [analytic.strain(analytic.circle(phase + orientation*2*math.pi*i/2048, 0.7),
                                 final['time'], (1, 0, 0)) for i in range(2048)]
        geometry_error = analytic.geometric_error([n['x'] for n in curve['nodes']], dense,
                                                   closed=True, subdivisions=8)
        require(geometry_error <= 0.001, 'adaptive ellipse missed chord-error budget')
        time = final['time']
        exact_length = analytic.integrate(lambda a: math.hypot(
            0.7*math.exp(time)*math.sin(a), 0.7*math.cos(a)), 0, 2*math.pi)
        relative_length = analytic.length_error(curve['length'], exact_length)
        require(relative_length <= 0.001, 'adaptive ellipse missed relative length budget')
        curvature_error = max(abs(node['curvature']-
            analytic.curvature(expected[i-1], expected[i], expected[(i+1) % len(expected)]))
            for i, node in enumerate(curve['nodes']))
        records.append({'case': 'D04', 'name': 'cubic_vs_linear_exact_ellipse',
                        'nodes': len(curve['nodes']), 'inserted': curve['inserted'],
                        'post_birth_cubic_nodes': len(actual_error),
                        'position_error': error, 'position_budget': analytic.POSITION_BUDGET,
                        'cubic_max_error': max(actual_error), 'linear_max_error': max(linear_error),
                        'geometric_error': geometry_error, 'geometric_budget': 0.001,
                        'relative_length_error': relative_length, 'length_budget': 0.001,
                        'curvature_error_at_exact_labels': curvature_error,
                        'curvature_coverage': 'reported against same-label exact circumcircles'})
        # An initial kink is deliberately unresolved and must retain that status even
        # after adding samples. Curvature output must not erase historical quality.
        initial, final, _ = case('persistent_kink', {'curve1/type': 'polyline',
            'curve1/coordinates': 'lifted', 'curve1/points': '-0.4 0 0;0 0 0;0 0.4 0',
            'material_curves/adaptive': 'true', 'material_curves/max_spacing': 0.1,
            'material_curves/max_angle': 0.2, 'material_curves/chord_error': 0.001,
            'problem/v1': 0, 'problem/v2': 0, 'time/nlim': 1})
        require(initial['curves'][0]['reason'] == final['curves'][0]['reason'] == 4
                and final['curves'][0]['unresolved_since'] == 0,
                'insertion erased initial nonsmooth unresolved status')
    if suite in ('regression', 'all'):
        initial, final, _ = case('shifted_ulp_seam', {
            'mesh/x1min': 1, 'mesh/x1max': 4, 'mesh/x2min': -0.7, 'mesh/x2max': 1.3,
            'mesh/x3min': -4, 'mesh/x3max': 6, 'curve1/x1': math.nextafter(1.0, -math.inf),
            'curve1/end1': 1.2, 'curve1/x2': 0.25, 'curve1/end2': 0.35,
            'problem/v1': 0.2, 'problem/v2': -0.1})
        positions(initial, final, lambda x, t: analytic.translation(x, t, (0.2, -0.1, 0)),
                  analytic.ROUNDING_TOL)
        require(initial['curves'][0]['nodes'][0]['image'][0] == -1,
                'shifted seam fixture did not start below lower periodic boundary')
        diagnostic_match(final)
        for name, arguments, message in (
                ('unsupported_rk3', {'time/integrator': 'rk3'}, 'only dynamic HD with RK2'),
                ('unsupported_reflecting', {'mesh/ix1_bc': 'reflect', 'mesh/ox1_bc': 'reflect'},
                 'active dimensions must have periodic'),
                ('invalid_seed', {'curve1/count': 1}, 'invalid seed sample count')):
            case(name, arguments, reject=message)
        fluid = {'problem/user_srcs': 'false', 'hydro_srcterms/const_accel': 'true',
                 'hydro_srcterms/const_accel_val': 0.3, 'hydro_srcterms/const_accel_dir': 1,
                 'output1/file_type': 'bin', 'output1/variable': 'hydro_u',
                 'output1/dt': 0.4}
        case('fluid_curves_on', fluid)
        case('fluid_curves_off', dict(fluid, **{'material_curves/enabled': 'false'}),
             expect_curves=False)
        def binary_data(path):
            with path.open('rb') as stream:
                for _ in range(32):
                    line = stream.readline()
                    if line.startswith(b'  header offset='):
                        offset = int(line.split(b'=')[1])
                        stream.seek(offset, 1)
                        return stream.read()
                raise ValueError('missing native binary output header offset')
        first = sorted((output_dir/'fluid_curves_on/bin').glob('*.bin'))
        second = sorted((output_dir/'fluid_curves_off/bin').glob('*.bin'))
        require(first and [p.name for p in first] == [p.name for p in second],
                'fluid output times/counts changed when curves enabled')
        for p, q in zip(first, second):
            require(binary_data(p) == binary_data(q), 'curves changed fluid binary payload')
        records.append({'case': 'R01', 'name': 'fluid_no_feedback',
                        'snapshots': len(first), 'comparison': 'bitwise float32 native output'})
    if suite in ('unsupported', 'all'):
        cooling = {'cooling/enabled': 'true', 'cooling/units': 'code',
            'cooling/cooling_model': 'powerlaw', 'cooling/heating_model': 'none',
            'cooling/cooling_density': 'mass_density',
            'cooling/cooling_powerlaw_axis': 'temperature',
            'cooling/cooling_powerlaw_axis_units': 'code',
            'cooling/cooling_powerlaw_value_units': 'code',
            'cooling/cooling_reference_axis': 1, 'cooling/cooling_reference_value': 0.01,
            'cooling/cooling_slope': -0.7}
        unsupported = [('standalone_cooling', cooling),
            ('frame_tracking', {'frame_tracking/enabled': 'true',
                                'frame_tracking/axes': 'x1'}),
            ('orbital_shearing_objects', {'shearing_box/qshear': 1.5,
                                         'shearing_box/omega0': 1}),
            ('rk3', {'time/integrator': 'rk3'}),
            ('kinematic', {'time/evolution': 'kinematic', 'hydro/rsolver': 'advect'}),
            ('restart_per_rank', {'output1/file_type': 'rst', 'output1/dt': 0.1,
                                  'output1/single_file_per_rank': 'true'}),
            ('restart_per_node', {'output1/file_type': 'rst', 'output1/dt': 0.1,
                                  'output1/single_file_per_node': 'true'})]
        messages = {'standalone_cooling': 'cooling source stages are unqualified',
            'frame_tracking': 'frame tracking is unqualified',
            'orbital_shearing_objects': 'shearing/orbital advection stages are unqualified',
            'rk3': 'only dynamic HD with RK2 is supported',
            'kinematic': 'only dynamic HD with RK2 is supported',
            'restart_per_rank': 'curves require shared fluid checkpoints',
            'restart_per_node': 'curves require shared fluid checkpoints'}
        for name, parameters in unsupported:
            # The native configuration must construct successfully with curves off.
            # A zero-cycle control isolates configuration validity from extra physics.
            control = dict(parameters, **{'material_curves/enabled': 'false', 'time/nlim': 0})
            case('native_control_'+name, control, expect_curves=False)
            case('target_unsupported_'+name, parameters,
                 reject='FATAL ERROR: material curves: '+messages[name])
    if suite in ('amr_restart', 'all'):
        parameters = {'mesh_refinement/refinement': 'adaptive',
            'mesh_refinement/num_levels': 2, 'mesh_refinement/refinement_interval': 5,
            'mesh_refinement/ncycle_check': 3, 'amr_criterion0/method': 'user',
            'mesh_refinement/max_nmb_per_rank': 128, 'problem/refine_time': 0.1,
            'problem/coarsen_time': 0.25, 'material_curves/adaptive': 'true',
            'material_curves/max_spacing': 0, 'material_curves/spacing_cells': 1,
            'material_curves/max_angle': 0.2, 'material_curves/chord_error': 0,
            'material_curves/output_dcycle': 1,
            'curve1/count': 8, 'output1/file_type': 'rst', 'output1/dt': -1,
            'output1/dcycle': 3}

        def history(name):
            paths = sorted((output_dir/name).glob('material.curves.[0-9]*.txt'))
            states = [reader.read(path) for path in paths]
            require(states, 'missing every-cycle AMR snapshots')
            require(len({s['cycle'] for s in states}) == len(states),
                    'duplicate cycle snapshots')
            return {s['cycle']: s for s in states}

        def compare_state(first, second):
            require(first['cycle'] == second['cycle'] and
                    abs(first['time']-second['time']) <= 2e-12,
                    'AMR restart changed cycle/time history')
            require(first['amr'] == second['amr'],
                    'AMR restart changed logical leaf/refinement clock history at cycle '+
                    str(first['cycle']))
            pcurve, qcurve = first['curves'][0], second['curves'][0]
            require(all(pcurve[key] == qcurve[key] for key in
                        ('id', 'closed', 'winding', 'period', 'reason', 'unresolved_since',
                         'next_sequence', 'inserted', 'samples')),
                    'AMR restart changed curve counters or persistent quality')
            require([n['id'] for n in pcurve['nodes']] ==
                    [n['id'] for n in qcurve['nodes']], 'AMR restart changed inserted IDs')
            for p, q in zip(pcurve['nodes'], qcurve['nodes']):
                require(all(p[key] == q[key] for key in
                            ('id', 'prev', 'next', 'parent_left', 'parent_right',
                             'a', 'anchor', 'placement')), 'AMR restart changed exact metadata')
                require(abs(p['birth']-q['birth']) <= 2e-12,
                        'AMR restart changed birth history for node {}: {} versus {}'.format(
                            p['id'], p['birth'], q['birth']))
            return analytic.check_positions([n['x'] for n in pcurve['nodes']],
                                            [n['x'] for n in qcurve['nodes']], 2e-11)

        for source_ranks, resume_ranks in ((2, 4), (1, 1)):
            prefix = 'amr_restart_'+str(source_ranks)+'_'
            _, reference, _ = case(prefix+'reference', parameters, ranks=source_ranks)
            _, interrupted, _ = case(prefix+'interrupted',
                dict(parameters, **{'time/nlim': 33}), ranks=source_ranks)
            require(interrupted['cycle'] == 33, 'interrupted AMR run missed selected cycle')
            reference_history = history(prefix+'reference')
            leaf_counts = [len(s['amr']['leaves']) for s in reference_history.values()]
            require(max(leaf_counts) > min(leaf_counts) and leaf_counts[-1] == 4,
                    'AMR history did not actually refine/coarsen')
            files = sorted((output_dir/(prefix+'interrupted')/'rst').glob('*.rst'))
            checkpoints = {reader.read(str(p)+'.curves', fluid_path=p)['cycle']: p
                           for p in files}
            for checkpoint_cycle in (30, 33):
                require(checkpoint_cycle in checkpoints, 'missing selected AMR checkpoint')
                checkpoint = checkpoints[checkpoint_cycle]
                saved = reader.read(str(checkpoint)+'.curves', fluid_path=checkpoint)
                require((saved['time'] < 0.25) == (checkpoint_cycle == 30),
                        'checkpoint did not straddle coarsening eligibility')
                error = compare_state(reference_history[checkpoint_cycle], saved)
                name = prefix+'resume_'+str(checkpoint_cycle)
                _, resumed, _ = case(name, dict(parameters, **{'time/nlim': -1}),
                                     ranks=resume_ranks, restart=checkpoint)
                resumed_history = history(name)
                expected_cycles = {c for c in reference_history if c > checkpoint_cycle}
                require(set(resumed_history) == expected_cycles,
                        'AMR restart lost/added per-cycle history')
                for cycle, state in resumed_history.items():
                    error = max(error, compare_state(reference_history[cycle], state))
                error = max(error, compare_state(reference, resumed))
                records.append({'case': 'Q10', 'name': name,
                    'source_ranks': source_ranks, 'resume_ranks': resume_ranks,
                    'checkpoint_cycle': checkpoint_cycle, 'checkpoint_time': saved['time'],
                    'refinement_interval': 5, 'ncycle_check': 3,
                    'history_cycles_compared': len(resumed_history)+1,
                    'position_error': error, 'budget': 2e-11,
                    'final_meshblocks': len(resumed['amr']['leaves']),
                    'logical_leaves_and_clocks': 'exact every-cycle equality'})
    require(records, 'zero production tests executed')
    for row in records:
        row['status'] = 'pass'
    return records


def run(**kwargs):
    global _results
    _results = exercise(Path.cwd()/'build/src/athena', Path.cwd()/'material-curve-results',
                        load_reader(ROOT/'vis/python/material_curves.py'))


def analyze():
    return bool(_results) and all(row['status'] == 'pass' for row in _results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--reader', type=Path, default=ROOT/'vis/python/material_curves.py')
    parser.add_argument('--suite', choices=('core', 'extended', 'mpi', 'amr',
                                           'adaptation', 'regression', 'unsupported', 'amr_restart', 'all'),
                        default='core')
    parser.add_argument('--launcher', default=os.environ.get('CURVES_MPI_LAUNCHER',
                                                            'mpiexec'))
    parser.add_argument('--allow-unqualified-accelerator', action='store_true',
                        help='opt in to accelerator experiments with unchanged accuracy gates')
    parser.add_argument('--replay', action='store_true',
                        help='reanalyze matching saved inputs/commands; do not execute AthenaK')
    args = parser.parse_args()
    result = {'schema': 'athenak-material-curves-production-v1',
              'implementation_execution': not args.replay, 'analysis_replay': args.replay,
              'suite': args.suite,
              'allow_unqualified_accelerator': args.allow_unqualified_accelerator}
    try:
        result['tests'] = exercise(args.executable.resolve(), args.output_dir.resolve(),
                                   load_reader(args.reader), args.suite, args.launcher,
                                   args.replay, args.allow_unqualified_accelerator)
        result['status'] = 'pass'
    except (ValueError, OSError, subprocess.TimeoutExpired) as error:
        result.update(status='fail', error=str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.output_dir/('replay-'+args.suite+'.json' if args.replay
                                   else 'results.json')
    result_path.write_text(json.dumps(result, indent=2,
                                                        allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'suite': args.suite,
                      'results': str(result_path)}))
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
