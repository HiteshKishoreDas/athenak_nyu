#!/usr/bin/env python3
"""Kokkos detector parity plus independently constructed physical/smooth fixtures.

Usage: .venv/bin/python parity.py executable --reference frozen-package-parent
All O(1) float comparisons are predeclared rtol=atol=1e-12. Masks are exact.
This is an engineering port check, not human scientific approval.
"""
import argparse
import csv
import io
import json
import hashlib
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

RTOL = 1e-12
ATOL = 1e-12
# Exact dirty-source contract frozen during M0, rather than its unrelated HEAD.
REFERENCE_SHA256 = {
    'local_pass.py': '6f5ec8bc89c90ecb1ad2e0135fb8abe7a79f2b43badfc73356fd4ae3f1cc3857',
    'config.py': '993746b2dd5b878f4e8c455262c68044168338c25a83d837b07a7574e168972f',
    'derived.py': '3508dc08f1e87199d8c1f6e129cefdd2c6e6bd37674f697bdbae1f786740c18d',
    'gradients.py': 'adafc584ebf4e329d7df1db921e7b6f3de559d32028cc3ce8d1fd3e22678e4e4',
    'sampling.py': '05b5f152efa4243c44e8ef83d10625969c39b77d1082d433a02ac094cf5fbfe4',
    'mach.py': 'f96e05e4f19717fbf39be9b866f6235fb37e20c1cb218b32ee95d053c191ca9d',
    'masks.py': 'aadabb59ca8670220eda518f7d25795ca3139506fedb5d7141d870e5a8bdbc40',
    'criteria.py': '055ac7a9ce15c098c5128a1fa5e42484f9045f9157b67bf46c07d00904b12543',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if not (args.reference / 'shockit' / 'local_pass.py').is_file():
        parser.error('--reference must contain the exact shockit package')
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.unlink(missing_ok=True)
    for file, expected_hash in REFERENCE_SHA256.items():
        data = (args.reference / 'shockit' / file).read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        if actual != expected_hash:
            parser.error(f'reference changed: {file}, expected {expected_hash}, '
                         f'got {actual}')
    sys.path.insert(0, str(args.reference.resolve()))
    from shockit.config import ShockFinderConfig
    from shockit.fields import FluidCube
    from shockit.local_pass import analyze_cube_pass
    from shockit.derived import compute_gradients

    records = []
    executable = str(args.executable.resolve())
    subprocess.run([executable, '--self-check'], check=True)

    def run(name, rho, pressure, velocity, gamma=1.4, spacing=(1., 1., 1.),
            width=1, min_mach=1.1, mach_max=1000., grad_floor=1e-30,
            div_limit=0., density_floor=0., temperature_floor=0.,
            physical_mach=None, analytic=None, invalid=False,
            reference_error=None, energy_override=None):
        shape = rho.shape
        energy = pressure / (gamma - 1.) if energy_override is None else energy_override
        # Compare identical representable EOS pressure, not nominal fixture intent.
        p = (gamma - 1.) * energy
        cfg = ShockFinderConfig(
            min_mach=min_mach, mach_max=mach_max, shock_width_cells=width,
            grad_T_min=grad_floor, div_v_threshold=div_limit,
            upstream_density_floor=density_floor or None,
            upstream_temperature_floor=temperature_floor or None,
            reduce_to_centers=False)
        with tempfile.TemporaryDirectory(prefix='shock-parity-') as tmp:
            fixture = Path(tmp) / 'fixture.txt'
            with fixture.open('w') as f:
                f.write(' '.join(map(str, shape)) + '\n')
                f.write(f'{gamma:.17g} {min_mach:.17g} {mach_max:.17g} '
                        f'{grad_floor:.17g} {div_limit:.17g} {width} '
                        f'{density_floor:.17g} {temperature_floor:.17g} '
                        + ' '.join(f'{x:.17g}' for x in spacing) + '\n')
                for idx in np.ndindex(shape):
                    values = [rho[idx]] + [v[idx] for v in velocity] + [energy[idx]]
                    f.write(' '.join(f'{x:.17g}' for x in values) + '\n')
            fixture_hash = hashlib.sha256(fixture.read_bytes()).hexdigest()
            if args.report:
                retained = args.report.parent / 'fixtures'
                retained.mkdir(exist_ok=True)
                shutil.copyfile(fixture, retained / f'{name}.txt')
            proc = subprocess.run([executable, str(fixture)], text=True,
                                  capture_output=True, check=False)
            if proc.returncode:
                raise AssertionError(f'{name}: native exit {proc.returncode}: '
                                     f'{proc.stderr}')
            lines = proc.stdout.splitlines()
            header = next(line for line in lines if line.startswith('# i,j,k,'))
            lines = lines[lines.index(header) + 1:]
            rows = list(csv.DictReader(io.StringIO('\n'.join(lines)),
                                       fieldnames=header[2:].split(',')))
            native = {key: np.array([float(row[key]) for row in rows]).reshape(shape)
                      for key in rows[0]}
            if invalid:
                assert np.any(native['valid'] == 0), name
                assert not np.any(native['inside']), name
                if reference_error:
                    cube = FluidCube(rho, p, *velocity, *spacing, gamma=gamma)
                    with np.errstate(over='ignore', invalid='ignore'):
                        try:
                            analyze_cube_pass(cube, cfg)
                        except ValueError as error:
                            assert reference_error in str(error), str(error)
                        else:
                            raise AssertionError(f'{name}: reference must reject')
                records.append(dict(case=name, status='pass', fixture_sha256=fixture_hash,
                                    reference_error=reference_error, invalid_cells=int(
                    np.count_nonzero(native['valid'] == 0))))
                return
            cube = FluidCube(rho, p, *velocity, *spacing, gamma=gamma)
            ref = analyze_cube_pass(cube, cfg)
            grads = compute_gradients(cube, ref.derived.temperature, ref.derived.entropy)
            expected = dict(
                valid=np.ones(shape), candidate=ref.shock_zone_mask,
                inside=ref.full_shock_mask, rho=rho, pressure=p,
                temperature=ref.derived.temperature, divergence=ref.derived.div_v,
                normal_x=ref.normals.normal_x, normal_y=ref.normals.normal_y,
                normal_z=ref.normals.normal_z, rho_up=ref.jumps.density_up,
                rho_down=ref.jumps.density_down, pressure_up=ref.jumps.pressure_up,
                pressure_down=ref.jumps.pressure_down,
                temperature_up=ref.jumps.temperature_up,
                temperature_down=ref.jumps.temperature_down,
                mach_pressure=ref.mach_fields.mach_pressure,
                mach_temperature=ref.mach_fields.mach_temperature)
            for axis in 'xyz':
                for field in 'tsr':
                    field_name = f'grad_{field}{axis}'
                    expected[field_name] = getattr(grads, field_name)
            maxima = {}
            for key, truth in expected.items():
                got = native[key]
                if key in ('valid', 'candidate', 'inside'):
                    okay = got == truth
                else:
                    okay = np.isclose(got, truth, rtol=RTOL, atol=ATOL, equal_nan=True)
                if not np.all(okay):
                    idx = tuple(np.argwhere(~okay)[0])
                    raise AssertionError(f'{name}: {key} first cell {idx}: '
                                         f'native={got[idx]} reference={truth[idx]}')
                finite = np.isfinite(truth) & np.isfinite(got)
                maxima[key] = float(np.max(np.abs(got[finite] - truth[finite]))) \
                    if finite.any() else None
            if physical_mach is not None:
                cells = native['inside'].astype(bool)
                assert cells.any(), name
                np.testing.assert_allclose(native['mach_pressure'][cells],
                                           physical_mach, rtol=RTOL, atol=ATOL)
                valid_t = cells & np.isfinite(native['mach_temperature'])
                if physical_mach <= mach_max:
                    assert np.all(valid_t[cells]), name
                    np.testing.assert_allclose(native['mach_temperature'][valid_t],
                                               physical_mach, rtol=RTOL, atol=ATOL)
                if physical_mach > mach_max:
                    assert not valid_t.any(), name
            if analytic is not None:
                analytic(native)
            records.append(dict(case=name, status='pass', fixture_sha256=fixture_hash,
                                cells=int(rho.size),
                                inside_cells=int(native['inside'].sum()),
                                max_absolute_errors=maxima))
        print(f'PASS {name}: {rho.size} cells; inside={int(native["inside"].sum())}')

    def rh(shape=(16, 4, 1), mach=3., gamma=1.4, reverse=False):
        # Independent RH forward states; do not call reference forward/inverse helpers.
        r = (gamma + 1.) * mach**2 / ((gamma - 1.) * mach**2 + 2.)
        q = (2. * gamma * mach**2 - (gamma - 1.)) / (gamma + 1.)
        rho, p = np.ones(shape), np.ones(shape)
        vx = np.full(shape, mach * np.sqrt(gamma))
        selection = np.arange(shape[0]) >= shape[0] // 2
        if reverse:
            selection = ~selection
        rho[selection], p[selection] = r, q
        vx[selection] /= r
        if reverse:
            vx = -vx
        return rho, p, [vx, np.zeros(shape), np.zeros(shape)]

    shape = (12, 7, 1)
    one = np.ones(shape)
    run('NUM01_constant', one, one, [one * .3, one * -.7, one * 2.])
    for gamma in (1.4, 5./3.):
        for mach in (1., 1.1, 1.5, 3., 10., 1000., 2000.):
            run(f'NUM10_RH_g{gamma:g}_M{mach:g}', *rh(gamma=gamma, mach=mach),
                gamma=gamma, min_mach=1., physical_mach=mach if mach > 1 else None)
    run('NUM12_min_mach_equal', *rh(gamma=5./3.), gamma=5./3., min_mach=3.,
        analytic=lambda a: np.testing.assert_equal(a['inside'].sum(), 8))
    run('NUM12_min_mach_above', *rh(gamma=5./3.), gamma=5./3.,
        min_mach=np.nextafter(3., np.inf),
        analytic=lambda a: np.testing.assert_array_equal(a['inside'], 0))
    run('NUM07_negative_x', *rh(reverse=True), physical_mach=3.)
    run('NUM08_width2_anisotropic', *rh(), width=2, spacing=(.2, .7, 1.3),
        physical_mach=3.)
    for axis in (1, 2):
        rho, p, vel = rh(shape=(16, 5, 3))
        rho, p = np.swapaxes(rho, 0, axis), np.swapaxes(p, 0, axis)
        v = np.swapaxes(vel[0], 0, axis)
        velocity = [np.zeros_like(v) for _ in range(3)]
        velocity[axis] = v
        run(f'NUM04_RH_axis{axis}', rho, p, velocity, physical_mach=3.)
    shape = (13, 7, 5)
    x, y, z = np.indices(shape)
    rho = 1. + .05 * np.sin(2 * np.pi * x / shape[0])
    p = 1. + .1 * np.cos(2 * np.pi * y / shape[1]) + .02 * z
    run('NUM03_04_noncubic_smooth', rho, p,
        [.01 * x**2, .03 * y, .04 * z], spacing=(.3, .7, 1.1))
    shape = (14, 5, 1)
    x, y, z = np.indices(shape)
    rho = 1. + .03125 * x
    t = 1. + .125 * x
    vel = [-.5 * x, np.zeros(shape), 12. * x + 3. * y]

    def linear_check(result):
        sl = (slice(1, -1), slice(1, -1), slice(None))
        np.testing.assert_allclose(result['grad_tx'][sl], .125, rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(result['grad_rx'][sl], .03125,
                                   rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(result['divergence'][sl], -.5,
                                   rtol=RTOL, atol=ATOL)
        assert not result['grad_tz'].any()
        assert not result['normal_z'].any()

    run('NUM02_04_linear_inactive_z', rho, rho*t, vel, analytic=linear_check)
    run('NUM06_equal_divergence', rho, rho*t, vel, div_limit=-.5,
        analytic=lambda a: np.testing.assert_array_equal(a['candidate'][1:-1], 0))
    run('NUM06_equal_gradient', rho, rho*t, vel, grad_floor=.125,
        analytic=lambda a: np.testing.assert_array_equal(a['candidate'][1:-1], 0))
    rho_orth = 1. + .03125*y
    run('NUM05_orthogonal', rho_orth, rho_orth*t, vel,
        analytic=lambda a: np.testing.assert_array_equal(a['candidate'][1:-1, 1:-1], 0))
    run('NUM05_tiny_gradient', rho, rho*t*1e-30, vel,
        analytic=lambda a: np.testing.assert_array_equal(a['candidate'][1:-1], 0))
    tied_t = 1. + .125*x + .125*y
    run('NUM07_dominant_axis_tie', rho, rho*tied_t, vel, gamma=2.,
        analytic=lambda a: np.testing.assert_allclose(
            a['rho_up'][1:-1, 1:-1], (rho-.03125)[1:-1, 1:-1],
            rtol=RTOL, atol=ATOL))
    derivative_errors = []
    for count in (16, 32, 64):
        phase = 2*np.pi*np.arange(count)[:, None, None]/count
        field = np.broadcast_to(2.+.2*np.sin(phase), (count, 4, 1)).copy()
        spacing_x = 2*np.pi/count

        def sine_check(a, count=count, phase=phase, h=spacing_x):
            exact_stencil = np.broadcast_to(.2*np.cos(phase)*np.sin(h)/h,
                                           (count, 4, 1))
            np.testing.assert_allclose(a['grad_tx'], exact_stencil, rtol=RTOL, atol=ATOL)
            truth = np.broadcast_to(.2*np.cos(phase), (count, 4, 1))
            derivative_errors.append(float(np.max(np.abs(a['grad_tx']-truth))))

        run(f'NUM03_sine_n{count}', np.ones_like(field), field,
            [np.zeros_like(field) for _ in range(3)], spacing=(spacing_x, .7, 1.),
            analytic=sine_check)
    assert derivative_errors[0]/derivative_errors[1] > 3.8
    assert derivative_errors[1]/derivative_errors[2] > 3.8
    for label in ('below', 'equal', 'above'):
        # gamma=2 gives an exactly representable physical saturation limit 3.
        limit = 3.
        r = dict(below=np.nextafter(limit, -np.inf), equal=limit,
                 above=np.nextafter(limit, np.inf))[label]
        rho, p, velocity = rh(gamma=2.)
        rho[8:] = r
        run(f'NUM13_density_bound_{label}', rho, p, velocity, gamma=2.,
            analytic=lambda a, label=label: np.testing.assert_equal(
                a['inside'].sum(), 0 if label == 'above' else 8))
    for floor in (.9, 1., 1.01):
        run(f'NUM14_density_floor_{floor}', *rh(), density_floor=floor)
        run(f'NUM14_temperature_floor_{floor}', *rh(), temperature_floor=floor)
    rho, p, vel = rh()
    run('NUM05_contact', rho, np.ones_like(p), vel,
        analytic=lambda a: np.testing.assert_array_equal(a['inside'], 0))
    # The periodic seam is an independent compression. Only the declared interior
    # is the expansion negative control; full-array port parity is still exact.
    run('NUM05_expansion_interior', rho, p, [-vel[0], vel[1], vel[2]],
        analytic=lambda a: np.testing.assert_array_equal(a['inside'][1:-1], 0))
    shape = (20, 4, 1)
    x, y, z = np.indices(shape)
    r = (1.4+1)*9/((1.4-1)*9+2)
    q = (2*1.4*9-(1.4-1))/(1.4+1)
    level = np.clip((x-6)/6., 0., 1.)
    rho, p = 1+(r-1)*level, 1+(q-1)*level
    vx = 3*np.sqrt(1.4)/rho
    run('NUM17_thick_full_region', rho, p, [vx, np.zeros(shape), np.zeros(shape)],
        width=3, analytic=lambda a: np.testing.assert_array_less(8, a['inside'].sum()))
    rho, p, vel = rh(shape=(20, 4, 1))
    level = ((x >= 4) & (x < 8)) | ((x >= 12) & (x < 16))
    rho, p = np.where(level, r, 1.), np.where(level, q, 1.)
    vx = np.where(level, 3*np.sqrt(1.4)/r, 3*np.sqrt(1.4))
    run('NUM17_two_sheets', rho, p, [vx, np.zeros(shape), np.zeros(shape)],
        analytic=lambda a: np.testing.assert_equal(a['inside'].sum(), 16))
    rng = np.random.default_rng(72941)
    shape = (11, 6, 3)
    rho, p = rng.uniform(.7, 1.4, size=shape), rng.uniform(.6, 2., size=shape)
    run('NUM16_random_seed72941', rho, p,
        [rng.normal(0., .2, size=shape) for _ in range(3)], spacing=(.3, .8, 1.2))
    for value in (0., -1., float('nan'), float('inf')):
        shape = (8, 4, 1)
        rho, p = np.ones(shape), np.ones(shape)
        rho[4, 2, 0] = value
        run(f'NUM15_invalid_rho_{value}', rho, p,
            [np.zeros(shape) for _ in range(3)], invalid=True)
    shape = (8, 4, 1)
    for label, rho, p in [('temperature_overflow', np.full(shape, 1e-200),
                           np.full(shape, 1e200)),
                          ('negative_pressure', np.ones(shape), -np.ones(shape))]:
        run(f'NUM15_{label}', rho, p, [np.zeros(shape) for _ in range(3)],
            invalid=True)
    extreme_path = Path(__file__).resolve().parents[2] / (
        'evidence/shocks/detector/pressure-gradient-overflow.txt')
    lines = extreme_path.read_text().splitlines()
    shape = tuple(map(int, lines[0].split()))
    fields = np.array([[float(x) for x in line.split()]
                       for line in lines[2:]]).reshape(*shape, 5)
    run('NUM15_pressure_gradient_overflow', fields[..., 0], .4*fields[..., 4],
        [fields[..., n] for n in (1, 2, 3)], spacing=(1e-10, 1., 1.),
        energy_override=fields[..., 4], invalid=True, reference_error='grad_px')
    shape = (8, 4, 1)
    rho, p = np.full(shape, 1e150), np.ones(shape)
    p[5, 2, 0] = p[4, 3, 0] = 1.6e308
    run('NUM15_pressure_gradient_magnitude_overflow', rho, p,
        [np.zeros(shape) for _ in range(3)], gamma=2., spacing=(.5, .5, 1.),
        invalid=True, reference_error='grad_p_mag')
    report = dict(status='pass', rtol=RTOL, atol=ATOL, cases=records,
                  reference_sha256=REFERENCE_SHA256,
                  executable_sha256=hashlib.sha256(
                      args.executable.read_bytes()).hexdigest(),
                  sine_derivative_errors=derivative_errors,
                  limitations=['NUM09 unsupported', 'NUM13 consistency-boundary '
                               'not independently injectable through positive EOS states',
                               'NUM19 GPU/single precision not run',
                               'No scientific sign-off'])
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f'{len(records)} cases PASS; exact masks; rtol=atol=1e-12')


if __name__ == '__main__':
    main()
