#!/usr/bin/env python3
"""Native C2/full-mixed production: affine core and exact periodic volume oracle.

Uses the existing test-only initial-state loader; no production instrumentation.
Affine data wrap discontinuously at periodic edges: only the central core has
an analytic local oracle. The smooth periodic fixture checks global exchange.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from verify import Verification

DT = 1e-7
RATE_TOL = 2e-5  # Same frozen O(dt) tolerance as verify_production.py.


def components(a, rho, K, dx, c1=.02, c2=.7):
    """Paper Eq. 8 contracted analytically: shear, nonlinear deviator, compression."""
    strain = (a+a.T)/2
    divergence = np.trace(a)
    norm = np.sum(a*a)
    shear = 2*c1*dx*math.sqrt(2*rho*K)*(
        np.sum(strain*strain)-divergence**2/3)
    nonlinear = -2*c2*K*(np.sum((a@a.T)*strain)/norm-divergence/3)
    compression = -2*K*divergence/3
    return np.array([shear, nonlinear, compression])


def periodic_rate(c2, kg):
    return -2*((1-c2)/3+c2)*kg


def reference_check():
    m = np.array([.3, .2, .1])
    K, c2 = .02, .7
    a = np.outer(m, [1., 0., 0.])
    # Independent rank-one traction includes transverse velocity work.
    tau = -2*c2*K*np.outer(m, m)/np.dot(m, m)+2*c2*K*np.eye(3)/3
    expected = np.dot(tau[:, 0], m)
    assert np.isclose(components(a, 1, K, .1)[1], expected, rtol=1e-13)
    assert np.isclose(periodic_rate(c2, K*m[0])-periodic_rate(0, K*m[0]),
                      expected, rtol=1e-13)
    rotation = np.roll(np.eye(3), 1, axis=0)
    assert np.allclose(components(rotation@a@rotation.T, 1, K, .1),
                       components(a, 1, K, .1), rtol=1e-13, atol=1e-14)
    assert np.allclose(components(a, 3, 3*K, .1),
                       3*components(a, 1, K, .1), rtol=1e-13)
    print('Rank-one transverse work, rotations and density scaling: PASS')


def initial(path, n, a=None):
    x = (np.arange(n)+.5)/n-.5
    z, y, xx = np.meshgrid(x, x, x, indexing='ij')
    if a is None:
        velocity = np.array([.3, .2, .1])[:, None, None, None]*np.sin(2*np.pi*xx)
        rho, K = 1., .02*(1+.4*np.cos(2*np.pi*xx))
    else:
        velocity = np.einsum('ij,jzyx->izyx', a, np.array([xx, y, z]))
        rho, K = (3., .24) if abs(np.trace(a)) > 1e-12 else (1., .01)
    data = np.empty((n, n, n, 8))
    data[..., 0] = rho
    data[..., 1:4] = np.moveaxis(rho*velocity, 0, -1)
    data[..., 4] = 1/(5/3-1)+rho*np.sum(velocity**2, axis=0)/2
    data[..., 5] = .2*rho
    data[..., 6] = .3*rho
    data[..., 7] = K
    with path.open('w') as handle:
        handle.write(f'{n} {n} {n} 8\n')
        np.savetxt(handle, data.reshape(-1, 8), fmt='%.17e')
    return rho, K


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    reference_check()
    if args.self_test:
        return 0
    if args.exe is None or args.output is None:
        parser.error('--exe and --output are required')
    suite = Verification(args.exe, args.output)
    a = np.array([[.3, .2, -.1], [.1, -.2, .15], [-.25, .05, -.1]])
    compression = a-np.eye(3)*.2
    matrices = dict(compression=compression, expansion=-compression,
                    traceless=a, backscatter=-a)
    protocol = dict(dt=DT, rate_relative_tolerance=RATE_TOL,
                    grids=[16, 32], reconstructions=['plm', 'ppm4'],
                    matrices={k:v.tolist() for k, v in matrices.items()},
                    coefficients=dict(c1=.02, c2=.7, ceps=1.5, ckappa=.65),
                    core='abs(x),abs(y),abs(z)<.125; avoids periodic seams',
                    periodic=dict(velocity=[.3, .2, .1], K0=.02, modulation=.4,
                                  c1=0, ceps=0, ckappa=0),
                    scope=__doc__, paper='https://arxiv.org/html/1010.4492')
    (suite.output/'mixed-protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    checks = []
    comparisons = []

    def check(case, label, actual, expected):
        actual = np.asarray(actual)
        scale = max(abs(float(expected)), 1e-12)
        error = float(np.max(np.abs(actual-expected))/scale)
        checks.append(dict(case=case, check=label, actual_min=float(actual.min()),
                           actual_max=float(actual.max()), expected=float(expected),
                           relative_error=error, threshold=RATE_TOL,
                           passed=error <= RATE_TOL))

    for n in protocol['grids']:
        for fixture, matrix in [*matrices.items(), ('periodic', None)]:
            file = suite.output/f'{fixture}{n}.txt'
            rho, K = initial(file, n, matrix)
            for recon in protocol['reconstructions']:
                case = f'{fixture}{n}_{recon}'
                pair = []
                for c2 in (0, .7):
                    name = case+('_base' if c2 == 0 else '_mixed')
                    overrides = {f'mesh/nx{d}':str(n) for d in (1, 2, 3)}
                    overrides.update({'mesh/nghost':'3', 'hydro/reconstruct':recon,
                                      'problem/initial_state':str(file), 'sgs/c2':str(c2),
                                      'time/tlim':str(DT), 'time/nlim':'2',
                                      'output1/dt':str(DT), 'output4/dt':str(DT)})
                    if matrix is None:
                        overrides.update({'sgs/c1':'0', 'sgs/ceps':'0', 'sgs/ckappa':'0'})
                    row = suite.run(name, 'box', 3, block=n//2, overrides=overrides,
                                    grid_cap=32, timeout=30)
                    if row['status'] != 'pass':
                        checks.append(dict(case=name, check='native verification', passed=False))
                        continue
                    first, last, snapshots, hist = suite.states.pop(name)
                    structural = (first['time'] == 0 and last['time'] == DT and
                                  first['cycle'] == 0 and last['cycle'] == 1 and
                                  len(snapshots) == len(hist) == 2)
                    for label in ('sgs-fallback-stages', 'sgs-stress-limited-integral',
                                  'sgs-diffusion-limited-integral'):
                        structural = structural and hist[-1][label] == hist[0][label] == 0
                    checks.append(dict(case=name, check='one unmodified RK2 step',
                                       passed=bool(structural)))
                    pair.append((first, last, hist))
                    if matrix is not None:
                        core = [key for key in first['field']
                                if all(abs((i+.5)/n-.5) < .125 for i in key)]
                        parts = components(matrix, rho, K, 1/n, c2=c2)
                        # Constant K core: conservative advection contributes -K div(u).
                        rate = sum(parts)-K*np.trace(matrix)-1.5*K**1.5/math.sqrt(rho)*n
                        actual = [(last['field'][key][-1]-first['field'][key][-1])/DT
                                  for key in core]
                        check(name, 'full local K rate, including advection and dissipation',
                              actual, rate)
                    else:
                        m = np.array([.3, .2, .1])
                        # Rank-one A=m tensor e_x: all three traction/velocity terms
                        # contract to div(u), not just tau_xx * du_x/dx.
                        # Exact shared-face sum, not a continuum approximation.
                        kg = .02*.4*m[0]*math.cos(math.pi/n)*math.sin(math.pi/n)*n
                        expected = periodic_rate(c2, kg)
                        actual = (hist[-1]['sgs-production-integral']-
                                  hist[0]['sgs-production-integral'])/DT
                        check(name, 'exact periodic applied production', actual, expected)
                        check(name, 'independent periodic K inventory rate',
                              (last['totals'][-1]-first['totals'][-1])/DT, expected)
                if len(pair) != 2:
                    continue
                first, base, _ = pair[0]
                _, mixed, _ = pair[1]
                if matrix is not None:
                    expected = components(matrix, rho, K, 1/n)[1]
                    actual = [(mixed['field'][key][-1]-base['field'][key][-1])/DT
                              for key in core]
                    check(case, 'isolated C2 deviatoric production', actual, expected)
                else:
                    expected = periodic_rate(.7, kg)-periodic_rate(0, kg)
                    check(case, 'isolated C2 periodic K transfer',
                          (mixed['totals'][-1]-base['totals'][-1])/DT, expected)
                    check(case, 'paired resolved energy transfer',
                          (mixed['totals'][4]-base['totals'][4])/DT, -expected)
                    check(case, 'independent resolved kinetic stress work',
                          (mixed['kinetic']-base['kinetic'])/DT, -expected)
                comparisons.append(dict(case=case, expected_C2=float(expected)))
                (suite.output/'mixed-results.json').write_text(json.dumps(checks, indent=2)+'\n')
    summary = dict(passed=sum(c['passed'] for c in checks), total=len(checks),
                   status='pass' if all(c['passed'] for c in checks) else 'fail',
                   maximum_relative_error=max(c.get('relative_error', 0) for c in checks),
                   native_records=len(suite.results), comparisons=comparisons,
                   qualification='Native amplitudes and exchange; no turbulent-reference claim')
    (suite.output/'mixed-results.json').write_text(json.dumps(checks, indent=2)+'\n')
    (suite.output/'mixed-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    return int(summary['status'] != 'pass')


if __name__ == '__main__':
    raise SystemExit(main())
