#!/usr/bin/env python3
"""Check signed Fourier support, solenoidal projection and exact restart.

Uses existing SGS forcing checks and native-double checkpoint reader. This is a
forcing-support check, not a stationary turbulence or closure validation.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess

from verify_forcing import ForcingVerification


def modes(lower, upper):
    return [n for n in itertools.product(range(lower, upper+1), repeat=3)
            if 1 <= sum(v*v for v in n) <= 9]


def covariance(vectors):
    cov = [[0.0]*3 for _ in range(3)]
    for n in vectors:
        r2 = sum(v*v for v in n)
        r = math.sqrt(r2)
        weight = max(0.0, 1.0-(r-2.0)**2)*4.0/r2
        for i in range(3):
            for j in range(3):
                cov[i][j] += weight*(int(i == j)-n[i]*n[j]/r2)
    trace = sum(cov[i][i] for i in range(3))
    return [[v/trace for v in row] for row in cov]


def rejected(suite, name, args, diagnostic):
    wd = suite.output/name
    wd.mkdir()
    command = [str(suite.exe)]+args
    with (wd/'run.log').open('w') as log:
        proc = subprocess.run(command, cwd=wd, stdout=log,
                              stderr=subprocess.STDOUT, timeout=20)
    row = dict(name=name, command=command, returncode=proc.returncode,
               expected_diagnostic=diagnostic,
               status='pass' if proc.returncode != 0 and
               diagnostic in (wd/'run.log').read_text() else 'fail')
    suite.results.append(row)
    suite.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = ForcingVerification(args.exe, args.output)
    signed = modes(-3, 3)
    signed_cov, octant_cov = covariance(signed), covariance(modes(0, 3))
    row = dict(name='analytic_signed_covariance', status='running',
               signed=signed_cov, octant=octant_cov, signed_mode_count=len(signed))
    suite.results.append(row)
    suite.check(row, 'signed.isotropic_covariance', max(
        abs(signed_cov[i][j]-int(i == j)/3) for i in range(3) for j in range(3)),
        1e-14, 'ensemble forcing covariance is isotropic on signed mode support')
    suite.check(row, 'signed.octant_anisotropy_oracle',
                abs(octant_cov[0][1]+0.074282817293), 1e-12,
                'independent known octant cross covariance detects missing signs')
    row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
    overrides = {f'turb_driving/min_k{d}':'-3' for d in ('x', 'y', 'z')}
    overrides.update({f'turb_driving/max_k{d}':'3' for d in ('x', 'y', 'z')})
    overrides.update({'turb_driving/npeak':'2', 'output1/dt':'0.025',
                      'output4/dt':'0.1'})
    suite.run_forcing('signed', True, 1.0, tlim=.3, overrides=overrides,
                      retain_snapshots=True)
    if 'signed' in suite.states:
        snapshots = suite.states['signed'][2]
        row = dict(name='native_signed_projection', status='running')
        suite.results.append(row)
        residual = 0.0
        for snap in snapshots:
            count = snap['forcing_metadata']['mode_count']
            suite.check(row, 'signed.native_mode_count', abs(count-len(signed)),
                        0, 'native mode count includes all signed orientations')
            for part in snap['forcing'].values():
                for m, n in enumerate(signed):
                    amp = [part[d*count+m] for d in range(3)]
                    norm = math.sqrt(sum(v*v for v in amp)*sum(v*v for v in n))
                    if norm:
                        residual = max(residual, abs(sum(n[d]*amp[d]
                            for d in range(3)))/norm)
        suite.check(row, 'signed.native_solenoidal_projection', residual, 1e-12,
                    'native OU amplitudes are transverse for every signed mode')
        row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
        suite.forcing_restart('signed_restart', 'signed', True, 1.0, .3)
        input_path = suite.output/'signed/input.athinput'
        checkpoint = sorted((suite.output/'signed/rst').glob('*.rst'))[1]
        diagnostic = 'invalid stage-consistent forcing range or duration'
        for d in ('x', 'y', 'z'):
            for key, value in [('min_k'+d, '-26755'), ('min_k'+d, '-2147483648'),
                               ('max_k'+d, '26755'), ('max_k'+d, '2147483647')]:
                rejected(suite, f'reject_{key}_{value}',
                         ['-i', str(input_path), f'turb_driving/{key}={value}'],
                         diagnostic)
            rejected(suite, 'restart_override_min_k'+d,
                     ['-r', str(checkpoint), f'turb_driving/min_k{d}=-2'],
                     'SGS restart mode, coefficients, forcing and scalar layout '
                     'must match')
    suite.run_forcing('default_octant', True, 1.0, tlim=.03)
    suite.save()
    summary = {status:sum(r['status'] == status for r in suite.results)
               for status in ('pass', 'fail', 'blocked')}
    print(json.dumps(summary, indent=2))
    return int(any(r['status'] != 'pass' for r in suite.results))


if __name__ == '__main__':
    raise SystemExit(main())
