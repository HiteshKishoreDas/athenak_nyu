#!/usr/bin/env python3
"""Native C1 production and dissipation amplitude oracle for cyclic 3D shear.

Uniform rho and k make initial diffusion zero. The cyclic velocity gradient
has diagonal A A^T and zero divergence, so C2/isotropic initial production is
zero. This checks applied face/stage amplitudes, not general mixed production.
"""
import argparse
import json
import math
from pathlib import Path

from verify import Verification


DT = 1e-7
RATE_TOL = 2e-5  # Frozen before execution; RK2 rates differ from t=0 by O(dt).


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    protocol = dict(dt=DT, rate_relative_tolerance=RATE_TOL, velocity=0.3,
                    c1=0.02, c2=0.7, ceps=1.5, ckappa=0.65,
                    grids=[16, 32], density_k_pairs=[[1, 0.01], [3, 0.08]],
                    reconstructions=['plm', 'ppm4'], timeout_seconds=30,
                    scope=__doc__)
    (suite.output/'production-protocol.json').write_text(
        json.dumps(protocol, indent=2)+'\n')
    checks = []
    for n in protocol['grids']:
        for rho, k in protocol['density_k_pairs']:
            dx = 1/n
            s = 2*protocol['velocity']*math.sin(math.pi/n)/dx
            expected_p = 3*protocol['c1']*dx*math.sqrt(2*rho*rho*k)*s*s/2
            expected_d = protocol['ceps']*rho*k**1.5/dx
            for reconstruction in protocol['reconstructions']:
                name = f'production_n{n}_rho{rho}_k{k}_{reconstruction}'
                overrides = {f'mesh/nx{axis}': str(n) for axis in (1, 2, 3)}
                overrides.update({
                    'mesh/nghost': '3', 'hydro/reconstruct': reconstruction,
                    'problem/amplitude': '0', 'problem/density0': str(rho),
                    'problem/k0': str(k), 'problem/velocity': '0.3',
                    'time/tlim': str(DT), 'time/nlim': '2',
                    'output1/dt': str(DT), 'output4/dt': str(DT)})
                row = suite.run(name, 'box', 3, block=n//2, overrides=overrides,
                                grid_cap=32, timeout=30)
                if row['status'] != 'pass':
                    checks.append(dict(case=name, passed=False, reason=row.get('error')))
                    continue
                first, last, snapshots, hist = suite.states.pop(name)
                elapsed = last['time']-first['time']
                structural = (first['time'] == 0 and abs(elapsed/DT-1) < 1e-12
                              and first['cycle'] == 0 and last['cycle'] == 1
                              and len(snapshots) == 2 and len(hist) == 2
                              and abs(first['totals'][-1]/(rho*k)-1) < 1e-12)
                for label in ('sgs-fallback-stages', 'sgs-stress-limited-integral',
                              'sgs-diffusion-limited-integral'):
                    structural = structural and hist[-1][label] == hist[0][label] == 0
                checks.append(dict(case=name, check='one unmodified RK2 step',
                                   passed=structural))
                for label, expected in (('production', expected_p),
                                        ('dissipation', expected_d)):
                    key = f'sgs-{label}-integral'
                    actual = (hist[-1][key]-hist[0][key])/elapsed
                    error = abs(actual/expected-1)
                    checks.append(dict(case=name, check=label, actual=actual,
                                       expected_initial_rate=expected,
                                       relative_error=error, threshold=RATE_TOL,
                                       passed=error <= RATE_TOL))
    (suite.output/'production-results.json').write_text(
        json.dumps(checks, indent=2)+'\n')
    passed = sum(check['passed'] for check in checks)
    print(f'Native production oracle: {passed}/{len(checks)} checks passed')
    return int(not all(check['passed'] for check in checks))


if __name__ == '__main__':
    raise SystemExit(main())
