#!/usr/bin/env python3
"""Native Eq. 11 amplitude and K/rho diffusion controls at varying density.

Compare diffusion-on/off after one short RK2 step; both retain isotropic SGS
pressure and the same Euler solver. Their difference isolates initial diffusion
up to O(dt), using an independent face flux formula. No turbulence-accuracy claim.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from verify import Verification

DT, TOL = 1e-7, 2e-5  # Frozen initial-rate versus RK2 O(dt) tolerance.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    checks = []
    for n in (16, 32):
        x = (np.arange(n) + .5) / n
        rho = 1 + .2*np.sin(2*np.pi*x)
        for varying_k in (False, True):
            specific = .01*(1 + .15*np.cos(2*np.pi*x)) if varying_k else np.full(n, .01)
            K = rho*specific
            # Independently assemble the published face coefficient and flux.
            rho_face = (rho + np.roll(rho, -1))/2
            K_face = (K + np.roll(K, -1))/2
            J = -.65/n*np.sqrt(rho_face*K_face)*(np.roll(specific, -1)-specific)*n
            expected = -(J-np.roll(J, 1))*n
            if not varying_k:
                assert np.max(abs(expected)) == 0
            else:
                assert np.max(abs(expected)) > 0 and abs(expected.sum()) < 1e-15
            ic = suite.output / f'initial_n{n}_vary{int(varying_k)}.dat'
            data = np.column_stack((rho, np.zeros((n, 3)), np.full(n, 1.5),
                                    .2*rho, .25*rho, K))
            with ic.open('w') as handle:
                handle.write(f'{n} 1 8\n')
                np.savetxt(handle, data, fmt='%.17e')
            states = []
            for ckappa in (0, .65):
                name = f'n{n}_vary{int(varying_k)}_ck{ckappa}'
                row = suite.run(name, 'uniform', 1, block=n//2, grid_cap=32,
                                timeout=30, overrides={
                    'mesh/nx1': str(n), 'mesh/nghost': '3',
                    'problem/initial_state': str(ic), 'hydro/reconstruct': 'ppm4',
                    'sgs/c1': '0', 'sgs/c2': '0', 'sgs/ceps': '0',
                    'sgs/ckappa': str(ckappa), 'time/tlim': str(DT), 'time/nlim': '2',
                    'output1/dt': str(DT), 'output4/dt': str(DT)})
                if row['status'] != 'pass':
                    checks.append(dict(case=name, check='completed native run', passed=False))
                    continue
                first, last, snapshots, hist = suite.states.pop(name)
                checks.append(dict(case=name, check='one unmodified step', passed=(
                    first['cycle'] == 0 and last['cycle'] == 1 and len(snapshots) == 2
                    and abs(last['time']/DT-1) < 1e-12
                    and row['fallback_stages'] == row['limited_integral'] ==
                    row['diffusion_limited_integral'] == 0)))
                states.append(np.array([last['field'][key][-1] for key in sorted(last['field'])]))
            if len(states) == 2:
                actual = (states[1]-states[0])/DT
                error = float(np.max(abs(actual-expected)))
                limit = TOL*float(np.max(abs(expected))) if varying_k else 1e-10
                checks.append(dict(n=n, varying_specific_energy=varying_k,
                    check='published diffusion amplitude' if varying_k else 'zero diffusion at constant K/rho',
                    maximum_absolute_rate_error=error, threshold=limit,
                    expected_rate=expected.tolist(), actual_rate=actual.tolist(), passed=error <= limit))
    (suite.output/'diffusion-results.json').write_text(json.dumps(checks, indent=2)+'\n')
    passed = all(c['passed'] for c in checks) and all(r['status']=='pass' for r in suite.results)
    print(f'Specific-energy diffusion: {sum(c["passed"] for c in checks)}/{len(checks)} checks')
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
