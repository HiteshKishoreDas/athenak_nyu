#!/usr/bin/env python3
"""Independent native pressure-work and numerical-heating amplitude oracles."""
import argparse
import json
import math
from pathlib import Path

from verify import Verification, restart, scaled_diff

DT, TOL = 1e-7, 2e-5  # Frozen before execution; initial-rate versus RK2 O(dt).
PRESSURE = 'sgs-pressure-work-integral'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-exe', type=Path)
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    protocol = dict(dt=DT, relative_rate_tolerance=TOL, pressure=1, velocity=.3,
                    density=1, amplitude=.1, gamma=5/3, grids=[16, 32],
                    timeout_seconds=30, no_limiting_or_fallback=True,
                    formulas=dict(compression='3*p0*amp*v*sin(2*pi/N)/(2*dx)',
                                  shear='sqrt(gamma*p0/rho+6*k)*rho*v^2*'
                                        'sin(pi/N)^2/dx'))
    (suite.output/'heating-protocol.json').write_text(json.dumps(protocol, indent=2))
    checks = []
    common = {'mesh/nghost': '3', 'hydro/reconstruct': 'dc',
              'time/tlim': str(DT), 'time/nlim': '2',
              'output1/dt': str(DT), 'output4/dt': str(DT)}
    for n in protocol['grids']:
        for sign in (1, -1):
            initial = suite.output/f'compression_n{n}_sign{sign}.txt'
            with initial.open('w') as handle:
                handle.write(f'{n} {n} {n} 8\n')
                for k in range(n):
                    for j in range(n):
                        for i in range(n):
                            phase = [2*math.pi*(a+.5)/n for a in (i, j, k)]
                            velocity = [-sign*.3*math.sin(a) for a in phase]
                            pressure = 1+.1*sum(math.cos(a) for a in phase)
                            energy = pressure/(5/3-1)+sum(v*v for v in velocity)/2
                            values = [1, *velocity, energy, .2, .25, 0]
                            handle.write(' '.join(format(v, '.17g') for v in values)+'\n')
            name = f'compression_n{n}_sign{sign}'
            row = suite.run(name, 'box', 3, block=n//2, grid_cap=32, timeout=30,
                            overrides={**common, 'problem/initial_state': str(initial),
                                       **{f'mesh/nx{d}': str(n) for d in (1, 2, 3)}})
            if row['status'] != 'pass':
                checks.append(dict(case=name, passed=False))
                continue
            first, last, _, hist = suite.states.pop(name)
            actual = (hist[-1][PRESSURE]-hist[0][PRESSURE])/(last['time']-first['time'])
            expected = sign*3*.1*.3*math.sin(2*math.pi/n)*n/2
            error = abs(actual/expected-1)
            checks.append(dict(case=name, expected_initial_rate=expected, actual=actual,
                               relative_error=error, passed=error <= TOL))
            checks.append(dict(case=name, check='one step', passed=last['cycle'] == 1))
        name = f'shear_n{n}'
        row = suite.run(name, 'shear', 2, block=n//2, grid_cap=32, timeout=30,
                        overrides={**common, 'problem/velocity': '.3',
                                   'mesh/nx1': str(n), 'mesh/nx2': str(n)})
        if row['status'] != 'pass':
            checks.append(dict(case=name, passed=False))
            continue
        first, last, _, hist = suite.states.pop(name)
        pressure = hist[-1][PRESSURE]-hist[0][PRESSURE]
        diss = hist[-1]['sgs-dissipation-integral']-hist[0]['sgs-dissipation-integral']
        actual = (last['thermal']-first['thermal']-pressure-diss)/DT
        expected = math.sqrt(5/3+6*.01)*.3**2*math.sin(math.pi/n)**2*n
        error = abs(actual/expected-1)
        checks.append(dict(case=name, expected_initial_rate=expected, actual=actual,
                           relative_error=error, passed=error <= TOL))
        checks.append(dict(case=name, check='zero pressure work', actual=pressure,
                           passed=abs(pressure) <= 1e-14 and last['cycle'] == 1))
    for row in suite.results:
        checks.append(dict(case=row['name'], check='unlimited',
                           passed=row.get('fallback_stages', 1) == 0
                           and row.get('limited_integral', 1) == 0
                           and row.get('diffusion_limited_integral', 1) == 0))
    temporal = dict(cfl=[.1, .05, .025], time=.02, order_range=[1.5, 2.5],
                    minimum_successive_difference=1e-14,
                    scope='fixed-space coupled evolution and ledger; numerical reference')
    (suite.output/'temporal-protocol.json').write_text(json.dumps(temporal, indent=2))
    values = []
    for cfl in temporal['cfl']:
        name = f'pressure_cfl{cfl}'
        row = suite.run(name, 'box', 3, timeout=30, overrides={
            **common, 'time/tlim': '.02', 'time/nlim': '1000',
            'time/cfl_number': str(cfl), 'output1/dt': '.005', 'output4/dt': '.005',
            'problem/initial_state': str(suite.output/'compression_n16_sign1.txt')})
        if row['status'] != 'pass':
            checks.append(dict(case=name, passed=False))
            continue
        first, last, _, hist = suite.states.pop(name)
        values.append(hist[-1][PRESSURE]-hist[0][PRESSURE])
        checks.append(dict(case=name, check='zero K, no limiting/fallback',
                           passed=last['totals'][-1] == 0
                           and row['fallback_stages'] == row['limited_integral'] == 0
                           and row['diffusion_limited_integral'] == 0))
    if len(values) == 3:
        a, b = values[0]-values[1], values[1]-values[2]
        order = math.log2(abs(a/b)) if b != 0 and a != 0 else None
        checks.append(dict(check='pressure ledger temporal refinement', values=values,
                           successive_differences=[a, b], observed_order=order,
                           passed=order is not None and a*b > 0
                           and min(abs(a), abs(b)) > 1e-14 and 1.5 <= order <= 2.5))
    suite.run('restart_base', 'box', 3, overrides={'mesh/nghost': '3'})
    chosen = next(p for p in sorted((suite.output/'restart_base/rst').glob('*.rst'))
                  if 0 < restart(p)['time'] < .004)
    suite.run('restart_continued', 'box', 3, rst=chosen)
    suite.compare('restart_base', 'restart_continued', 'I13', 1e-12)
    if args.baseline_exe:
        baseline = Verification(args.baseline_exe, suite.output/'baseline')
        baseline.run('box', 'box', 3, overrides={'mesh/nghost': '3'})
        old, new = baseline.states['box'][1], suite.states['restart_base'][1]
        error = max(scaled_diff(x, y) for key in old['field']
                    for x, y in zip(old['field'][key], new['field'][key]))
        checks.append(dict(check='diagnostic field passivity', maximum_error=error,
                           passed=error == 0))
    (suite.output/'heating-results.json').write_text(json.dumps(checks, indent=2))
    passed = all(c['passed'] for c in checks) and all(
        r['status'] == 'pass' for r in suite.results)
    print(f'Heating controls: {sum(c["passed"] for c in checks)}/{len(checks)} checks')
    return int(not passed)


if __name__ == '__main__':
    raise SystemExit(main())
