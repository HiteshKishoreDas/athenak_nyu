#!/usr/bin/env python3
"""Native SGS compression diagnostic against constant entropy K/rho**(5/3).

Uses the existing test-only IC loader and native inventories. This isolates
isotropic SGS pressure: c1=c2=ceps=ckappa=0; it does not retune the LES model.
Periodic shock pairs are a related test, not Semenov's cosmic-ray shock tube.
Scientific entropy errors are reported separately from execution/inventory checks.
"""
import argparse
import json
import math
from pathlib import Path

from verify import Verification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    protocol = dict(grids=[64, 128, 256], cfl=[.2, .1, .05],
                    entropy=1e-6, gamma_sgs=5/3, end_time=.01,
                    shock_states=dict(weak=[1., .1, 5/3],
                                      strong=[17.172, .05, 5/3],
                                      compressible=[17.172, .05, 1.1],
                                      near_isothermal=[17.172, .05, 1.001]),
                    smooth_velocity=.3, timeout_seconds=30,
                    coefficients=dict(c1=0, c2=0, ceps=0, ckappa=0),
                    scope=__doc__)
    (suite.output/'adiabatic-protocol.json').write_text(json.dumps(protocol, indent=2))
    rows = []
    checks = []
    for mode in [*protocol['shock_states'], 'smooth']:
        for n in protocol['grids']:
            pleft, pright, gamma = protocol['shock_states'].get(mode, [1., 1., 5/3])
            initial = suite.output/f'{mode}_n{n}.txt'
            with initial.open('w') as handle:
                handle.write(f'{n} 1 8\n')
                for i in range(n):
                    x = (i+.5)/n-.5
                    rho = 1. if mode == 'smooth' or x < 0 else .125
                    p = pleft if x < 0 else pright
                    v = .3*math.sin(2*math.pi*x) if mode == 'smooth' else 0.
                    K = protocol['entropy']*rho**(5/3)
                    values = [rho, rho*v, 0., 0., p/(gamma-1)+rho*v*v/2,
                              .2*rho, .25*rho, K]
                    handle.write(' '.join(format(v, '.17g') for v in values)+'\n')
            for cfl in protocol['cfl'] if mode != 'smooth' else [.1]:
                name = f'{mode}_n{n}_cfl{cfl}'
                row = suite.run(name, 'contact', 1, block=32, grid_cap=256, timeout=30,
                                overrides={'mesh/nx1': str(n), 'mesh/nghost': '3',
                                           'hydro/reconstruct': 'ppm4',
                                           'hydro/gamma': str(gamma),
                                           'problem/initial_state': str(initial),
                                           'sgs/c1': '0', 'sgs/c2': '0',
                                           'sgs/ceps': '0', 'sgs/ckappa': '0',
                                           'time/cfl_number': str(cfl),
                                           'time/tlim': '.01', 'time/nlim': '20000',
                                           'output1/dt': '.01', 'output4/dt': '.01'})
                if row['status'] != 'pass':
                    continue
                first, last, _, hist = suite.states.pop(name)
                cells = [last['field'][key] for key in sorted(last['field'])]
                rho = [cell[0] for cell in cells]
                K = [cell[-1] for cell in cells]
                expected = [protocol['entropy']*r**(5/3) for r in rho]
                ratio = [k/e for k, e in zip(K, expected)]
                # Uniform passive specific scalar is the same-code entropy-transport control.
                passive_error = max(abs(cell[5]/cell[0]/.2-1) for cell in cells)
                checks.append(dict(case=name, check='passive uniform entropy proxy',
                                   maximum_error=passive_error, passed=passive_error < 1e-12))
                checks.append(dict(case=name, check='initial constant SGS entropy',
                                   passed=max(abs(c[-1]/c[0]**(5/3)/1e-6-1)
                                              for c in first['field'].values()) < 1e-12))
                checks.append(dict(case=name, check='unlimited pressure-only evolution',
                                   passed=row['fallback_stages'] == row['limited_integral']
                                   == row['diffusion_limited_integral']
                                   == row['applied_dissipation'] == 0))
                rows.append(dict(case=name, mode=mode, grid=n, cfl=cfl,
                                 cycles=last['cycle'], mass_weighted_entropy_error=
                                 math.fsum(r*(a-1) for r, a in zip(rho, ratio))/sum(rho),
                                 mass_weighted_entropy_rms=math.sqrt(math.fsum(
                                     r*(a-1)**2 for r, a in zip(rho, ratio))/sum(rho)),
                                 maximum_entropy_excess=max(ratio)-1,
                                 maximum_entropy_deficit=1-min(ratio),
                                 total_K_over_adiabatic=sum(K)/sum(expected),
                                 passive_error=passive_error,
                                 x=[(i+.5)/n-.5 for i in range(n)], density=rho,
                                 entropy_ratio=ratio, sgs_energy=K,
                                 adiabatic_energy=expected))
    decay = []
    for cfl in [.2, .1, .05, .025]:
        name = f'decay_cfl{cfl}'
        row = suite.run(name, 'uniform', 1, timeout=30, overrides={
            'time/cfl_number': str(cfl), 'time/tlim': '.25', 'time/nlim': '20000',
            'output1/dt': '.25', 'output4/dt': '.25', 'problem/velocity': '0'})
        if row['status'] == 'pass':
            first, last, _, _ = suite.states.pop(name)
            exact = .01/(1+1.5*16*math.sqrt(.01)*.25/2)**2
            decay.append(dict(cfl=cfl, cycles=last['cycle'], native_K=last['totals'][-1],
                              analytic_K=exact, relative_error=last['totals'][-1]/exact-1))
    if len(decay) == 4:
        errors = [abs(d['relative_error']) for d in decay]
        ratios = [a/b for a, b in zip(errors[:-1], errors[1:])]
        checks.append(dict(check='analytic decay and temporal refinement',
                           errors=errors, successive_error_ratios=ratios,
                           passed=max(errors) < 1e-4 and all(3.5 < r < 4.5 for r in ratios)))
    result = dict(protocol=protocol, profiles=rows, decay=decay, checks=checks,
                  engineering_pass=all(c['passed'] for c in checks)
                  and len(rows) == 39 and len(decay) == 4
                  and all(r['status'] == 'pass' for r in suite.results))
    (suite.output/'adiabatic-results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"Compression profiles: {len(rows)}; decay controls: {len(decay)}; "
          f"engineering pass: {result['engineering_pass']}")
    return int(not result['engineering_pass'])


if __name__ == '__main__':
    raise SystemExit(main())
