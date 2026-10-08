#!/usr/bin/env python3
"""Independent sampled longitudinal mixed-stress acoustic branch oracle.

This is a linear wave with isolated C2 pressure, not a full mixed manufactured
solution. The phase avoids initially exact zero face gradients. Thresholds are
frozen before running; no result establishes general zero-gradient stability.
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
    gamma, k0, amplitude = 5/3, 0.3, 1e-5
    b = 2/3+4*0.7/3
    speed = math.sqrt(gamma+b*(1+b)*k0)
    limits = {'normalized_l1_fine_max': 0.1, 'fine_coarse_ratio_max': 0.4,
              'amplitude': amplitude, 'phase': 0.123, 'b': b,
              'sound_speed': speed, 'time': 0.25/speed,
              'scope': 'linear longitudinal branch, c1=ceps=ckappa=0'}
    (args.output/'wave-thresholds.json').write_text(json.dumps(limits, indent=2)+'\n')
    errors = {}
    checks = []
    for n in (16, 32):
        name = f'mixed_wave_n{n}'
        row = suite.run(
            name, 'wave', 1, block=n, overrides={
                'mesh/nx1': str(n), 'meshblock/nx1': str(n),
                'problem/wave_linear_mixed': 'true',
                'problem/amplitude': str(amplitude), 'problem/k0': str(k0),
                'problem/velocity': '0', 'sgs/c1': '0', 'sgs/ceps': '0',
                'sgs/ckappa': '0', 'time/tlim': str(limits['time']),
                'time/cfl_number': '0.05'})
        if row['status'] != 'pass':
            suite.save()
            return 1
        final = suite.states[name][1]
        scalar_error = max(abs(values[5+s]/values[0]-(0.2+0.05*s))
                           for values in final['field'].values() for s in (0, 1))
        fallbacks = int(final['config']['sgs/fallbacks'])
        checks.append({'grid': n, 'variable': 'passive fractions',
                       'max_error': scalar_error, 'threshold': 1e-12,
                       'passed': scalar_error <= 1e-12})
        checks.append({'grid': n, 'variable': 'fallback stages',
                       'count': fallbacks, 'passed': fallbacks == 0})
        measured = {name: [] for name in ('rho', 'velocity', 'pressure', 'K')}
        for (i, _, _), values in final['field'].items():
            x = (i+0.5)/n
            delta = amplitude*math.sin(2*math.pi*(x-speed*final['time'])+0.123)
            rho, momentum, energy, K = values[0], values[1], values[4], values[-1]
            pressure = (gamma-1)*(energy-sum(v*v for v in values[1:4])/(2*rho))
            expected = {'rho': 1+delta, 'velocity': speed*delta,
                        'pressure': 1+gamma*delta, 'K': k0*(1+(1+b)*delta)}
            actual = {'rho': rho, 'velocity': momentum/rho,
                      'pressure': pressure, 'K': K}
            scales = {'rho': amplitude, 'velocity': speed*amplitude,
                      'pressure': gamma*amplitude, 'K': k0*(1+b)*amplitude}
            for name in measured:
                measured[name].append(abs(actual[name]-expected[name])/scales[name])
        errors[n] = {name: math.fsum(values)/n for name, values in measured.items()}
    for name in errors[16]:
        ratio = errors[32][name]/errors[16][name]
        checks.append({'variable': name, 'coarse_l1': errors[16][name],
                       'fine_l1': errors[32][name], 'ratio': ratio,
                       'passed': errors[32][name] < limits['normalized_l1_fine_max']
                       and ratio < limits['fine_coarse_ratio_max']})
    (args.output/'wave-results.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(checks, indent=2))
    return int(not all(row['passed'] for row in checks))


if __name__ == '__main__':
    raise SystemExit(main())
