#!/usr/bin/env python3
"""Bounded PPM4 SGS checks, reusing independent native inventory oracles."""
import argparse
from pathlib import Path

import verify
from verify_forcing import ForcingVerification


def run(exe, output):
    suite = ForcingVerification(exe, output)
    common = {'hydro/reconstruct': 'ppm4', 'mesh/nghost': '3'}
    for name, mode, ndim, extra in [
        ('uniform', 'uniform', 3, {}),
        ('zero', 'uniform', 3, {'problem/k0': '0', 'problem/velocity': '0'}),
        ('tiny_contact', 'contact', 1, {'problem/k0': '1e-12'}),
        ('diffusion', 'diffusion', 3, {}),
        ('box', 'box', 3, {}),
        ('collision', 'collision', 1, {'problem/velocity': '5',
                                      'problem/pressure0': '0.01'}),
    ]:
        suite.run(name, mode=mode, ndim=ndim, overrides={**common, **extra})
    suite.run('box_single', mode='box', ndim=3, block=16, overrides=common)
    suite.compare('box', 'box_single', 'PPM.decomposition',
                  verify.LIMITS['decomposition_scaled_linf'])
    paths = sorted((output/'box/rst').glob('*.rst'))
    checkpoint = next(p for p in paths if 0 < verify.restart(p)['time'] < .004)
    suite.run('box_restart', mode='box', ndim=3, overrides=common, rst=checkpoint)
    suite.compare('box', 'box_restart', 'I13',
                  verify.LIMITS['checkpoint_restart_scaled_linf'])
    suite.run('plm_ng2', mode='box', ndim=3)
    suite.run('plm_ng3', mode='box', ndim=3, overrides={'mesh/nghost': '3'})
    suite.compare('plm_ng2', 'plm_ng3', 'PLM.ghost_layout',
                  verify.LIMITS['decomposition_scaled_linf'])
    suite.run_forcing('forced', True, 1., tlim=.05, overrides={
        **common, 'output4/dt': '.01', 'output1/dt': '.01'})
    suite.forcing_restart('forced_restart', 'forced', True, 1., .05)
    suite.save()
    print({status: sum(r['status'] == status for r in suite.results)
           for status in ('pass', 'fail', 'blocked')}, flush=True)
    return all(r['status'] == 'pass' for r in suite.results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(0 if run(args.exe, args.output) else 1)
