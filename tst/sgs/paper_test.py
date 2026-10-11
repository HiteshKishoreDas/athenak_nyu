#!/usr/bin/env python3
"""Bounded Schmidt--Federrath Section 4 attempt with the native HD solver.

64 cubed is the paper's lowest grid; smaller grids are cost/stability pilots.
The native forcing normalization and solver differ from the published experiment.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np

import verify
from verify_forcing import ForcingVerification

L, V, T = .5, 1., .5
GAMMA, C0 = 1.001, 1 / 5.5
DEPARTURES = [
    'AthenaK RK2/LLF/PLM instead of Enzo 1.5',
    'Instantaneous physical accel_rms normalization instead of modal OU amplitude',
    'Native shell-compensated parabolic envelope; equivalence is unproven',
    'Assumed V/c0=5.5: Table 4 does not explicitly give this for solenoidal forcing',
    'Tiny positive SGS seed k0=1e-8 to avoid the discrete zero-K absorbing state',
    'One seed, no high-resolution filtered reference or convergence claim',
    'Discrete pressure-work ledger; residual combines spatial and RK effects',
]


def cell_diagnostics(data, gamma, delta, ceps=1.5):
    rho = data[:, 0]
    velocity = data[:, 1:4] / rho[:, None]
    thermal = data[:, 4] - .5 * rho * np.sum(velocity**2, axis=1)
    cs2 = gamma * (gamma - 1) * thermal / rho
    k = data[:, -1] / rho
    mach = np.sqrt(np.sum(velocity**2, axis=1) / cs2)
    sgs_mach = np.sqrt(2 * k / cs2)
    eps = ceps * k**1.5 / delta
    fluctuation = velocity - np.mean(velocity, axis=0)
    summary = dict(
        M_rms=float(np.sqrt(np.mean(mach**2))),
        M_sgs=float(np.sqrt(np.mean(sgs_mach**2))),
        K=float(np.mean(data[:, -1])), epsilon=float(np.mean(eps)),
        rho_epsilon=float(np.mean(rho * eps)), rho_mean=float(np.mean(rho)),
        cs2_mean=float(np.mean(cs2)),
        velocity_covariance=(fluctuation.T @ fluctuation / len(rho)).tolist(),
        density_percentiles=np.percentile(rho, [1, 10, 50, 90, 99]).tolist(),
        sgs_mach_percentiles=np.percentile(
            sgs_mach, [1, 10, 50, 90, 99]).tolist(),
    )
    return summary, sgs_mach


def self_check():
    # Two unequal-density cells: volume averages must not become mass averages.
    rho = np.array([1., 3.])
    u = np.array([[2., 0., 0.], [0., 1., 0.]])
    k = np.array([.5, 2.])
    cs2 = np.array([1., 4.])
    data = np.zeros((2, 6))
    data[:, 0] = rho
    data[:, 1:4] = rho[:, None] * u
    data[:, 4] = cs2 * rho / (GAMMA * (GAMMA - 1))
    data[:, 4] += .5 * rho * np.sum(u**2, axis=1)
    data[:, -1] = rho * k
    diag, _ = cell_diagnostics(data, GAMMA, .25)
    assert np.isclose(diag['M_rms'], np.sqrt(2.125), atol=1e-12)
    assert np.isclose(diag['M_sgs'], 1., atol=1e-12)
    assert np.isclose(diag['K'], 3.25)
    assert np.isclose(diag['epsilon'], 3 * (np.sqrt(.125) + np.sqrt(8)))
    assert diag['rho_epsilon'] != diag['epsilon']
    print('Cellwise volume Mach and dissipation diagnostics pass')


def run(args):
    started = time.monotonic()
    suite = ForcingVerification(args.exe, args.output)
    case = 'solenoidal' if args.sol_fraction == 1. else 'compressive'
    protocol = dict(
        L=L, V=V, T=T, gamma=GAMMA, c0=C0, n=args.n,
        amplitude=args.accel_rms, sol_fraction=args.sol_fraction, seed=args.seed,
        end_T=args.end_T, requested_window_T=[2., 10.],
        dt_update_T=.01, snapshot_cadence_T=.1,
        timeout_seconds=args.timeout, cycle_cap=args.nlim,
        reconstruction=args.reconstruct, rsolver="llf", cfl=args.cfl,
        nghost=args.nghost or (3 if args.reconstruct == "ppm4" else 2),
        departures=[DEPARTURES[0].replace("PLM", args.reconstruct.upper()),
                    *DEPARTURES[1:]], coefficients=dict(c1=.02, c2=.7, ceps=1.5,
                                                 ckappa=.65),
        paper='https://arxiv.org/html/1010.4492',
    )
    (suite.output/'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
    row = suite.run_forcing(
        case, True, args.sol_fraction, seed=args.seed,
        tlim=args.end_T*T, timeout=args.timeout,
        retain_snapshots=False,
        overrides={
            **{f'mesh/nx{d}': str(args.n) for d in (1, 2, 3)},
            **{f'meshblock/nx{d}': str(min(args.n, 16)) for d in (1, 2, 3)},
            **{f'turb_driving/min_k{d}': '-3' for d in ('x', 'y', 'z')},
            **{f'turb_driving/max_k{d}': '3' for d in ('x', 'y', 'z')},
            'hydro/reconstruct': args.reconstruct, 'time/cfl_number': str(args.cfl),
            'mesh/nghost': str(protocol['nghost']),
            'hydro/gamma': str(GAMMA), 'problem/pressure0': str(C0**2/GAMMA),
            'problem/k0': '1e-8', 'time/nlim': str(args.nlim), 'time/ndiag': '1000',
            'turb_driving/normalization': 'accel_rms', 'turb_driving/dedt': None,
            'turb_driving/accel_rms': str(args.accel_rms),
            'turb_driving/tcorr': str(T), 'turb_driving/dt_update': str(.01*T),
            'turb_driving/spectrum': 'parabolic', 'turb_driving/npeak': '2',
            'output1/dt': str(.1*T), 'output4/dt': str(.1*T),
        },
    )
    suite.states.clear()
    samples, final = [], {}
    for path in sorted((suite.output/case/'rst').glob('*.rst')):
        state = verify.restart(path)
        data = np.array(list(state['field'].values()))
        sample, sgs_mach = cell_diagnostics(data, GAMMA, 1/args.n)
        sample.update(time=state['time'], kinetic=state['kinetic'],
                      thermal=state['thermal'], represented=sum(
                          state['totals'][v] for v in (4, -1)))
        samples.append(sample)
        del state
    if samples:
        final = dict(density=data[:, 0].tolist(), sgs_mach=sgs_mach.tolist())
    hpath = suite.output/case/f'{case}.hydro.hst'
    diagnostics = dict(protocol=protocol, samples=samples, final=final,
                       history=verify.history(hpath) if hpath.exists() else [],
                       numerical_status=row['status'])
    (suite.output/'diagnostics.json').write_text(
        json.dumps(diagnostics, indent=2) + '\n')
    if samples and diagnostics['history']:
        from paper_plots import write_report
        write_report(suite.output)
    (suite.output/'workflow.json').write_text(json.dumps(dict(
        seconds=time.monotonic()-started, numerical_status=row['status'],
        simulation_and_verification_seconds=row['seconds']), indent=2) + '\n')
    return row['status'] == 'pass'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--nghost', type=int, choices=[2, 3])
    parser.add_argument('--reconstruct', choices=['plm', 'ppm4'], default='plm')
    parser.add_argument('--n', type=int, choices=[16, 32, 64, 128, 256], default=64)
    parser.add_argument('--nlim', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--sol-fraction', type=float, choices=[0., 1.], default=1.)
    parser.add_argument('--end-T', type=float, default=10.)
    parser.add_argument('--cfl', type=float, default=.1)
    parser.add_argument('--accel-rms', type=float, default=2.)
    parser.add_argument('--timeout', type=float, default=7200.)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_check()
        return
    if args.output is None:
        parser.error('--output is required')
    if args.plot_only:
        from paper_plots import write_report
        write_report(args.output)
        return
    if args.exe is None:
        parser.error('--exe is required')
    if (not 0 < args.end_T <= 10 or args.timeout <= 0 or args.accel_rms <= 0
            or not 0 < args.cfl <= .2 or args.nlim <= 0):
        parser.error('require 0<end-T<=10, positive timeout/accel-rms/nlim, 0<cfl<=.2')
    if args.reconstruct == 'ppm4' and args.nghost == 2:
        parser.error('ppm4 requires three ghost zones')
    raise SystemExit(0 if run(args) else 1)


if __name__ == '__main__':
    main()
