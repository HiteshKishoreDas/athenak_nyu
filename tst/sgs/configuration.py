#!/usr/bin/env python3
"""Check early SGS rejections and disabled-path layout using bounded inputs."""
import argparse
from contextlib import contextmanager
from pathlib import Path
import subprocess

from verify import SOURCE, Verification, input_text, parameters, restart, sha, state_diff


@contextmanager
def attempt(suite, row):
    row['status'] = 'running'
    suite.results.append(row)
    suite.save()
    try:
        yield row
    except Exception as exc:
        row.update(status='fail', error=f'{type(exc).__name__}: {exc}')
    finally:
        suite.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    config = parameters((SOURCE/'inputs/hydro/sgs_verify.athinput').read_text())
    config.update({'time/tlim': '0.001', 'time/nlim': '100'})
    cases = [
        ('C05_model', {'sgs/model': 'unknown'}, 'model must'),
        ('C06_negative', {'sgs/c1': '-1'}, 'coefficients/seed'),
        ('C06_nonfinite', {'sgs/seed': 'nan'}, 'coefficients/seed'),
        ('C07_mixed', {'sgs/c2': '1.1'}, 'c2 must'),
        ('C08_isothermal', {'hydro/eos': 'isothermal',
                           'hydro/iso_sound_speed': '1'}, 'Newtonian ideal'),
        ('C10_amr', {'mesh_refinement/refinement': 'adaptive',
                     'mesh_refinement/max_nmb_per_rank': '64'}, 'AMR/SMR'),
        ('C11_stretched', {'mesh/x2max': '0.6'}, 'cubic cells'),
        ('C12_ghost', {'mesh/nghost': '1'}, 'ghost zones'),
        ('C13_integrator', {'time/integrator': 'rk3'}, 'rk2 integrator'),
        ('C13_solver', {'hydro/rsolver': 'hllc'}, 'llf solver'),
        ('C13_reconstruction', {'hydro/reconstruct': 'wenoz'}, 'plm or dc'),
        ('C16_fofc', {'hydro/fofc': 'true'}, 'unsupported hydro option'),
        ('C16_viscosity', {'hydro/viscosity': '0.1'}, 'unsupported hydro option'),
        ('C16_legacy_cooling', {'hydro/ism_cooling': 'true'}, 'unsupported hydro option'),
        ('C16_cooling', {'cooling/enabled': 'true'}, 'unsupported'),
        ('C16_acceleration', {'hydro_srcterms/const_accel': 'true',
                             'hydro_srcterms/const_accel_val': '0.1',
                             'hydro_srcterms/const_accel_dir': '1'}, 'unsupported'),
        ('C16_frame_tracking', {'frame_tracking/enabled': 'true'}, 'unsupported'),
        ('C16_user_source', {'problem/user_source': 'true'}, 'user source'),
        ('R17_cfl', {'time/cfl_number': '0.3'}, 'cfl_number'),
        ('C16_history_capacity', {'hydro/nscalars': '25'}, 'Hydro history exceeds'),
        ('L05_many_scalars', {'hydro/nscalars': '12'}, None),
        ('L05_zero', {'problem/k0': '0'}, None),
        ('L06_velocity_moments', {'output2/variable': 'hydro_v_moments'}, None),
        ('L06_velocity_density_moments', {'output2/variable': 'hydro_vi_moments'}, None),
        ('L06_sgs_energy_output', {'output2/variable': 'hydro_u_sgs'}, None),
        ('L06_sgs_specific_output', {'output2/variable': 'hydro_w_sgs'}, None),
        ('C02_disabled', {'sgs/enabled': 'false'}, None),
        ('C01_omitted', {}, None),
    ]
    states = {}
    for name, changes, expected in cases:
        run = suite.output/name
        run.mkdir()
        effective = dict(config, **changes)
        if name == 'C01_omitted':
            effective = {k: v for k, v in effective.items()
                         if not k.startswith('sgs/')}
        effective['job/basename'] = name
        inp = run/'input.athinput'
        inp.write_text(input_text(effective))
        command = [str(suite.exe), '-i', str(inp)]
        with attempt(suite, dict(id=name, expected=expected, command=command,
                                 input_sha256=sha(inp))) as row:
            with (run/'run.log').open('w') as log:
                process = subprocess.run(command, cwd=run, stdout=log,
                                         stderr=subprocess.STDOUT, text=True, timeout=10)
            row['returncode'] = process.returncode
            log = (run/'run.log').read_text()
            passed = process.returncode != 0 and expected in log if expected else (
                process.returncode == 0)
            if expected is None and passed:
                snapshots = sorted((run/'rst').glob('*.rst'))
                state = restart(snapshots[-1])
                states[name] = state
                passed &= state['finite'] and state['time'] >= 0.001-1e-14
                if name == 'L05_zero':
                    passed &= all(row[-1] == 0 for row in state['field'].values())
                if name == 'L05_many_scalars':
                    for choice, fields, prefix, extra in (
                            ('hydro_u', ['dens', 'mom1', 'mom2', 'mom3', 'ener'],
                             'r_', 'sgs_energy'),
                            ('hydro_w', ['dens', 'velx', 'vely', 'velz', 'eint'],
                             's_', 'sgs_specific')):
                        expected_fields = fields+[f'{prefix}{n:02d}' for n in range(12)]+[extra]
                        files = sorted((run/'bin').glob(f'*.{choice}.*.bin'))
                        passed &= bool(files)
                        for snapshot in files:
                            with snapshot.open('rb') as output:
                                header = next(line for line in output
                                              if line.startswith(b'  variables:'))
                            passed &= header.decode().split()[1:] == expected_fields
            row['status'] = 'pass' if passed else 'fail'
    with attempt(suite, dict(id='C02_disabled_vs_omitted')) as row:
        difference = state_diff(states['C02_disabled'], states['C01_omitted'])
        row.update(status='pass' if difference == 0 else 'fail', scaled_linf=difference)
    # Mode changes on an SGS checkpoint must fail before array data are loaded.
    for key, value in [('sgs/enabled', 'false'), ('sgs/c1', '0.03'),
                       ('hydro/nscalars', '1')]:
        with attempt(suite, dict(id='C15_'+key)) as row:
            rst = sorted((suite.output/'L05_zero/rst').glob('*.rst'))[-1]
            command = [str(suite.exe), '-r', str(rst), f'{key}={value}']
            row['command'] = command
            suite.save()
            log_path = suite.output/f'C15_{key.replace("/", "_")}.log'
            with log_path.open('w') as log:
                process = subprocess.run(command, cwd=suite.output, stdout=log,
                                         stderr=subprocess.STDOUT, text=True, timeout=10)
            passed = process.returncode != 0 and 'must match' in log_path.read_text()
            row.update(status='pass' if passed else 'fail', returncode=process.returncode)
    print(f'{sum(r["status"] == "pass" for r in suite.results)}/{len(suite.results)} passed')
    return int(any(r['status'] != 'pass' for r in suite.results))


if __name__ == '__main__':
    raise SystemExit(main())
