#!/usr/bin/env python3
"""Small serial-double SGS checks; native-double restart inventories include modal OU state.

Build with -DPROBLEM=sgs_verify. Output directories must be new. Every executed
attempt, including failure/time-limit failure, is kept in results.json. The reader
is deliberately restricted to unrefined serial Newtonian HD native-double files;
it validates metadata, byte counts, cell count and positions before summing.
"""
import argparse
import array
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys
import time

SOURCE = Path(__file__).resolve().parents[2]
LIMITS = json.loads(Path(__file__).with_name('thresholds.json').read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameters(text):
    result = {}
    section = None
    for line in text.splitlines():
        line = line.split('#', 1)[0].strip()
        if line.startswith('<'):
            section = line.strip('<>')
        elif '=' in line and section:
            key, value = line.split('=', 1)
            result[f'{section}/{key.strip()}'] = value.strip()
    return result


def input_text(config):
    blocks = {}
    for key, value in config.items():
        block, name = key.split('/')
        blocks.setdefault(block, []).append(f'{name} = {value}')
    return '\n'.join(f'<{block}>\n'+'\n'.join(values)
                     for block, values in blocks.items())+'\n'


def restart(path):
    raw = path.read_bytes()
    pos = raw.index(b'<par_end>\n') + len(b'<par_end>\n')
    config = parameters(raw[:pos].decode())

    def unpack(fmt):
        nonlocal pos
        values = struct.unpack_from('=' + fmt, raw, pos)
        pos += struct.calcsize('=' + fmt)
        return values

    nmb, root_level = unpack('ii')
    region = unpack('9d')
    mesh = unpack('19i')
    block = unpack('19i')
    t, dt, cycle, ranks = unpack('ddii')
    if (ranks != 1 or config.get('mesh/refinement', 'none') != 'none' or
            config.get('mesh_refinement/refinement', 'none') != 'none'):
        raise ValueError('reader supports serial unrefined restart only')
    locations = [unpack('4i') for _ in range(nmb)]
    if any(loc[3] != root_level for loc in locations):
        raise ValueError('unexpected refinement level')
    pos += 8*nmb + 8*ranks  # per-block float cost/int rank; rank gids/counts
    forced = any(key.startswith('turb_driving/') for key in config)
    rng_state = None
    forcing_metadata = None
    forcing = {}
    if forced:
        # Native C++ record: 20 int fields followed by 18 Real (double) fields.
        # The ints end on an 8-byte boundary, so no internal/tail padding occurs.
        fmt = '@20i18d'
        values = struct.unpack_from(fmt, raw, pos)
        pos += struct.calcsize(fmt)
        version, mode_count, n_updates = values[:3]
        if version != 1 or not 0 < mode_count <= 1000000 or n_updates < 0:
            raise ValueError('unsupported or invalid modal OU metadata')
        forcing_metadata = dict(version=version, mode_count=mode_count,
                                n_updates=n_updates, fields=values)
        rng_size = struct.calcsize('@35qid')
        rng_state = raw[pos:pos+rng_size].hex()
        pos += rng_size
        for part in ('real', 'imag'):
            forcing[part] = unpack(f'{3*mode_count}d')
        if any(not math.isfinite(x) for part in forcing.values() for x in part):
            raise ValueError('nonfinite Fourier OU amplitudes')
    chunk, = unpack('Q')
    ng, nx, ny, nz, il, ih, jl, jh, kl, kh = block[:10]
    shapes = (nx+2*ng, ny+2*ng if ny > 1 else 1, nz+2*ng if nz > 1 else 1)
    cells = math.prod(shapes)
    nscalars = int(config.get('hydro/nscalars', '0'))
    enabled = config.get('sgs/enabled', 'false') in ('true','1')
    nvars = 5+nscalars+int(enabled)
    if chunk != 8*cells*nvars or len(raw) != pos+nmb*chunk:
        raise ValueError('unexpected variable count/layout/file size')
    field = {}
    for loc in locations:
        data = array.array('d')
        data.frombytes(raw[pos:pos+chunk])
        pos += chunk
        if sys.byteorder != 'little':
            raise ValueError('reader currently qualified on little-endian CPU only')
        for k in range(kl, kh+1):
            for j in range(jl, jh+1):
                for i in range(il, ih+1):
                    cell = (k*shapes[1]+j)*shapes[0]+i
                    key = (loc[0]*nx+i-il, loc[1]*ny+j-jl, loc[2]*nz+k-kl)
                    if key in field:
                        raise ValueError('duplicate active cell')
                    field[key] = tuple(data[v*cells+cell] for v in range(nvars))
    if len(field) != math.prod(mesh[1:4]):
        raise ValueError('active-cell count does not match mesh')
    vol = math.prod((region[d+3]-region[d])/mesh[d+1] for d in range(3))
    ordered = [field[key] for key in sorted(field)]
    totals = [math.fsum(row[v] for row in ordered)*vol for v in range(nvars)]
    kinetic = [sum(x*x for x in row[1:4])/(2*row[0]) for row in ordered]
    thermal = [row[4]-ke for row,ke in zip(ordered,kinetic)]
    return dict(time=t, dt=dt, cycle=cycle, config=config, field=field,
                forcing=forcing, forcing_metadata=forcing_metadata, rng_state=rng_state,
                totals=totals, thermal=math.fsum(thermal)*vol, kinetic=math.fsum(kinetic)*vol,
                min_thermal=min(thermal), min_density=min(row[0] for row in ordered),
                min_sgs=min(row[-1] for row in ordered) if enabled else None,
                finite=all(math.isfinite(x) for row in ordered for x in row))


def history(path):
    lines = path.read_text().splitlines()
    names = re.findall(r'\[\d+\]=([^\s]+)', lines[1])
    canonical = ('represented-E', 'sgs-production-integral',
                 'sgs-dissipation-integral', 'sgs-stress-limited-integral',
                 'sgs-fallback-stages', 'sgs-diffusion-limited-integral')
    names = [next((full for full in canonical if full[:10] == name), name)
             for name in names]
    return [dict(zip(names, map(float, line.split())))
            for line in lines if line.strip() and not line.startswith('#')]


def scaled_diff(a, b):
    return abs(a-b)/max(1., abs(a), abs(b))


def state_diff(a, b):
    if a['field'].keys() != b['field'].keys() or len(a['totals']) != len(b['totals']):
        raise ValueError('different global active-cell coordinates')
    return max(scaled_diff(x, y) for key in a['field']
               for x, y in zip(a['field'][key], b['field'][key]))


class Verification:
    def __init__(self, exe, output):
        self.exe = exe.resolve()
        self.output = output.resolve()
        self.output.mkdir(parents=True, exist_ok=False)
        self.results = []
        self.states = {}
        manifest = dict(executable=str(self.exe), executable_sha256=sha(self.exe),
                        source_commit=subprocess.check_output(
                            ['git', 'rev-parse', 'HEAD'], cwd=SOURCE, text=True).strip(),
                        diff_sha256=hashlib.sha256(subprocess.check_output(
                            ['git', 'diff'], cwd=SOURCE)).hexdigest(),
                        thresholds=LIMITS, started_utc=time.strftime('%FT%TZ', time.gmtime()))
        manifest['source_sha256'] = {str(p.relative_to(SOURCE)): sha(p)
                                     for p in (SOURCE/'src').rglob('*') if p.is_file()}
        manifest['harness_sha256'] = {p.name:sha(p) for p in Path(__file__).parent.glob('*')
                                     if p.is_file() and p.suffix in ('.py','.json')}
        (self.output/'manifest.json').write_text(json.dumps(manifest, indent=2))

    def save(self):
        (self.output/'results.json').write_text(json.dumps(self.results, indent=2))

    def check(self, row, test_id, value, threshold, description):
        passed = value <= threshold
        row.setdefault('checks', []).append(dict(id=test_id, passed=passed,
            measured=value, threshold=threshold, description=description))
        return passed

    def run(self, name, mode='uniform', ndim=1, block=8, overrides=None, rst=None,
            grid_cap=None):
        config = parameters((SOURCE/'inputs/hydro/sgs_verify.athinput').read_text())
        config.update({'job/basename': name, 'problem/mode': mode,
                       'time/tlim': '0.004' if ndim == 3 else '0.01', 'time/nlim': str(LIMITS['maximum_cycles']),
                       'output1/data_format': '%24.17e', 'output1/dt': '0.0025',
                       'output2/dt': '-1', 'output3/dt': '-1', 'output4/dt': '0.0025'})
        for d in range(1, 4):
            config[f'mesh/nx{d}'] = '16' if d <= ndim else '1'
            config[f'meshblock/nx{d}'] = str(block) if d <= ndim else '1'
        config.update(overrides or {})
        grid_cap = LIMITS['maximum_grid_axis'] if grid_cap is None else grid_cap
        if any(int(config[f'mesh/nx{d}']) > grid_cap for d in (1, 2, 3)):
            raise ValueError('grid cap exceeded')
        wd = self.output/name
        wd.mkdir()
        inp = wd/'input.athinput'
        inp.write_text(input_text(config))
        cmd = [str(self.exe), '-i', str(inp)]
        if rst:
            # Restart input comes from checkpoint; preserve its integrator/stage state.
            cmd = [str(self.exe), '-r', str(rst), f'job/basename={name}',
                   f"time/tlim={config['time/tlim']}"]
        row = dict(name=name, mode=mode, ndim=ndim, command=cmd, input_sha256=sha(inp),
                   status='running', grid_cap=grid_cap)
        self.results.append(row)
        self.save()
        started = time.monotonic()
        try:
            with (wd/'run.log').open('w') as log:
                proc = subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=LIMITS['run_timeout_seconds'])
            row['returncode'] = proc.returncode
            if proc.returncode:
                raise RuntimeError(f'execution returned {proc.returncode}')
            snapshots = [restart(p) for p in sorted((wd/'rst').glob('*.rst'))]
            first, last = snapshots[0], snapshots[-1]
            row['actual_time'] = last['time']
            row['cycles'] = last['cycle']
            if scaled_diff(last['time'], float(config['time/tlim'])) > 1e-12:
                raise RuntimeError('stopped before requested tlim; nlim is not a pass')
            if any(not snap['finite'] or snap['min_density'] <= 0 or
                   snap['min_thermal'] <= 0 for snap in snapshots):
                raise RuntimeError('nonfinite or nonpositive hydro reservoir')
            if any(snap['min_sgs'] is not None and snap['min_sgs'] < 0 for snap in snapshots):
                raise RuntimeError('negative SGS reservoir')
            row['minima'] = {k: last[k] for k in ('min_density','min_thermal','min_sgs')}
            row['checks'] = [dict(id='R.admissibility', passed=True,
                                 description=LIMITS['positivity'])]
            hist = history(wd/f'{name}.hydro.hst')
            initial_hist, final_hist = hist[0], hist[-1]
            for index, label in enumerate(('mass', '1-mom', '2-mom', '3-mom')):
                self.check(row, 'P22.'+label, scaled_diff(first['totals'][index],
                           last['totals'][index]), LIMITS['mass_momentum_scaled_absolute'],
                           'independent active-cell periodic integral')
                self.check(row, 'P23.'+label, scaled_diff(last['totals'][index],
                           final_hist[label]), LIMITS['history_vs_restart_scaled_absolute'],
                           'native-double restart versus online history')
            if last['min_sgs'] is not None:
                E0, E1 = first['totals'][4], last['totals'][4]
                K0, K1 = first['totals'][-1], last['totals'][-1]
                self.check(row, 'P22.represented', scaled_diff(E0+K0, E1+K1),
                           LIMITS['represented_energy_scaled_absolute'], 'periodic E+K')
                applied = (final_hist['sgs-production-integral']-
                           final_hist['sgs-dissipation-integral'])
                initial_applied = (initial_hist['sgs-production-integral']-
                                   initial_hist['sgs-dissipation-integral'])
                self.check(row, 'P21.paired_transfer', scaled_diff(K1-K0,
                           applied-initial_applied), LIMITS['paired_transfer_scaled_absolute'],
                           'global K change equals integrated applied stage production-dissipation')
                for label, value in [('tot-E',E1),('sgs-E',K1),
                                     ('represented-E',E1+K1),('thermal-E',last['thermal'])]:
                    self.check(row, 'P23.'+label, scaled_diff(value, final_hist[label]),
                               LIMITS['history_vs_restart_scaled_absolute'],
                               'independent active-cell inventory versus history')
                self.check(row, 'P23.kinetic-E', scaled_diff(last['kinetic'],
                           sum(final_hist[f'{v}-KE'] for v in (1,2,3))),
                           LIMITS['history_vs_restart_scaled_absolute'],
                           'independent cellwise kinetic inventory versus history')
                row['applied_production'] = final_hist['sgs-production-integral']
                row['applied_dissipation'] = final_hist['sgs-dissipation-integral']
                row['fallback_stages'] = final_hist['sgs-fallback-stages']
                row['limited_integral'] = final_hist['sgs-stress-limited-integral']
                row['diffusion_limited_integral'] = final_hist['sgs-diffusion-limited-integral']
                for n in range(5, len(last['totals'])-1):
                    self.check(row, f'L02.scalar{n-5}', scaled_diff(first['totals'][n],
                               last['totals'][n]), LIMITS['scalar_mass_scaled_absolute'],
                               'distinct ordinary tracer integral, SGS stored separately')
            if mode == 'uniform' and last['min_sgs'] is not None:
                rho = float(config['problem/density0'])
                k0 = float(config['problem/k0'])
                ceps = float(config['sgs/ceps'])
                dx = 1/int(config['mesh/nx1'])
                expected = rho*k0/(1+ceps*math.sqrt(k0)*last['time']/(2*dx))**2
                row['decay_error'] = abs(last['totals'][-1]-expected)
                row['thermal_oracle_error'] = abs(last['thermal']-(
                    float(config['problem/pressure0'])/(float(config['hydro/gamma'])-1)
                    +rho*k0-expected))
            row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
            self.states[name] = (first, last, snapshots, hist)
        except subprocess.TimeoutExpired:
            row['status'] = 'fail'
            row['error'] = 'pilot exceeded fixed timeout; no completed-solution claim'
        except Exception as exc:
            row['status'] = 'fail'
            row['error'] = str(exc)
        row['seconds'] = time.monotonic()-started
        self.save()
        print(name, row['status'], row.get('error',''), flush=True)
        return row

    def compare(self, left, right, test_id, threshold):
        row = dict(name=f'{test_id}:{left}:{right}', status='running')
        self.results.append(row)
        if left not in self.states or right not in self.states:
            row.update(status='blocked', error='required completed cases unavailable')
        else:
            a, b = self.states[left][1], self.states[right][1]
            self.check(row, test_id, state_diff(a, b), threshold, 'global native-double cell norm')
            if test_id in ('I13','I22'):
                for label in ('sgs-production-integral', 'sgs-dissipation-integral',
                              'sgs-stress-limited-integral', 'sgs-diffusion-limited-integral',
                              'sgs-fallback-stages'):
                    self.check(row, test_id+'.'+label, scaled_diff(
                        self.states[left][3][-1][label], self.states[right][3][-1][label]),
                        threshold, 'restart must retain integrated diagnostics')
            row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
        self.save()


def stress_refinement(suite):
    shear_names = []
    for cfl in ('0.1','0.05','0.025','0.0125'):
        name = 'stress_shear_cfl'+cfl
        shear_names.append(name)
        suite.run(name, mode='shear', ndim=2, overrides={
            'time/cfl_number':cfl, 'time/tlim':'0.04', 'hydro/nscalars':'0',
            'sgs/c2':'0', 'sgs/ceps':'0', 'sgs/ckappa':'0'})
    component = dict(name='P17:isolated_stress_component_time_refinement',
                     status='blocked', description='C2=ceps=ckappa=0 shear. '
                     'Fixed-space full-solver Richardson check; LLF damping remains. '
                     'Finest CFL0.0125 is a numerical reference, not an exact PDE oracle.')
    suite.results.append(component)
    if all(name in suite.states for name in shear_names):
        states = [suite.states[name][1] for name in shear_names]
        errors = [state_diff(states[n], states[n+1]) for n in (0,1,2)]
        component['field_successive_linf'] = errors
        orders = [math.log(errors[n]/errors[n+1],2) for n in (0,1)]
        component['field_orders'] = orders.copy()
        component['reservoirs'] = {}
        for label, values in (
                ('kinetic', [x['totals'][4]-x['thermal'] for x in states]),
                ('thermal', [x['thermal'] for x in states]),
                ('sgs', [x['totals'][-1] for x in states])):
            diffs = [abs(values[n]-values[n+1]) for n in (0,1,2)]
            observed = [math.log(diffs[n]/diffs[n+1],2) for n in (0,1)]
            component['reservoirs'][label] = dict(values=values,
                successive_absolute=diffs, observed_orders=observed)
            orders.extend(observed)
        component['status'] = 'pass' if all(LIMITS['decay_order_min'] <= p <=
            LIMITS['decay_order_max'] for p in orders) else 'fail'
    suite.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--stress-only', action='store_true')
    args = parser.parse_args()
    suite = Verification(args.exe, args.output)
    if args.stress_only:
        stress_refinement(suite)
        return int(any(row['status']=='fail' for row in suite.results))
    for mode, ndim in [('uniform',1),('shear',2),('diffusion',1),('wave',1),
                       ('contact',1),('collision',1),('box',3)]:
        suite.run(f'{mode}_{ndim}d', mode=mode, ndim=ndim)
    if not args.pilot:
        for mode, ndim in [('uniform',3),('shear',3),('diffusion',2),('diffusion',3),
                           ('wave',3),('contact',2),('contact',3),('collision',3),
                           ('box',1),('box',2)]:
            suite.run(f'{mode}_{ndim}d', mode=mode, ndim=ndim)
        suite.run('box_3d_single', mode='box', ndim=3, block=16)
        suite.compare('box_3d','box_3d_single','I04',LIMITS['decomposition_scaled_linf'])
        suite.run('zero_3d', ndim=3, overrides={'problem/k0':'0','problem/velocity':'0'})
        decay_rows = []
        for cfl in ('0.2','0.1','0.05'):
            decay_rows.append(suite.run('decay_cfl'+cfl, overrides={
                'time/cfl_number':cfl,'time/tlim':'0.05'}))
        order = dict(name='P03-P04:decay_temporal_order', status='blocked')
        suite.results.append(order)
        if all(row['status']=='pass' for row in decay_rows):
            errors = [row['decay_error'] for row in decay_rows]
            orders = [math.log(errors[n]/errors[n+1],2) for n in (0,1)]
            order.update(status='pass' if all(LIMITS['decay_order_min'] <= p <=
                LIMITS['decay_order_max'] for p in orders) else 'fail', errors=errors, orders=orders)
        for cfl in ('0.1','0.05','0.025'):
            for mode in ('diffusion','contact'):
                suite.run(mode+'_cfl'+cfl, mode=mode, overrides={'time/cfl_number':cfl})
        for boost in ('0','0.3','100'):
            suite.run('contact_tinyK_boost'+boost, mode='contact', overrides={
                'problem/k0':'1e-12', 'problem/velocity':boost})
        suite.run('highK_lowP_wave', mode='wave', overrides={
            'problem/k0':'0.1', 'problem/pressure0':'0.001'})
        suite.run('hard_collision', mode='collision', overrides={
            'problem/velocity':'5', 'problem/pressure0':'0.01'})
        suite.run('hard_collision_cfl02', mode='collision', overrides={
            'problem/velocity':'20', 'problem/pressure0':'0.0001',
            'problem/k0':'1e-10', 'time/cfl_number':'0.2'})
        limiter = suite.run('backscatter_coarse4_highK', mode='collision', block=4,
            overrides={'mesh/nx1':'4', 'problem/velocity':'20', 'problem/k0':'10',
                       'problem/pressure0':'0.1', 'time/cfl_number':'0.2'})
        if limiter['status'] == 'pass':
            suite.check(limiter, 'R01.active_stress_limiter',
                        int(limiter['limited_integral'] <= 0), 0,
                        'deliberately coarse admissible case must activate stress limiter')
            if not limiter['checks'][-1]['passed']:
                limiter['status'] = 'fail'
        stress_refinement(suite)
        if 'box_3d' in suite.states:
            checkpoints = sorted((suite.output/'box_3d/rst').glob('*.rst'))
            chosen = next((p for p in checkpoints if 0 < restart(p)['time'] < .004), None)
            if chosen:
                suite.run('box_restart', mode='box', ndim=3, rst=chosen)
                suite.compare('box_3d','box_restart','I13',LIMITS['checkpoint_restart_scaled_linf'])
        if 'hard_collision' in suite.states:
            checkpoints = sorted((suite.output/'hard_collision/rst').glob('*.rst'))
            chosen = next((p for p in checkpoints if 0 < restart(p)['time'] < .01
                           and int(restart(p)['config']['sgs/fallbacks']) > 0), None)
            if chosen:
                suite.run('hard_collision_restart', mode='collision', rst=chosen)
                suite.compare('hard_collision','hard_collision_restart','I22',
                              LIMITS['checkpoint_restart_scaled_linf'])
        suite.save()
    print(json.dumps({status:sum(row['status']==status for row in suite.results)
                      for status in ('pass','fail','blocked')},indent=2))
    return int(any(row['status']=='fail' for row in suite.results))


if __name__ == '__main__':
    raise SystemExit(main())
