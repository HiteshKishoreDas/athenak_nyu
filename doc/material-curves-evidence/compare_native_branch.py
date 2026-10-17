"""Compare feature-off native outputs with an untouched multiphase_turb baseline."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--baseline', type=Path, required=True)
parser.add_argument('--port', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--launcher', default='/opt/homebrew/bin/prterun')
a = parser.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
records = {'executables': {'baseline': {'path': str(a.baseline), 'sha256': digest(a.baseline)},
                          'port': {'path': str(a.port), 'sha256': digest(a.port)}},
           'cases': []}
mc = (a.source/'inputs/particles/lagrangian_mc_thermo.athinput').read_text()
plain = mc.split('<particles>')[0] + '<problem>' + mc.split('<problem>')[1].split('<output1>')[0]
plain += '''<output1>
file_type = bin
variable = hydro_u
dt = 0.01
<output2>
file_type = rst
dt = 0.02
'''
configs = [
    ('plain_hydro', plain, 1, []),
    ('standalone_cooling', (a.source/'inputs/unit_tests/cooling_powerlaw.athinput').read_text(), 1, []),
    ('frame_tracking', plain + '''<frame_tracking>
enabled = true
axes = x1
mode = velocity
tau_vel = 0.01
diagnostic_every = 1
''', 1, ['problem/vx0=0.1', 'time/nlim=8', 'time/tlim=0.1']),
    ('forcing', plain + '''<turb_driving>
normalization = edot
dedt = 0.1
nlow = 1
nhigh = 2
tcorr = 0.5
''', 1, []),
    ('lagrangian_mc_serial', mc, 1, []),
    ('lagrangian_mc_mpi', mc, 2, ['meshblock/nx1=16', 'meshblock/nx2=16']),
]
for name, text, ranks, overrides in configs:
    case = a.output/name
    case.mkdir()
    input_file = case/'input.athinput'
    input_file.write_text(text + '\n<material_curves>\nenabled = false\n')
    record = {'name': name, 'ranks': ranks, 'input_sha256': digest(input_file), 'runs': []}
    for label, exe in [('baseline', a.baseline), ('port', a.port)]:
        out = case/label
        out.mkdir()
        command = [a.launcher, '-n', str(ranks), str(exe), '-i', str(input_file), '-d', str(out)] + overrides
        result = subprocess.run(command, cwd=a.source, capture_output=True, text=True, timeout=45)
        (out/'stdout.txt').write_text(result.stdout)
        (out/'stderr.txt').write_text(result.stderr)
        record['runs'].append({'command': command, 'returncode': result.returncode})
        assert result.returncode == 0, (name, label, result.stdout, result.stderr)
        assert 'time =' in result.stdout or 'time=' in result.stdout, (name, label)
        if name == 'frame_tracking':
            assert 'FrameTracker applied boost' in result.stdout, 'frame controller path did not execute'
        if name.startswith('lagrangian_mc'):
            assert list(out.rglob('*.thp')), 'missing thermodynamic tracer history'
        assert not list(out.glob('*.curves.*')), 'disabled curves wrote output'
    def outputs(path):
        return {str(p.relative_to(path)): digest(p) for p in path.rglob('*')
                if p.is_file() and p.name not in ('stdout.txt', 'stderr.txt')}
    before, after = outputs(case/'baseline'), outputs(case/'port')
    assert before, 'no native outputs compared'
    record['files'] = before
    record['identical'] = before == after
    record['mismatches'] = {p: [before.get(p), after.get(p)] for p in before.keys() | after.keys()
                            if before.get(p) != after.get(p)}
    assert record['identical'], (name, record['mismatches'])
    records['cases'].append(record)
    (a.output/'results.json').write_text(json.dumps(records, indent=2)+'\n')
    print(name, 'PASS', len(before), 'byte-identical native output files', flush=True)
