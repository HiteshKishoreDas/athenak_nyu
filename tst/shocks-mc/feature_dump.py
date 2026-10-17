"""Bounded merge interactions: passivity, disabled paths and explicit unsupported use.

Uses existing line and planar-shock fixtures; never assigns histories to curve nodes.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from common import ROOT, edit_input, provenance, sha
from integration import traces

sys.path.insert(0, str(ROOT/'vis/python'))
from bin_convert import read_binary
import material_curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    binary = args.executable.resolve()
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    report = dict(provenance(binary), runs=[], checks=[], scientific_signoff='pending')

    def case(name, fixture, changes, error=None, restart=None):
        directory = work/name
        directory.mkdir()
        text = edit_input((ROOT/fixture).read_text(), dict({'time/nlim': 3}, **changes))
        path = directory/'input.athinput'
        path.write_text(text)
        command = ([str(binary), '-i', str(path)] if restart is None else
                   [str(binary), '-r', str(restart)] +
                   [key+'='+str(value) for key, value in changes.items()])
        with (directory/'run.log').open('w') as log:
            code = subprocess.run(command, cwd=directory, stdout=log,
                                  stderr=subprocess.STDOUT, timeout=30).returncode
        report['runs'].append(dict(name=name, command=command, exit=code,
                                   expected_error=error,
                                   input_sha256=sha(path) if restart is None else None,
                                   template_input_sha256=sha(path) if restart else None,
                                   log_sha256=sha(directory/'run.log'),
                                   restart_sha256=sha(restart) if restart else None))
        if error is None:
            assert code == 0, name
        else:
            assert code != 0 and error in (directory/'run.log').read_text(), name
            assert not (directory/'shock_events').exists(), name
        return directory

    def fluid(directory):
        result = {}
        for path in sorted((directory/'bin').glob('*.bin')):
            data = read_binary(path)
            fields = {k: np.asarray(v) for k, v in data['mb_data'].items()}
            assert all(np.isfinite(v).all() and v.dtype == np.float64
                       for v in fields.values())
            value = (data['time'], fields)
            cycle = data['cycle']
            if cycle in result:
                assert value[0] == result[cycle][0]
                assert all(np.array_equal(v, result[cycle][1][k]) for k, v in fields.items())
            result[cycle] = value
        assert len(result) == 4  # initial and every one of three steps
        return result

    def same_fluid(left, right):
        a, b = fluid(left), fluid(right)
        assert a.keys() == b.keys()
        for generation, (time, fields) in a.items():
            assert time == b[generation][0]
            assert fields.keys() == b[generation][1].keys()
            assert all(np.array_equal(v, b[generation][1][k]) for k, v in fields.items())

    try:
        line = 'inputs/material_curves/line.athinput'
        outputs = {'output1/file_type': 'bin', 'output1/variable': 'hydro_u',
                   'output1/data_precision': 'real', 'output1/dcycle': 1}
        curves = case('curves', line, outputs)
        off = case('curves_shocks_off', line,
                   dict(outputs, **{'shock_tracking/enabled': 'false'}))
        gas = case('fluid_only', line, dict(outputs, **{
            'material_curves/enabled': 'false', 'shock_tracking/enabled': 'false'}))
        same_fluid(curves, off)
        same_fluid(curves, gas)
        for directory in (curves, off):
            material_curves.read(directory/'material_line.curves.final.txt')
        assert (curves/'material_line.curves.final.txt').read_bytes() == \
            (off/'material_line.curves.final.txt').read_bytes()
        assert not (curves/'shock_events').exists()
        report['checks'].append('float64 fluid passivity and exact curve disabled-shock equivalence')
        checkpoint_output = dict(outputs, **{'output2/file_type': 'rst', 'output2/dt': 1})
        reference = case('curve_restart_reference', line,
                         dict(checkpoint_output, **{'time/nlim': 6}))
        source = case('curve_restart_source', line, checkpoint_output)
        checkpoint = sorted((source/'rst').glob('*.rst'))[-1]
        resumed = case('curve_restart_resumed', line, {'time/nlim': 6}, restart=checkpoint)
        material_curves.read(str(checkpoint)+'.curves', fluid_path=checkpoint)
        assert (reference/'material_line.curves.final.txt').read_bytes() == \
            (resumed/'material_line.curves.final.txt').read_bytes()
        report['checks'].append('native paired curve restart yields exact final snapshot')
        missing = 'shock tracking requires lagrangian_mc particles'
        case('curves_missing_particles', line, {'shock_tracking/enabled': 'true'}, missing)
        case('hydro_missing_particles', line, {
            'material_curves/enabled': 'false', 'shock_tracking/enabled': 'true'}, missing)
        report['checks'].append('orphan shock requests rejected for both curves and hydro')
        shock = 'inputs/tests/shock_mc.athinput'
        enabled = case('shocks', shock, {})
        disabled = case('shocks_curves_off', shock, {'material_curves/enabled': 'false'})
        a, b = traces(enabled/'shock_events'), traces(disabled/'shock_events')
        assert len(a) == len(b) == 512
        assert a.keys() == b.keys()
        assert all(json.dumps(a[k], sort_keys=True) == json.dumps(b[k], sort_keys=True)
                   for k in a)
        report['checks'].append('shock observations unchanged when curves explicitly disabled')
        case('shocks_curves_on', shock, {'material_curves/enabled': 'true'},
             'requires passive Newtonian single-fluid HD without particles/forcing')
        report['checks'].append('curve plus MC shock history remains explicitly unsupported')
        report['status'] = 'passed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        (work/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'report': str(work/'report.json')}))


if __name__ == '__main__':
    main()
