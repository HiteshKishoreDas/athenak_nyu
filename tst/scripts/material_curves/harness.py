"""Check build-tool selection and MPI argv without compiling or launching MPI."""

import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

if __package__:
    from ..utils import athena
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from scripts.utils import athena

_results = []


def run(**kwargs):
    global _results
    _results = []
    original = os.getcwd()
    with tempfile.TemporaryDirectory(prefix='athenak-harness-') as directory:
        directory = str(Path(directory).resolve())
        build = Path(directory)/'build'
        (build/'src').mkdir(parents=True)
        os.chdir(directory)
        try:
            for preferred, available in (
                    ('/tools/cmake3', {'cmake3': '/tools/cmake3', 'cmake': '/tools/cmake'}),
                    ('/tools/cmake', {'cmake': '/tools/cmake'})):
                with patch.object(athena.shutil, 'which', side_effect=available.get), \
                        patch.object(athena.subprocess, 'check_call') as call:
                    athena.make(['-DAthena_ENABLE_MPI=OFF'])
                    commands = [entry.args[0] for entry in call.call_args_list]
                    assert commands[1] == [preferred, '../../', '-DAthena_ENABLE_MPI=OFF']
                    assert os.getcwd() == directory
                    _results.append({'case': 'R05', 'name': preferred.split('/')[-1],
                                     'status': 'pass', 'command': commands[1]})
            with patch.object(athena.shutil, 'which', return_value=None), \
                    patch.object(athena.subprocess, 'check_call'):
                try:
                    athena.make([])
                except athena.AthenaError as error:
                    assert 'Neither cmake3 nor cmake' in str(error)
                else:
                    raise AssertionError('absent build dependency was accepted')
                assert os.getcwd() == directory
                _results.append({'case': 'R05', 'name': 'missing_cmake_rejected',
                                 'status': 'pass'})
            with patch.object(athena.subprocess, 'check_call') as call:
                athena.mpirun(2, 'tests/example.athinput', ['time/nlim=1'])
                command = call.call_args.args[0]
                assert command == ['mpiexec', '-n', '2', './athena', '-i',
                                   '../../../inputs/tests/example.athinput', 'time/nlim=1']
                assert os.getcwd() == directory
                _results.append({'case': 'R05', 'name': 'mpi_separate_argv',
                                 'status': 'pass', 'command': command})
        finally:
            os.chdir(original)
    return 'mocked command checks; no build or MPI process executed'


def analyze():
    return len(_results) == 4 and all(row['status'] == 'pass' for row in _results)


if __name__ == '__main__':
    run()
    print(json.dumps({'schema': 'athenak-material-curves-harness-v1',
                      'status': 'pass' if analyze() else 'fail', 'tests': _results},
                     indent=2, allow_nan=False))
    raise SystemExit(0 if analyze() else 1)
