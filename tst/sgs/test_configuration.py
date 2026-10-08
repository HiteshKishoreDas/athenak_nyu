#!/usr/bin/env python3
"""Exercise persistent configuration results without running a simulation."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

import configuration
from verify import input_text, parameters


def main():
    config = {'mesh/nx1': '8', 'hydro/eos': 'ideal', 'mesh/nx2': '8'}
    assert parameters(input_text(config)) == config
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        exe = root/'athena'
        exe.write_bytes(b'fake executable')
        output = root/'results'

        real_run = subprocess.run

        def run(command, **kwargs):
            if command[0] == 'git':
                return real_run(command, **kwargs)
            rows = json.loads((output/'results.json').read_text())
            assert rows[-1]['status'] == 'running'
            assert rows[-1]['command'] == command
            name = Path(kwargs['cwd']).name
            log = kwargs['stdout']
            if name == 'C05_model':
                log.write('model must be supported\n')
                return subprocess.CompletedProcess(command, 1)
            if name == 'C06_negative' or 'sgs/enabled=false' in command:
                log.write('partial output before timeout\n')
                raise subprocess.TimeoutExpired(command, 10)
            if name == 'L05_zero':
                (Path(kwargs['cwd'])/'rst').mkdir()
                (Path(kwargs['cwd'])/'rst/bad.rst').write_bytes(b'<par_end>\n')
            if '-r' in command:
                log.write('must match checkpoint\n')
                return subprocess.CompletedProcess(command, 1)
            return subprocess.CompletedProcess(command, 0)

        argv = ['configuration.py', '--exe', str(exe), '--output', str(output)]
        with patch.object(sys, 'argv', argv), patch.object(
                configuration.subprocess, 'run', side_effect=run):
            assert configuration.main() == 1
        rows = json.loads((output/'results.json').read_text())
        results = {row['id']: row for row in rows}
        assert len(rows) == 32 and all(row['status'] != 'running' for row in rows)
        assert results['C05_model']['status'] == 'pass'
        for name in ('C06_negative', 'C15_sgs/enabled'):
            assert results[name]['status'] == 'fail'
            assert 'TimeoutExpired' in results[name]['error']
        assert 'error' in results['L05_zero']  # Malformed checkpoint is recorded.
        assert 'error' in results['C02_disabled_vs_omitted']  # Missing states too.
        assert results['C15_sgs/c1']['status'] == 'pass'
        assert results['C15_hydro/nscalars']['status'] == 'pass'
        for log in (output/'C06_negative/run.log', output/'C15_sgs_enabled.log'):
            assert 'partial output before timeout' in log.read_text()
    print('Configuration failure-recording regression passed')


if __name__ == '__main__':
    main()
