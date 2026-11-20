"""Compile small temporary variants of all four real legacy programs.

Only embedded run settings and checkpoint input location change in the copies;
the numerical update, convergence test and failure guard run unmodified.
"""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
COMPILER = sys.argv[1]
EIGEN_INCLUDE = sys.argv[2]
SETTINGS = {'m_tilde': '0.125', 'G_tilde': '-1.5', 'gamma': '0.1',
            'delta_t': '0.001', 'delta_x': '0.125', 'Length': '16',
            'K_fin': '1', 'err_iter': 'std::stod(std::getenv("PICARD_TEST_TOLERANCE"))',
            'err_config': '0.0'}
GRID_POINTS = 128
LENGTH = 16
DX = 0.125
FAILURE_TOLERANCE = '1e-30'
SUCCESS_TOLERANCE = '1e6'

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    for index, original in enumerate(sorted((ROOT / 'legacy').glob('*.cpp'))):
        source = original.read_text()
        continued = '_contd_' in original.name
        start = '0.01' if continued else '0.0'
        for name, value in {**SETTINGS, 't_init': start, 't_fin': start}.items():
            source, count = re.subn(r'\b' + name + r'\s*=\s*[^,;\n]+',
                                    f'{name}={value}', source, count=1)
            assert count == 1, (original, name)
        if continued:
            source, count = re.subn(r'input_Phi.open\([^\n]+;',
                                    'input_Phi.open("checkpoint.dat");', source, count=1)
            assert count == 1
        variant = root / f'legacy_{index}.cpp'
        variant.write_text(source)
        executable = root / f'legacy_{index}'
        subprocess.run([COMPILER, '-std=c++17', '-O0',
                        *['-I' + path for path in EIGEN_INCLUDE.split(';')],
                        str(variant), '-o', str(executable)], check=True)
        for mode, tolerance in [('failure', FAILURE_TOLERANCE), ('success', SUCCESS_TOLERANCE)]:
            output = root / f'{index}_{mode}'
            output.mkdir()
            if continued:
                (output / 'checkpoint.dat').write_text(''.join(
                    f'{i * DX - LENGTH / 2} 1 0\n' for i in range(GRID_POINTS)))
            result = subprocess.run([str(executable)], cwd=output, capture_output=True,
                                    text=True, env={**os.environ, 'PICARD_TEST_TOLERANCE': tolerance})
            if mode == 'failure':
                assert result.returncode != 0, original
                assert 'after 1 iterations; residual ' in result.stderr, result.stderr
                assert ', tolerance ' in result.stderr
                assert not list(output.glob('Phi_solution*')), 'Failed field was saved'
                velocity = next(output.glob('impurity_velocity*'))
                rows = velocity.read_text().splitlines()
                assert len(rows) == 1 and float(rows[0].split()[0]) == float(start)
            else:
                assert result.returncode == 0, result.stdout + result.stderr
                assert 'Nonlinear iteration failed' not in result.stderr
        print(original.name + ': failure stops before output; convergence on final iteration succeeds')
