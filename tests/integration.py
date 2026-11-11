"""Small numerical and CLI regression checks; Python standard library only."""
from pathlib import Path
import math
import subprocess
import sys
import tempfile

PROGRAM = str(Path(sys.argv[1]).resolve())
ABS_TOL = 1e-12
RESTART_TOL = 2e-13
NORM_TOL = 1e-9
CHECKPOINT_STEP = 10
GRID_POINTS = 128


def rows(path):
    return [[float(value) for value in line.split()]
            for line in path.read_text().splitlines() if not line.startswith('#')]


def run(root, name, *options, success=True):
    result = subprocess.run([PROGRAM, '--output', str(root / name), *options],
                            capture_output=True, text=True)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return root / name


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    for label, velocity in [('positive', '2'), ('negative', '-2')]:
        full = run(root, label, '--velocity', velocity)
        half = run(root, label + '_half', '--velocity', velocity, '--end', '0.01')
        restart_v = rows(half / 'velocity.dat')[-1][1]
        resumed = run(root, label + '_restart', '--velocity', velocity,
                      '--restart', str(half / f'phi_{CHECKPOINT_STEP}.dat'),
                      '--start', '0.01', '--end', '0.02',
                      '--restart-velocity', repr(restart_v))
        expected = rows(full / 'phi_20.dat')
        actual = rows(resumed / 'phi_10.dat')
        assert len(actual) == GRID_POINTS
        assert max(abs(a-b) for left, right in zip(expected, actual)
                   for a, b in zip(left, right)) < RESTART_TOL
        for filename in ['velocity.dat', 'norm.dat', 'convergence.dat']:
            assert all(math.isfinite(v) for row in rows(full / filename) for v in row)
        print(label, 'restart matches continuous run')
    # V=0 and G=0 has a spatially uniform solution of unit modulus.
    uniform = run(root, 'uniform', '--velocity', '0', '--coupling', '0')
    assert max(abs(re*re + im*im - 1) for _, re, im in rows(uniform / 'phi_20.dat')) < NORM_TOL
    norms = rows(uniform / 'norm.dat')
    assert max(abs(row[2] - 16) for row in norms) < NORM_TOL
    print('uniform solution and particle number preserved')
    # Report interacting norm drift without asserting an unestablished bound.
    data = rows(root / 'positive' / 'norm.dat')
    print('interacting relative particle-number drift:',
          max(abs(row[2] / data[0][2] - 1) for row in data))
    run(root, 'positive', success=False)  # Never append to existing results.
    run(root, 'invalid', '--dx', '0', success=False)
    run(root, 'fractional', '--end', '0.0205', success=False)
    run(root, 'unknown', '--unknown', '1', success=False)
    run(root, 'missing', '--restart', 'missing.dat', '--start', '0.01',
        '--restart-velocity', '2', success=False)
    malformed = root / 'bad.dat'
    malformed.write_text('0 1 0\n')
    run(root, 'bad', '--restart', str(malformed), '--start', '0.01',
        '--restart-velocity', '2', success=False)
    run(root, 'nonconvergent', '--iterations', '1', success=False)
    print('invalid inputs and nonlinear nonconvergence rejected')
