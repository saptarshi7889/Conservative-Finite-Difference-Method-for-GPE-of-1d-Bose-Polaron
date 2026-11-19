"""Check library archives and compare evolution with the CLI."""
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time

LIBRARY_DRIVER, CLI = map(str, map(Path.resolve, map(Path, sys.argv[1:])))
EXPECTED_TIMES = [0, 0.005, 0.01, 0.015, 0.02]
LIVE_TIMEOUT = 10.0
POLL_INTERVAL = 0.01
EXPECTED_SNAPSHOTS = ['phi_0.dat', 'phi_10.dat', 'phi_20.dat']

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    subprocess.run([LIBRARY_DRIVER, directory], check=True)
    for name, velocity in [('positive', '2'), ('negative', '-2')]:
        library = root / name
        cli = root / (name + '_cli')
        subprocess.run([CLI, '--output', str(cli), '--velocity', velocity,
                        '--velocity-every', '5', '--force', '7'], check=True)
        assert (library / 'velocity.dat').read_bytes() == (cli / 'velocity.dat').read_bytes()
        rows = [line.split() for line in (library / 'velocity.dat').read_text().splitlines()
                if not line.startswith('#')]
        assert [float(row[0]) for row in rows] == EXPECTED_TIMES
        with tarfile.open(library / 'configurations.tar.gz', 'r:gz') as archive:
            assert archive.getnames() == EXPECTED_SNAPSHOTS
            for member in archive:
                assert archive.extractfile(member).read() == (cli / member.name).read_bytes()
        assert not list(library.glob('*.tmp'))
    # Read records while evolution is still running: output must not wait for return.
    process = subprocess.Popen([LIBRARY_DRIVER, directory, 'live'])
    try:
        deadline = time.monotonic() + LIVE_TIMEOUT
        while time.monotonic() < deadline:
            velocity_file = root / 'live' / 'velocity.dat'
            if velocity_file.exists():
                lines = velocity_file.read_text().splitlines()
                if len(lines) >= len(EXPECTED_TIMES) + 1:
                    assert process.poll() is None, 'Simulation finished before streaming check'
                    observed = [float(line.split()[0]) for line in lines[1:len(EXPECTED_TIMES)+1]]
                    assert observed == EXPECTED_TIMES
                    break
            time.sleep(POLL_INTERVAL)
        else:
            raise AssertionError('Velocity records were not visible during evolution')
    finally:
        if process.poll() is None:
            process.terminate()
        process.wait(timeout=LIVE_TIMEOUT)
    print('Library/CLI equivalence, archive integrity, velocity cadence, and restart verified')
