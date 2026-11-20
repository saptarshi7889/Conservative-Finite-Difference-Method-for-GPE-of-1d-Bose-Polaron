"""All four library interfaces reject and never record failed Picard steps."""
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

CHECKPOINT_TIME = 0.01
INITIAL_TIME = 0.0

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    subprocess.run([str(Path(sys.argv[1]).resolve()), directory], check=True)
    for case in ('quenched', 'driven'):
        for name, start in [('fresh_failed', INITIAL_TIME), ('restart_failed', CHECKPOINT_TIME)]:
            output = root / case / name
            for filename in ('velocity.dat', 'norm.dat'):
                rows = [line for line in (output / filename).read_text().splitlines()
                        if not line.startswith('#')]
                assert len(rows) == 1 and float(rows[0].split()[0]) == start
            assert len((output / 'convergence.dat').read_text().splitlines()) == 1
            with tarfile.open(output / 'configurations.tar.gz') as archive:
                assert archive.getnames() == ['phi_0.dat']
print('All four interfaces stop before recording failed steps; convergence at the iteration cap succeeds.')
