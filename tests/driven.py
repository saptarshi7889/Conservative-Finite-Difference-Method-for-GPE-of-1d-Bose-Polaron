"""Run driven acceleration and restart checks in an isolated output directory."""
from pathlib import Path
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory() as directory:
    subprocess.run([str(Path(sys.argv[1]).resolve()), directory], check=True)
print('Driven acceleration and compressed/in-memory restarts verified for both force and coupling signs.')
