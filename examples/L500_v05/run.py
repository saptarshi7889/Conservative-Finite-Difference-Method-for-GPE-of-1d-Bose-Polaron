"""Run only to t=10, compress completed configurations, and make review plots."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results_L500_v05"
EXAMPLE = Path(__file__).resolve().parent
POLL_SECONDS = 0.25
COMPRESSION_WORKERS = 1
REVIEW_TIME = 10
COMMAND = [str(ROOT / "build/gpe"), "--length", "500", "--dx", "0.001",
           "--velocity", "0.5", "--end", str(REVIEW_TIME), "--dt", "0.01",
           "--mass-ratio", "0.333", "--coupling", "1", "--gamma", "0.1",
           "--iterations", "10", "--tolerance", "1e-8", "--save-every", "50",
           "--velocity-every", "10", "--output", str(OUTPUT)]


def compress(path):
    # Only final .dat snapshots are submitted; the solver writes .dat.tmp first.
    subprocess.run(["xz", "-9e", "--threads=1", "--check=sha256", "--keep", str(path)], check=True)
    subprocess.run(["xz", "--test", str(path) + ".xz"], check=True)
    path.unlink()
    print(f"Compressed and verified {path.name}; raw file removed", flush=True)


def main():
    with (ROOT / "build/L500_v05.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if OUTPUT.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {OUTPUT}")
        pending = {}
        with (EXAMPLE / "solver.log").open("w") as log, ThreadPoolExecutor(max_workers=COMPRESSION_WORKERS) as worker:
            process = subprocess.Popen(COMMAND, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            try:
                while True:
                    if OUTPUT.exists():
                        for path in sorted(OUTPUT.glob("phi_*.dat")):
                            if path not in pending:
                                pending[path] = worker.submit(compress, path)
                    try:
                        code = process.wait(timeout=POLL_SECONDS)
                        break
                    except subprocess.TimeoutExpired:
                        pass
                # Capture a final rename that happened between the last scan and exit.
                for path in sorted(OUTPUT.glob("phi_*.dat")):
                    if path not in pending:
                        pending[path] = worker.submit(compress, path)
            finally:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
            errors = []
            for path, future in pending.items():
                try:
                    future.result()
                except Exception as error:
                    errors.append(f"{path.name}: {error}")
                    print(f"Compression failed; raw data retained: {path}: {error}", flush=True)
        if code:
            raise SystemExit(f"Solver failed (exit {code}); partial output retained. See solver.log.")
        env = dict(os.environ, MPLCONFIGDIR="/tmp/gpe-L500-mpl")
        subprocess.run([sys.executable, str(EXAMPLE / "plot.py")], env=env, check=True)
        if errors:
            raise SystemExit("Compression errors: " + "; ".join(errors))
        print("Stopped at t=10. No further evolution has been launched.", flush=True)


if __name__ == "__main__":
    main()
