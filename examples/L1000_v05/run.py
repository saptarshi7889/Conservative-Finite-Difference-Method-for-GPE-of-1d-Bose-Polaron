"""Run the local test, compress closed output with XZ -9e, then plot it (Linux)."""
import concurrent.futures
import ctypes
import os
from pathlib import Path
import select
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results_L1000_v05"
ARCHIVES = (ROOT / "results_L1000_v05_previous_s10",
            ROOT / "results_L1000_v05_interrupted_s50",
            Path("/tmp/bose-polaron-L1000-pilot"),
            Path("/tmp/bose-polaron-L1000-pilot-50"))
IN_CLOSE_WRITE = 0x00000008
IN_Q_OVERFLOW = 0x00004000
EVENT = struct.Struct("iIII")
READ_BYTES = 65536
POLL_SECONDS = 1
COMPRESSION_WORKERS = 1
END_TIME = 25
COMMAND = [str(ROOT / "build/gpe"), "--length", "1000", "--dx", "0.001",
           "--dt", "0.01", "--end", str(END_TIME), "--mass-ratio", "0.333",
           "--coupling", "1", "--gamma", "0.1", "--velocity", "0.5",
           "--tolerance", "1e-8", "--iterations", "10", "--save-every", "50",
           "--output", str(OUTPUT)]


def compress(path):
    if path.stat().st_size:
        compressed = path.with_suffix(path.suffix + ".xz")
        if compressed.exists():
            # An interrupted compressor can leave a partial archive. The raw
            # source is still intact, so discard only the invalid replacement.
            check = subprocess.run(["xz", "-t", str(compressed)], capture_output=True)
            if check.returncode:
                compressed.unlink()
            else:
                raise RuntimeError(f"Both raw and valid compressed files exist: {path}")
        subprocess.run(["xz", "-9e", "--threads=1", "--check=sha256", str(path)],
                       check=True)
        print(f"Compressed {path.name}", flush=True)


def main():
    # Watch the parent before the solver creates its directory. The solver must
    # create OUTPUT itself; attach its watch as soon as it exists, before output
    # records close. A final scan also covers any close event before attachment.
    libc = ctypes.CDLL(None, use_errno=True)
    fd = libc.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
    if fd < 0:
        raise OSError(ctypes.get_errno(), "inotify_init1")
    pending, seen = [], set()
    with concurrent.futures.ThreadPoolExecutor(max_workers=COMPRESSION_WORKERS) as pool, \
         concurrent.futures.ThreadPoolExecutor(max_workers=COMPRESSION_WORKERS) as archive_pool:
        def queue(path, archive=False):
            if path not in seen and path.suffix in (".dat", ".txt"):
                seen.add(path)
                pending.append((archive_pool if archive else pool).submit(compress, path))

        for archive in ARCHIVES:
            if archive.exists():
                for path in sorted(archive.iterdir()):
                    queue(path, archive=True)
        process = subprocess.Popen(COMMAND, cwd=ROOT)
        while not OUTPUT.exists() and process.poll() is None:
            time.sleep(POLL_SECONDS)
        if OUTPUT.exists():
            watch = libc.inotify_add_watch(fd, os.fsencode(OUTPUT), IN_CLOSE_WRITE)
            if watch < 0:
                process.terminate()
                process.wait()
                raise OSError(ctypes.get_errno(), "inotify_add_watch")
            # Recover files already closed before watch attachment. Never touch
            # a file still open by the solver.
            opened = set()
            for descriptor in Path(f"/proc/{process.pid}/fd").iterdir():
                try:
                    opened.add(Path(os.readlink(descriptor)))
                except FileNotFoundError:
                    pass
            for path in OUTPUT.iterdir():
                if path not in opened:
                    queue(path)
        while process.poll() is None:
            if select.select([fd], [], [], POLL_SECONDS)[0]:
                events = os.read(fd, READ_BYTES)
                offset = 0
                while offset < len(events):
                    _, mask, _, length = EVENT.unpack_from(events, offset)
                    name = os.fsdecode(events[offset + EVENT.size:offset + EVENT.size + length].split(b"\0")[0])
                    offset += EVENT.size + length
                    if mask & IN_Q_OVERFLOW:
                        print("Watch queue overflow; final scan will recover remaining files.", flush=True)
                    if mask & IN_CLOSE_WRITE:
                        queue(OUTPUT / name)
            for future in pending:
                if future.done():
                    future.result()
        code = process.wait()
        os.close(fd)
        if OUTPUT.exists():
            for path in sorted(OUTPUT.iterdir()):
                queue(path)
        for future in pending:
            future.result()
    if code:
        raise SystemExit(f"Solver failed (exit {code}); compressed partial output preserved.")
    env = dict(os.environ, MPLCONFIGDIR="/tmp/bose-polaron-mpl")
    result = subprocess.run([sys.executable, str(Path(__file__).with_name("plot.py")), str(OUTPUT)],
                            env=env, capture_output=True, text=True)
    Path(__file__).with_name("plot.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise SystemExit(result.returncode)
    report = Path(__file__).with_name("README.md")
    report.write_text(report.read_text().replace("Status: simulation running.",
        "Status: review simulation completed to t=25; awaiting approval to continue.\n\n```text\n" + result.stdout + "```\n\n"
        "![Velocity](velocity.png)\n\n![Density evolution](density_evolution.png)\n\n"
        "![Density profiles](density_profiles.png)\n\n![Synchronized animation](evolution.gif)"))
    print(result.stdout, flush=True)


if __name__ == "__main__":
    main()
