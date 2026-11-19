"""Run the t=10 capped-Picard experiment; update velocity plots and archive outputs."""
from pathlib import Path
import hashlib
import os
import subprocess
import tarfile
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gpe-F7-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results_F7_t10_picard12"
POLL_SECONDS = 10
FINAL_TIME = 10.0
TOLERANCE = 1e-8
FIGSIZE = (7, 4)
DPI = 170
CHUNK_SIZE = 1024 * 1024


def plot_velocity(complete=False):
    path = OUTPUT / "velocity.dat"
    if not path.exists():
        return
    rows = [line.split() for line in path.read_text().splitlines()
            if line and not line.startswith("#")]
    if not rows:
        return
    times = [float(row[0]) for row in rows]
    velocities = [float(row[1]) for row in rows]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(times, velocities)
    ax.set(xlabel="Time", ylabel="Impurity velocity", xlim=(0, FINAL_TIME),
           title="F=7, G=15, M=3m; max 12 Picard iterations")
    ax.grid(alpha=0.25)
    fig.text(0.5, 0.01, "Unconverged steps accepted" +
             ("" if complete else f"; running, saved through t={times[-1]:g}"),
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUTPUT / "velocity.png", dpi=DPI)
    fig.savefig(OUTPUT / "velocity.pdf")
    plt.close(fig)


def main():
    if OUTPUT.exists():
        raise SystemExit(f"Refusing to overwrite {OUTPUT}")
    with (ROOT / "build/F7_t10_solver.log").open("w") as log:
        process = subprocess.Popen([str(ROOT / "build/F7"), str(OUTPUT)],
                                   cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            plot_velocity()
            time.sleep(POLL_SECONDS)
        code = process.wait()
    plot_velocity(complete=code == 0)
    rows = [line.split() for line in (OUTPUT / "convergence.dat").read_text().splitlines()
            if line and not line.startswith("#")]
    violations = sum(float(row[2]) > TOLERANCE for row in rows)
    status = f"exit_code {code}\ntime_steps {len(rows)}\nunconverged_steps {violations}\n"
    (OUTPUT / "status.txt").write_text(status)
    print(status, flush=True)
    archive_path = ROOT / (OUTPUT.name + ".tar.gz")
    with tarfile.open(archive_path, "w:gz") as archive:
        archive.add(OUTPUT, arcname=OUTPUT.name)
    with tarfile.open(archive_path) as archive:
        for member in archive:
            if member.isfile():
                stream = archive.extractfile(member)
                while stream.read(CHUNK_SIZE):
                    pass
    digest = hashlib.sha256()
    with archive_path.open("rb") as stream:
        while chunk := stream.read(CHUNK_SIZE):
            digest.update(chunk)
    archive_path.with_suffix(".sha256").write_text(digest.hexdigest()+"  "+archive_path.name+"\n")
    print(f"Verified archive: {archive_path}", flush=True)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
