"""Plot the solver velocity log and the t=10 density checkpoint."""
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results_L500_v05"
OUTPUT = Path(__file__).resolve().parent
REVIEW_TIME = 10.0
SNAPSHOT_STEP = 1000
GRID_POINTS = 500000
LENGTH = 500.0
SPACING = 0.001
GRID_TOLERANCE = 1e-9
X_LIMIT = 20.0
FIGURE_SIZE = (8, 4.5)
FIGURE_DPI = 160
EXPECTED_VELOCITY_ROWS = 101
EXPECTED_SNAPSHOTS = 21
MAX_ITERATIONS = 10
TOLERANCE = 1e-8

velocity = np.loadtxt(RESULTS / "velocity.dat", ndmin=2)
assert len(velocity) == EXPECTED_VELOCITY_ROWS
assert velocity[-1, 0] == REVIEW_TIME
assert np.isfinite(velocity).all()
assert np.allclose(np.diff(velocity[:, 0]), 0.1, rtol=0, atol=GRID_TOLERANCE)
with lzma.open(RESULTS / f"phi_{SNAPSHOT_STEP}.dat.xz", "rt") as stream:
    field = np.loadtxt(stream)
assert field.shape == (GRID_POINTS, 3)
assert np.isfinite(field).all()
assert np.max(abs(field[:, 0] - (np.arange(GRID_POINTS)*SPACING-LENGTH/2))) < GRID_TOLERANCE
selected = abs(field[:, 0]) < X_LIMIT
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
ax.plot(velocity[:, 0], velocity[:, 1])
ax.set(xlabel="Time", ylabel="Impurity velocity", title="L=500: velocity from velocity.dat")
fig.tight_layout()
fig.savefig(OUTPUT / "velocity_t10.png", dpi=FIGURE_DPI)
plt.close(fig)
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
ax.plot(field[selected, 0], field[selected, 1]**2 + field[selected, 2]**2)
ax.set(xlabel="Position in impurity frame", ylabel=r"$|\psi(x)|^2$",
       title="L=500: boson density at t=10", xlim=(-X_LIMIT, X_LIMIT))
fig.tight_layout()
fig.savefig(OUTPUT / "density_t10.png", dpi=FIGURE_DPI)
plt.close(fig)
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
ax.plot(field[:, 0], field[:, 1]**2 + field[:, 2]**2)
ax.set(xlabel="Position in impurity frame", ylabel=r"$|\psi(x)|^2$",
       title="L=500: full-system boson density at t=10", xlim=(-LENGTH/2, LENGTH/2))
fig.tight_layout()
fig.savefig(OUTPUT / "density_t10_full_system.png", dpi=FIGURE_DPI)
plt.close(fig)
convergence = np.loadtxt(RESULTS / "convergence.dat", ndmin=2)
norm = np.loadtxt(RESULTS / "norm.dat", ndmin=2)
assert convergence[-1, 0] == REVIEW_TIME
assert convergence[:, 1].max() <= MAX_ITERATIONS
assert convergence[:, 2].max() <= TOLERANCE
assert len(list(RESULTS.glob("phi_*.dat.xz"))) == EXPECTED_SNAPSHOTS
assert not list(RESULTS.glob("phi_*.dat"))
assert not list(RESULTS.glob("*.tmp"))
summary = (f"Stopped at t={REVIEW_TIME:g}; no further evolution launched.\n"
           f"Final velocity: {velocity[-1, 1]:.12g}\n"
           f"Velocity records: {len(velocity)} (every 0.1 time units, including t=0)\n"
           f"Compressed configurations: {EXPECTED_SNAPSHOTS} (every 0.5, including t=0)\n"
           f"Max Picard iterations: {int(convergence[:, 1].max())}\n"
           f"Max accepted iterate difference: {convergence[:, 2].max():.12g}\n"
           f"Max sampled relative particle-number drift: {np.max(abs(norm[:, 2]/norm[0, 2]-1)):.12g}\n")
(OUTPUT / "summary.txt").write_text(summary)
print(summary)
