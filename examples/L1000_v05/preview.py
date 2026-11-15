"""Plot the saved t=4 density and, when available, the solver velocity log."""
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results_L1000_v05"
RECOVERED = ROOT / "results_L1000_v05_recovered"
OUTPUT = Path(__file__).resolve().parent
X_LIMIT = 20
FINAL_STEP = 400
SPACING = 0.001
LENGTH = 1000
POINTS = 1000000
GRID_TOLERANCE = 1e-9
PLOT_STRIDE = 10
FIGURE_SIZE = (8, 4.5)
FIGURE_DPI = 160
FINAL_TIME = 4

path = RESULTS / f"phi_{FINAL_STEP}.dat"
if not path.exists():
    path = path.with_suffix(path.suffix + ".xz")
opener = lzma.open if path.suffix == ".xz" else open
with opener(path, "rt") as stream:
    field = np.loadtxt(stream)
assert field.shape == (POINTS, 3)
assert np.isfinite(field).all()
assert np.max(abs(field[:, 0] - (np.arange(POINTS) * SPACING - LENGTH / 2))) < GRID_TOLERANCE
selected = np.flatnonzero(abs(field[:, 0]) < X_LIMIT)[::PLOT_STRIDE]
x = field[selected, 0]
density = field[selected, 1]**2 + field[selected, 2]**2
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
ax.plot(x, density)
ax.set(xlabel="Position in impurity frame", ylabel=r"$|\psi(x)|^2$",
       title="Boson density at t=4", xlim=(-X_LIMIT, X_LIMIT))
fig.tight_layout()
fig.savefig(OUTPUT / "density_t4.png", dpi=FIGURE_DPI)
plt.close(fig)
print("Validated t=4 field and saved density_t4.png", flush=True)
velocity_path = RECOVERED / "velocity.dat"
if not velocity_path.exists():
    velocity_path = RESULTS / "velocity.dat"
if velocity_path.stat().st_size:
    velocity = np.loadtxt(velocity_path, ndmin=2)
    velocity = velocity[velocity[:, 0] <= FINAL_TIME]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.plot(velocity[:, 0], velocity[:, 1], "o-")
    ax.set(xlabel="Time", ylabel="Impurity velocity", title="Velocity from velocity.dat")
    fig.tight_layout()
    fig.savefig(OUTPUT / "velocity_current.png", dpi=FIGURE_DPI)
    print("Saved velocity_current.png from", velocity_path, flush=True)
else:
    print("velocity.dat is empty; no velocity plot generated.", flush=True)
