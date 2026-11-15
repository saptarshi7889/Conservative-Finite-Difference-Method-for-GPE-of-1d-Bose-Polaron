"""Plot saved solver output; requires existing NumPy and Matplotlib installations."""
import argparse
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

X_LIMIT = 40.0
SPACE_STRIDE = 100  # Plot spacing 0.1 for simulation dx=0.001.
TIME_STEP = 0.01
SPATIAL_STEP = 0.001
FIGURE_DPI = 160
PROFILE_TIMES = (0, 1, 5, 10, 20, 25)
FIGURE_SIZE = (8, 5)
ANIMATION_SIZE = (12, 4.5)
ANIMATION_FPS = 10
ANIMATION_DPI = 100
AXIS_PADDING = 0.05


def open_text(path):
    return lzma.open(path, "rt") if path.suffix == ".xz" else path.open()


def load_diagnostic(name):
    path = args.results / name
    if not path.exists():
        path = path.with_suffix(path.suffix + ".xz")
    with open_text(path) as stream:
        return np.loadtxt(stream, ndmin=2)


def snapshot_step(path):
    return int(path.name.split("_")[1].split(".")[0])

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results", type=Path)
args = parser.parse_args()
output = Path(__file__).resolve().parent
velocity = load_diagnostic("velocity.dat")
norm = load_diagnostic("norm.dat")
convergence = load_diagnostic("convergence.dat")
with lzma.open(output / "velocity.csv.xz", "wt", preset=9 | lzma.PRESET_EXTREME) as stream:
    np.savetxt(stream, velocity, delimiter=",",
               header="time,real_velocity,imaginary_velocity", comments="")
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
ax.plot(velocity[:, 0], velocity[:, 1])
ax.set(xlabel="Time", ylabel="Impurity velocity", title="Impurity velocity (L=1000)")
fig.tight_layout()
fig.savefig(output / "velocity.png", dpi=FIGURE_DPI)
plt.close(fig)

snapshots = sorted([*args.results.glob("phi_*.dat"), *args.results.glob("phi_*.dat.xz")],
                   key=snapshot_step)
densities, times = [], []
for path in snapshots:
    # Stream only the requested window; avoid loading a million-point snapshot.
    rows = []
    with open_text(path) as stream:
        for index, line in enumerate(stream):
            if index % SPACE_STRIDE:
                continue
            x, re, im = map(float, line.split())
            if -X_LIMIT < x < X_LIMIT:
                rows.append((x, re * re + im * im))
    data = np.asarray(rows)
    positions = data[:, 0]
    densities.append(data[:, 1])
    times.append(snapshot_step(path) * TIME_STEP)
densities, times = np.asarray(densities), np.asarray(times)
np.savez_compressed(output / "density_window.npz", x=positions, time=times,
                    density=densities)
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
mesh = ax.pcolormesh(positions, times, densities, shading="auto", cmap="viridis")
ax.set_ylim(times[0], times[-1])
fig.colorbar(mesh, ax=ax, label=r"Boson density $|\phi|^2$")
ax.set(xlabel="Position", ylabel="Time", title="Boson density evolution")
fig.tight_layout()
fig.savefig(output / "density_evolution.png", dpi=FIGURE_DPI)
plt.close(fig)
fig, ax = plt.subplots(figsize=FIGURE_SIZE)
for time in PROFILE_TIMES:
    index = np.argmin(abs(times - time))
    ax.plot(positions, densities[index], label=f"t={times[index]:g}")
ax.set(xlabel="Position", ylabel=r"Boson density $|\phi|^2$",
       title="Boson density profiles")
ax.legend()
fig.tight_layout()
fig.savefig(output / "density_profiles.png", dpi=FIGURE_DPI)
plt.close(fig)
# Both panels advance at exactly the saved configuration times. The density
# is in the impurity frame; the impurity is stationary at the potential site.
fig, (velocity_ax, density_ax) = plt.subplots(1, 2, figsize=ANIMATION_SIZE)
velocity_ax.plot(velocity[:, 0], velocity[:, 1], color="0.85")
history, = velocity_ax.plot([], [], color="tab:blue")
marker, = velocity_ax.plot([], [], "o", color="tab:red")
velocity_ax.set(xlabel="Time", ylabel="Impurity velocity", xlim=(times[0], times[-1]))
profile, = density_ax.plot(positions, densities[0], color="tab:purple")
density_min, density_max = densities.min(), densities.max()
padding = AXIS_PADDING * (density_max - density_min)
density_ax.set(xlabel="Position relative to impurity", ylabel=r"$|\phi|^2$",
               xlim=(-X_LIMIT, X_LIMIT), ylim=(density_min-padding, density_max+padding))
density_ax.axvline(-SPATIAL_STEP, color="0.5", linestyle=":", label="Impurity (x=-dx)")
density_ax.legend()
clock = fig.suptitle("")
fig.tight_layout()


def update(frame):
    time = times[frame]
    selected = velocity[:, 0] <= time + np.finfo(float).eps * max(1, abs(time))
    history.set_data(velocity[selected, 0], velocity[selected, 1])
    marker.set_data([time], [np.interp(time, velocity[:, 0], velocity[:, 1])])
    profile.set_ydata(densities[frame])
    clock.set_text(f"Bose polaron: t = {time:.2f}")
    return history, marker, profile, clock


animation = FuncAnimation(fig, update, frames=len(times), blit=False)
animation.save(output / "evolution.gif", writer=PillowWriter(fps=ANIMATION_FPS), dpi=ANIMATION_DPI)
plt.close(fig)
print(f"Final time: {velocity[-1, 0]:g}")
print(f"Max Picard iterations: {int(convergence[:, 1].max())}")
print(f"Max accepted iterate difference: {convergence[:, 2].max():.12g}")
print(f"Max sampled relative particle-number drift: {np.max(abs(norm[:, 2] / norm[0, 2] - 1)):.12g}")
print(f"Final velocity: {velocity[-1, 1]:.12g}")
