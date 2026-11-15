# L=1000 Bose-polaron test run

Status: stopped at the user’s request. All generated simulation data and previews were deleted.

The simulation runs locally on the laptop. Data are written under
`results_L1000_v05/` in this local repository, excluded from Git. Nothing is
automatically uploaded to GitHub.

Parameters: L=1000, dx=0.001, dt=0.01, review time=25 (stop here for review before continuing), m/M=0.333,
G=1, gamma=0.1, V0=0.5, absolute Picard tolerance=1e-8,
maximum Picard iterations=10. Save interval: **50 time steps (0.5 time units)**.
The initial field is uniform. The original six-iteration cap failed at the first
step; a separate diagnostic converged in seven iterations there.

The earlier run saving every 10 steps was stopped and preserved under
`results_L1000_v05_previous_s10/`. Its output is also being compressed.
The replacement run stops at t=25 for review. It starts at t=0, avoiding checkpoint ambiguity.
The interrupted save-every-50 startup is preserved in `results_L1000_v05_interrupted_s50/`.
Pilot outputs in /tmp are included in background compression.

Run from the repository root, with a new output directory:

```sh
python3 examples/L1000_v05/run.py
```

`run.py` launches the existing solver and watches Linux file-close events.
Completed simulation text files are queued for **lossless XZ extreme compression**
(`xz -9e --threads=1 --check=sha256`). Separate compression workers handle new and earlier output;
files may wait in the compression queue. The original is removed by XZ only after
successful compression. Files still open by the solver are compressed after
closure, including diagnostics and parameters at solver exit. On solver failure,
partial output is preserved and no completed-run figures are generated.
The existing XZ, NumPy, Matplotlib, and Pillow installations are reused.

At t=25, the solver stops and plots and an animation are generated here. To
regenerate from compressed data:

```sh
MPLCONFIGDIR=/tmp/bose-polaron-mpl python3 examples/L1000_v05/plot.py results_L1000_v05
```

- `velocity.png`: impurity velocity versus time.
- `density_evolution.png`: boson density as a function of position and time.
- `density_profiles.png`: density profiles at selected times.
- `evolution.gif`: synchronized velocity history and density profile, 10 frames/s.
- `velocity.csv.xz`: compressed saved velocity values.
- `density_window.npz`: compressed density samples, positions, and times.

Density is calculated as Re(psi)^2+Im(psi)^2 from configuration files, sampled
every 100 grid points (plot spacing 0.1) in the strict window -40<x<40.
This does not change the simulation grid. Animation frames use saved physical
times; playback is accelerated. Axes stay fixed throughout the animation.
The density is in the impurity frame; the inherited potential is at x=-dx,
as explained in the main README.

PNG/GIF images and the NPZ density archive already use compressed formats.
The completed-run report includes iteration and sampled particle-number drift
diagnostics. This single run does not establish spatial or temporal convergence.
