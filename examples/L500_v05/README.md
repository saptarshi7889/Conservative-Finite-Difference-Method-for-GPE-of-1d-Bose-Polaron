# L=500 research comparisons

For a first simulation, use the [main quick start](../../README.md) and
[complete C++ example](../library/main.cpp). This directory contains advanced
research and regression scripts with experiment-specific settings.

## Reproduce the comparison

The two supported comparison runners use L=500, dx=0.001, dt=0.01, V0=0.5,
m/M=0.333, G=1, gamma=0.1, tolerance=1e-8, and at most 10 Picard iterations.
Velocity is saved every 10 time steps and configurations every 50 steps.

1. `compare_library.py`: original fresh-run update versus `gpe::GPE_solver`,
   from t=0 to t=25 (251 saved velocities, 51 configurations per run).
2. `compare_continuation.py`: original continuation update versus
   `gpe::GPE_solver_contd`, from the first comparison's t=20 checkpoint to
   t=25 (51 saved velocities, 11 configurations per run).

See [COMPARISON.md](COMPARISON.md) for commands, inputs, provenance, and the
measured equivalence results. Existing Eigen, zlib, libarchive, NumPy, and
Matplotlib installations are needed; the scripts install nothing. These are
large research runs, separate from the small quick-start example.

Generated result directories, field archives, plots, videos, logs, and build
outputs are excluded from the source tree. Run the fresh comparison before the
continuation comparison to recreate its input checkpoint. Each runner requires
a new output directory and can regenerate reports with `--compare-only`.

The remaining `run.py`, `continue_run.py`, `prepare_continuation.py`, `plot*.py`,
and `animate.py` files are historical helpers for the original experiment;
retain their experiment-specific settings when using them for reproduction.
A pre-rendered [animation](../../docs/media/velocity_density_animation.gif) and
[video](../../docs/media/velocity_density_animation.mp4) are kept under `docs/media/`.
