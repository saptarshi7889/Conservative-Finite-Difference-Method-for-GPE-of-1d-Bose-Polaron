# L=100 run, review at t=10

The generated L=100 data and plots were deleted at the user’s request.

Run locally with `python3 examples/L100_v05/run.py` from the repository root.
The runner refuses existing output and takes an exclusive lock to prevent
simultaneous launches. Data are saved under `results_L100_v05/` (excluded from Git).

Parameters: L=100, dx=0.001, dt=0.01, V0=0.5, m/M=0.333, G=1,
gamma=0.1, Picard limit=10, absolute tolerance=1e-8. Initial field: psi=1.
Velocity is written directly to `velocity.dat` and flushed every 10 time steps.
Configurations and norm are saved every 50 steps, including initial and final states.
Snapshots are written to `.dat.tmp`, closed, and atomically renamed to `.dat`.
One compression worker compresses completed files with XZ -9e and SHA-256,
tests each archive, then removes its raw source. Compression jobs may queue.
Diagnostics remain readable during the run; compression failures retain raw files.

The solver ends at t=10 and generates the two plots below. No process to continue
to t=25 is launched. The final field and logged velocity are the restart checkpoint.
Existing Eigen, NumPy, Matplotlib, and XZ installations are reused.

![Velocity directly from velocity.dat](velocity_t10.png)

![Density from the t=10 configuration, -20<x<20](density_t10.png)

See `summary.txt` after completion for convergence and norm-drift diagnostics.
This single run does not establish spatial or temporal convergence.
