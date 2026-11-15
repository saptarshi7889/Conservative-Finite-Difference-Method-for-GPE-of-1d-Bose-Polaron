# L=500 run, stop at t=10

Local output: `results_L500_v05/`, excluded from Git.
The simulation runs only to t=10; no continuation is launched.

Parameters: L=500, dx=0.001, dt=0.01, V0=0.5, m/M=0.333,
G=1, gamma=0.1, Picard limit=10, absolute tolerance=1e-8.
Initial field: psi=1. Velocity is saved and flushed every 10 steps;
configurations and norm are saved every 50 steps, including initial/final states.

Completed configurations are published by atomic rename. One worker compresses
with XZ -9e and SHA-256, tests each archive, then removes its raw source.
The runner prevents simultaneous launches and refuses existing output.
Compression failures preserve raw files. Existing dependencies are reused.

This folder was restored after accidental deletion while the simulation was
paused. Data in `results_L500_v05/` were retained; the existing process resumes.
The original solver.log was unlinked by the deletion; velocity, norm, convergence,
and parameter files remain available in the output directory.

At completion, the plotting script reads velocity.dat and the t=10 checkpoint.
See summary.txt for convergence and norm-drift diagnostics. A single run does
not establish elimination of boundary effects or spatial/temporal convergence.

![Velocity](velocity_t10.png)

![Density, -20<x<20](density_t10.png)

![Full-system density](density_t10_full_system.png)
