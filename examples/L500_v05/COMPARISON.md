Run `compare_library.py --eigen /path/to/include/eigen3` from a Python
installation with existing NumPy and Matplotlib. The script installs nothing.
It refuses to overwrite `results_L500_v05_library_comparison/`.

This repeats the documented L=500, V0=0.5 simulation through t=25:
dx=0.001, dt=0.01, m/M=0.333, G=1, gamma=0.1, uniform initial field,
periodic boundaries, Picard tolerance=1e-8 and iteration limit=10.
Configurations from both runs are saved every 50 steps as before. Velocity is
saved every 10 time steps in both implementations: 251 time points, including
t=0 and t=25. Comparison uses only these saved velocities.

The reference executable is generated from the original fresh-run research
source. Its numerical routines and evolution block are retained; only the
parameters, stopping conditions, and output handling are adapted. Integer
step counting avoids an extra step or accumulated time-label rounding.
Nonconvergence causes an error. The original tracked file is unchanged.
The other executable calls the public `gpe::GPE_solver` function, compiled
from the current library source. Both use the same g++, optimization flags,
and Eigen headers, without fast-math.

Outputs in the comparison directory:

- `manifest.json`: parameters, compiler, Eigen path, source hashes and adaptations.
- `reference.cpp`, `library.cpp`: generated drivers used for the runs.
- `reference/velocity.dat`, `library/velocity.dat`: full precision histories.
- `comparison.csv`: every time, both complex velocities, differences and equality.
- `summary.json`: exact equality counts, maximum differences and final field check.
- `velocity_comparison.png`, `comparison_report.pdf`: overlay, difference plot and report.
- `reference/configurations.tar.gz`, `library/configurations.tar.gz`: 51 field
  snapshots each, with archive integrity verified after completion.

The comparison checks aligned times and finite values before assessing exact
equality of parsed double values. It also checks all 500,000 complex final
field values using binary outputs. A completed comparison reports the measured
outcome, including mismatches; equality is never assumed. This is an equivalence
test for these parameters, not a convergence or conservation assessment.

Regenerate the comparison artifacts without rerunning simulations using
`compare_library.py --compare-only`.

Completed local check, 2026-11-18: all 251 saved complex velocities through
t=25 matched exactly at round-trip double precision. Maximum real and imaginary
velocity differences were both zero; both final velocities were
0.28207112249696997. All 500,000 complex final field values also matched exactly.
Both retained gzip archives passed integrity checks and contained the expected
51 configurations, at steps 0,50,...,2500. The plot, CSV, PDF and JSON summary
are stored alongside the two runs. This result applies to the stated parameters
and compiler recorded in `manifest.json`.

Proposed commit: `Add V0=0.5 original-versus-library simulation comparison`.
Summary: add a reproducible comparison runner, retain both compressed field
trajectories and velocity histories, document the completed equivalence check,
and exclude generated run data from Git. No commit was created.
