# Validation and limits

## Local checks

Built with GCC and an existing Eigen installation using the Makefile, and
rechecked in the publication checkout. No
packages were installed. `make test` passes:

- Fresh uniform initialization and finite interacting evolution.
- Split-run restart agreement with continuous evolution for V0=+2 and V0=-2,
  with maximum snapshot difference below 2e-13.
- Unit modulus and particle-number preservation within 1e-9 for the spatially
  uniform uncoupled case (V0=0, G=0).
- Rejection of existing output directories, invalid grids, fractional time-step
  counts, unknown options, missing/truncated checkpoint data, and nonlinear
  iteration exhaustion.

A separate temporary build of the original fresh solver was run with the same
small grid and duration, the equilibrium stop disabled, increased output
precision, and output every 10 steps. The final field at t=0.02 matched the new
program exactly (maximum absolute difference 0). The original source files in
the workspace were not modified. This comparison was performed locally; it is
not part of the automated integration suite.

## L=500 original-versus-library comparison

On 2026-11-18, both the original fresh-run update and the public library were
run from a uniform field through t=25 with L=500, dx=0.001, dt=0.01, V0=0.5,
m/M=0.333, G=1, gamma=0.1, Picard tolerance=1e-8, and at most 10 iterations.
Velocity was saved every 10 time steps and configurations every 50 steps.
Both were compiled with GCC 15.2.0, `-O3 -DNDEBUG`, and Eigen 3.4.0.

All 251 saved complex velocities matched exactly at round-trip double
precision: maximum real and imaginary differences were both zero. The final
velocity in both runs was 0.28207112249696997, and all 500,000 complex final
field values matched exactly. Both retained compressed archives passed gzip
integrity checks and contained the expected 51 configurations.

The reproducible runner is `examples/L500_v05/compare_library.py`; the original
reference source is now under `legacy/`. Full trajectories and comparison
artifacts are local data excluded from Git and release packages. This checks
implementation equivalence for the stated run, not conservation or convergence.

## Conservation finding

For the default interacting test (L=16, dx=0.125, dt=0.001, T=0.02,
m/M=0.125, G=-1.5, gamma=0.1, V0=2), the maximum sampled relative drift of
`dx * sum(abs(phi)^2)` was approximately 5.4938e-5. This is much larger than
machine precision. The new driver reproduces the original solver's evolution;
it does not fix or explain this drift. The interaction of the one-sided
velocity stencil, nonlinear iteration, and conservation should be investigated
before claiming exact conservation for interacting evolution. Reducing the
Picard tolerance alone is not an established solution.

The tests establish packaging and restart behavior on small examples. They do
not establish spatial/time convergence, stability over long times, conservation
for general parameters, or reproduction of the paper's results. Production
parameter sweeps remain the responsibility of the scientific user.

## v0.2.0 release checks

On 2026-11-18, a fresh GCC 15.2.0 Release build passed all three CTest suites
(numerics, library, integration). The new README C++ and CMake example was
compiled and run as a separate consumer project. Source tar/zip contents and
SHA-256 checksums were verified, and the extracted source archive passed
`make test`. Both programs moved into `legacy/` were checked byte-for-byte
against their prior committed versions. Existing dependencies were reused;
none were installed for these checks.

Hosted GitHub Actions have not been run for v0.2.0 during local preparation.
The release workflow validates GCC and Clang builds before publication;
hosted results must be checked after the tag is pushed.
