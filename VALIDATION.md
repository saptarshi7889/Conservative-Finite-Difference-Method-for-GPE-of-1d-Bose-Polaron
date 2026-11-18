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

## Compressed-file continuation

On 2026-11-18, the new `GPE_solver_contd` interface passed split-run tests using
single gzip/XZ fields and gzip/XZ tar archives, with no files extracted by the
reader. Tests cover both velocity signs, input spanning multiple read buffers,
grid/row-count validation, missing or duplicate archive members, truncated
streams, nonfinite data, and missing continuation parameters. All four CTest
suites passed; `make test` also passed from an extracted source package. The
installed-library example successfully ran a fresh simulation and continued
directly from its compressed archive.

A separate temporary driver retained the original continuation program's
numerical routines and update block. For L=16, dx=0.125, dt=0.001, m/M=0.125,
G=-1.5, gamma=0.1, tolerance=1e-8 and iteration limit=12, evolution from t=0.01
to t=0.02 agreed exactly in the final field with the new interface, for both
V0=+2 and V0=-2. Both began from the same checkpoint field and velocity. This
original-source check was local and is not part of the automated test suite.
These checks do not establish long-time convergence or conservation.

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

## Checks not run

CMake/CTest passed locally for the compressed-continuation changes described
above. Hosted GitHub Actions have not been run for those changes. The CI
definition covers GCC and Clang on Ubuntu once published.
The publication checkout was inspected and fetched from the configured GitHub
remote before integrating the package. Existing tracked files were retained,
and both supplied research sources are byte-for-byte identical to the prepared
package. Hosted CI results are recorded by GitHub Actions after publication;
they are not a local validation result.
