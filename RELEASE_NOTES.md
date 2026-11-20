# v1.3.1 — Python support alongside the C++ library

Python users can now import `bose_polaron_gpe` and call the same C++ numerical
implementation with ordinary Python syntax. Both quenched and driven impurity
simulations support fresh runs and compressed-checkpoint continuation.

## Python interfaces

- `SimulationParameters`, `DrivenSimulationParameters` and `SimulationResult`
  are exposed alongside all four solver functions: `GPE_solver`,
  `GPE_solver_contd`, `GPE_solver_driven` and `GPE_solver_driven_contd`.
- Complex fields use NumPy arrays, paths accept strings or `pathlib.Path`, and
  C++ errors become Python exceptions. Strict Picard stopping is preserved.
- The final field is a read-only NumPy view that retains its C++ owner. Initial
  fields and run settings are copied before evolution; the GIL is released
  while the compiled solver runs.
- `environment-python.yml`, installation instructions, a complete Python
  example and installed-package tests are included. Build dependencies and
  runtime native libraries are required for source installation.

## C++ usage and portability

C++ users can include the four supplied public headers unchanged and call the
library from their own programs by linking to `gpe::solver`. These are interfaces
to a compiled library, not standalone header-only implementations: the full
source and implementation under `src/` are required. Python bindings remain
optional for C++ users.

Eigen 1×1 complex products now use explicit scalar extraction, fixing the
macOS Clang/libc++ compilation error while retaining transpose semantics.
Release packaging follows the reorganized animation folders.

## Downloads and checks

The release provides complete `.tar.gz` and `.zip` source archives, a smaller
Python source distribution (`bose_polaron_gpe-1.3.1.tar.gz`), and `SHA256SUMS`.
This is a source release: no portable wheels, precompiled executables, or PyPI
upload are included. See README.md and docs/PYTHON.md for installation.

Publication waits for GCC/Clang C++ validation, extracted-archive tests, and
Python 3.10/3.13 package tests on Linux and macOS. Tests cover Python/C++
agreement, NumPy ownership, checkpoint formats, restart consistency and Picard
failure handling. These are implementation checks; scientific conservation and
grid/time convergence still need to be established for each study.

License: MIT.

---

# v0.3.1 — New driven-impurity codes and strict Picard convergence

This release adds codes for a **constant-force (driven) impurity**, alongside
those for a constant-momentum (quenched) impurity.

## New driven-impurity library interfaces

- `<gpe/GPE_solver_driven.hpp>` exposes `gpe::GPE_solver_driven` for a fresh run
  or a field supplied in memory.
- `<gpe/GPE_solver_driven_contd.hpp>` exposes `gpe::GPE_solver_driven_contd` for
  continuation directly from a compressed checkpoint.
- `gpe::DrivenSimulationParameters` uses the dimensionless force `F` in place
  of the injection velocity `V_0`. The initial drive velocity is zero and the
  drive contribution is `(m/M)*F*t`, evaluated at absolute simulation time.
- Both interfaces link to the same compiled `gpe::solver` library. A complete
  fresh-run and compressed-continuation example is in `examples/driven/`.
- The two original driven research programs now reside in `legacy/`, alongside
  the two original quenched programs.

## Convergence safeguard

All four library interfaces and all four legacy programs stop if the requested
Picard tolerance is not reached at the iteration limit, or if the final
iterate/error is nonfinite. A failed step is never accepted or recorded.
Convergence on the final allowed iteration is accepted. Library calls throw
an exception; the command-line and legacy programs exit with failure.

**API change:** `accept_unconverged_steps` has been removed from both parameter
structs. Callers using that setting must remove it and handle nonconvergence.

## Documentation and recorded examples

The README has separate quenched and driven sections. The README and project
website include driven-impurity animations for G=0.5, F=0.5 and G=15, F=4,
both with gamma=1/400 and M=3m, plus the associated research papers.

## Download and validation

Download `bose-polaron-gpe-0.3.1.tar.gz` or `.zip`; these are source archives,
not precompiled binaries. Build instructions are in README.md. Requirements
remain C++17, Eigen 3.3+, zlib, libarchive with gzip/XZ support, and Make or
CMake 3.17+; Python 3 runs the tests. Check downloads with `SHA256SUMS`.

Tests cover driven acceleration, both force/coupling signs, in-memory and
compressed restarts, strict nonconvergence handling, and convergence on the
last allowed iteration. Small temporary variants exercise all four legacy
programs. Release CI tests GCC and Clang and the extracted source archive.
These checks do not establish full driven research-run equivalence, exact
interacting conservation, or grid/time convergence. See VALIDATION.md.

License: MIT.

---

# v0.2.0 — Public C++ library and legacy source organization

Use `<gpe/solver.hpp>` and link to `gpe::solver` to run your own simulations.
The README now leads with a minimal C++ example and CMake integration.

- The public compiled library provides parameter configuration, uniform or
  supplied initial fields, in-memory continuation, streamed velocities, and
  incrementally compressed configuration archives.
- Both original C++ research programs are moved unchanged into `legacy/`.
  Comparison/continuation scripts and source packaging use the new paths.
- The completed L=500, V0=0.5 comparison through t=25 found exact agreement at
  all 251 saved complex velocities and in the complete final field.
- The optional CLI, Make/CMake builds, tests, and animation remain available.

This release changes source organization and documentation; it does not change
the numerical update. Code that directly uses the historical source paths must
use their new `legacy/` locations. The public header path remains
`gpe/solver.hpp`; its implementation is compiled from `src/`.

## Download and build

Download `bose-polaron-gpe-0.2.0.tar.gz` or `.zip` and extract it. With an
existing C++17 compiler, Eigen 3.3+, zlib, and libarchive development headers/libraries, run:

```sh
make EIGEN_INCLUDE=/path/to/directory/containing/Eigen
make test
./build/gpe --help
./build/gpe --output results
```

On systems with Eigen under `/usr/include/eigen3`, plain `make` is sufficient.
CMake 3.17+ is an alternative; see README.md. Python 3 is needed for tests;
no Python packages are required for the solver or integration tests.
Eigen, zlib, and libarchive are not bundled. This release contains source, not precompiled executables.
`SHA256SUMS` covers the two named release archives:

```sh
sha256sum -c SHA256SUMS
```

## Validation and limitations

GCC/Make integration tests cover fresh evolution, positive/negative-velocity
restarts, independent output intervals, and invalid-input rejection. Exact
particle-number conservation is not established for interacting drift dynamics:
the default short interacting test has relative drift around 5.5e-5.
Spatial/time convergence and general long-time stability are not established.
Read VALIDATION.md before scientific use. The animation illustrates a run;
it is not a convergence benchmark.

License: MIT. Citation information is provided in CITATION.cff.
