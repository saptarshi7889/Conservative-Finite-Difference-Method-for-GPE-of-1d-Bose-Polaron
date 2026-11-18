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
