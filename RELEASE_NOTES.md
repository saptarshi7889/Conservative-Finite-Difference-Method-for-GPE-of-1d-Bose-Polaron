# v0.1.0 — Bose-polaron GPE solver

Initial source release of the C++17 solver for a one-dimensional Bose polaron.

- Periodic finite differences, cyclic Thomas solves, and Picard iteration.
- Uniform initial fields and restart from saved configurations and velocities.
- Independent configuration and velocity save intervals, flushed diagnostics,
  and complete configuration files published by atomic rename.
- Make/CMake builds, integration tests, and documented numerical equations.
- Both original research programs and the synchronized simulation animation.

## Download and build

Download `bose-polaron-gpe-0.1.0.tar.gz` or `.zip` and extract it. Install a
C++17 compiler and Eigen 3.3+ using your preferred environment, then run:

```sh
make EIGEN_INCLUDE=/path/to/directory/containing/Eigen
make test
./build/gpe --help
./build/gpe --output results
```

On systems with Eigen under `/usr/include/eigen3`, plain `make` is sufficient.
CMake 3.16+ is an alternative; see README.md. Python 3 is needed for tests;
no Python packages are required for the solver or integration tests.
Eigen is not bundled. This release contains source, not precompiled executables.
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
