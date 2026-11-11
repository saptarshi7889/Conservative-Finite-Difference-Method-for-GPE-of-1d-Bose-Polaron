# Finite difference solver for a 1D Bose polaron

C++ solver for the coupled Gross–Pitaevskii equation and impurity velocity in
Eqs. (5) and (6) of [Majumdar and Petković, arXiv:2507.10402v2](https://arxiv.org/abs/2507.10402v2).
The method uses periodic boundaries, a cyclic Thomas solver, and Picard iteration
for the nonlinear update, with a velocity-dependent one-sided spatial stencil.

The command-line program supports both a uniform initial field, `phi(x)=1`, and
continuation from a saved field. The two supplied research programs remain
unchanged in the repository root. `src/numerics.hpp` preserves their shared
numerical routines (copied from the fresh-run program); `src/main.cpp` retains
the fresh-run update equations with a new input/output driver.

**Numerical status:** the packaged update matches the supplied fresh-run code
exactly in a small regression case. This does not establish exact particle-number
conservation: the short interacting example exhibits relative drift of about
`5.5e-5`. Check conservation and convergence for each scientific run. See
[VALIDATION.md](VALIDATION.md).

## Build

Requirements: a C++17 compiler, Eigen 3.3 or newer, and either Make or CMake 3.16+.
Python 3 is needed only for tests. No Python packages are required. Dependencies
must already be installed; neither local build command downloads anything.

```sh
make EIGEN_INCLUDE=/path/to/directory/containing/Eigen
make test
./build/gpe --help
```

On systems with Eigen under `/usr/include/eigen3`, plain `make` is sufficient.
An Eigen include directory contains `Eigen/Core`, rather than being the `Eigen`
directory itself. Paths containing spaces are supported when quoted.

Alternatively, with an installed Eigen CMake package:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
cmake --install build --prefix "$HOME/.local"
```

For a nonstandard Eigen installation, pass `-DEigen3_DIR=/path/to/share/eigen3/cmake`.
Make also supports `make install PREFIX="$HOME/.local"`.
Do not mix the Make and CMake build directories; remove your generated build
directory or choose another CMake directory when switching build systems.

## Fresh run

```sh
./build/gpe --output results
```

The small default example uses length 16, spacing 0.125, time step 0.001, and
20 time steps. It is a functionality check, not a converged research simulation.
The original defaults use 20 million spatial points and many large complex
vectors, requiring several GB of memory; they are deliberately not the CLI defaults.

| Option | Default | Meaning |
| --- | --- | --- |
| `--length` | 16 | Periodic system length |
| `--dx` | 0.125 | Spatial spacing |
| `--dt` | 0.001 | Time step |
| `--end` | 0.02 | Final time, reached exactly |
| `--mass-ratio` | 0.125 | Boson/impurity mass ratio, m/M |
| `--coupling` | -1.5 | Dimensionless impurity coupling |
| `--gamma` | 0.1 | Positive boson interaction parameter |
| `--velocity` | 2 | Original injection velocity V0; retain this for restarts |
| `--tolerance` | 1e-8 | Absolute Euclidean norm tolerance between Picard iterates |
| `--iterations` | 12 | Maximum Picard iterations per time step |
| `--save-every` | 10 | Steps between diagnostics and field snapshots |
| `--output` | results | New output directory; must not already exist |

Options use `--name value`. Length/dx and (end-start)/dt must be integers within
floating-point tolerance. At least three spatial points are required. The CLI
stops with a nonzero exit code on nonconvergence or invalid input; it does not
silently accept a nonconverged iterate. Partial output remains available after a
failed run. The simulation runs a fixed number of steps; the legacy equilibrium
and cluster CPU-time stop conditions are not part of the CLI.

## Continue a run

For the default example, `results/phi_10.dat` is the field at time 0.01.
Read the corresponding real velocity from `results/velocity.dat`:

```sh
awk '$1 == 0.01 { print $2 }' results/velocity.dat
```

Substitute that value for `V_AT_CHECKPOINT` below:

```sh
./build/gpe --restart results/phi_10.dat --start 0.01 \
  --restart-velocity V_AT_CHECKPOINT --end 0.02 --output continued
```

Retain the original physical parameters, grid, time step, injection velocity,
and iteration settings. `--velocity` is always V0, whereas
`--restart-velocity` is the instantaneous velocity at the saved time. The velocity
update depends on both adjacent fields, so the checkpoint velocity must be
provided separately. The checkpoint grid and row count are checked, but the
three-column field file contains no physical-parameter metadata: matching those
parameters and the checkpoint time is the user's responsibility.

Legacy three-column snapshots are accepted when the grid matches. Legacy files
use lower output precision, so restarting from them can introduce rounding error.
The new output uses round-trip double precision. Snapshot step numbers restart
at zero in each output directory.

## Output

| File | Columns / content |
| --- | --- |
| `velocity.dat` | time, Re(V), Im(V) |
| `norm.dat` | time, Euclidean norm of phi, dx × sum(abs(phi)^2) |
| `convergence.dat` | time, accepted iteration count, norm of last iterate difference |
| `phi_STEP.dat` | x, Re(phi), Im(phi), without a header |
| `parameters.txt` | Physical and numerical parameters and restart provenance |

Diagnostics have comment headers starting with `#`. Initial and final states
are always saved. The original norm output is `sqrt(sum(abs(phi)^2))`; the third
column of the new norm file gives the dimensionless particle-number integral.

The grid is `x_i = i*dx - length/2`. The supplied delta-potential routine places
its nonzero element at index `N/2-1`, which is `x=-dx` for an even grid. This
one-cell offset is retained for regression compatibility.

## Equation and units

The coefficient structure implements the dimensionless equations

$$
i\partial_t\phi = \left[-\frac{1+m/M}{2}\partial_x^2
 + |\phi|^2 + \widetilde G\delta(x) + i V(t)\partial_x\right]\phi,
\qquad
V(t)=V_0+\frac{i(m/M)}{\sqrt\gamma}\int\phi^*\partial_x\phi\,dx.
$$

Here x is measured in healing lengths, time in inverse chemical-potential units,
phi in square roots of the background density, and V in sound-speed units.
These conventions correspond to the coefficients of the supplied code.
The kinetic factor `1+m/M` selects Eq. (5), rather than the alternative Eq. (7)
discussed in the paper. The discrete velocity routine uses its original
symmetrization over both adjacent time slices. The inherited source comments
refer to separate numerical notes (Eqs. 4.6/4.7); those notes were not supplied.

## License and citation

MIT license; see [LICENSE](LICENSE). Eigen is an external dependency and is not
vendored. Citation metadata is in [CITATION.cff](CITATION.cff).

Saptarshi Majumdar and Aleksandra Petković, “Relaxation dynamics of a mobile
impurity injected into a one-dimensional Bose gas,” Physical Review A **112**,
053326 (2025), [doi:10.1103/tgl5-rdn5](https://doi.org/10.1103/tgl5-rdn5).
