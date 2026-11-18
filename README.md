# Finite difference solver for a 1D Bose polaron

C++ solver for the coupled Gross–Pitaevskii equation and impurity velocity in
Eqs. (5) and (6) of [Majumdar and Petković, arXiv:2507.10402v2](https://arxiv.org/abs/2507.10402v2).
The method uses periodic boundaries, a cyclic Thomas solver, and Picard iteration
for the nonlinear update, with a velocity-dependent one-sided spatial stencil.

For your own simulations, use the public C++ interface in
[`include/gpe/solver.hpp`](include/gpe/solver.hpp): set
`gpe::SimulationParameters` and call `gpe::GPE_solver`. The library supports
both a uniform initial field, `phi(x)=1`, and continuation from a saved field.
It is a **compiled library**: include the public header and link your program
to `gpe::solver`, which supplies the implementation and its dependencies.

The original research programs are preserved unchanged in [`legacy/`](legacy)
for reference and regression comparisons. The supported interface for new
simulation programs is the public header. `src/numerics.hpp` contains internal
numerical routines, `src/solver.cpp` implements the library, and `src/main.cpp`
provides the optional command-line driver.

**Numerical status:** the packaged update matches the supplied fresh-run code
exactly in a small regression case and at all 251 saved velocity points of the
L=500, V0=0.5 run through t=25. This does not establish exact particle-number
conservation: the short interacting example exhibits relative drift of about
`5.5e-5`. Check conservation and convergence for each scientific run. See
[VALIDATION.md](VALIDATION.md).

## Start your own simulation

Use this pattern in your C++ program:

```cpp
#include <gpe/solver.hpp>

int main() {
    gpe::SimulationParameters parameters;
    // Small demonstration; choose physical and numerical parameters for your study.
    parameters.output_directory = "my-simulation";
    const auto result = gpe::GPE_solver(parameters);
}
```

In a CMake project with this repository checked out under `bose-polaron-gpe/`:

```cmake
cmake_minimum_required(VERSION 3.17)
project(my_simulation LANGUAGES CXX)
add_subdirectory(bose-polaron-gpe)
add_executable(simulation main.cpp)
target_link_libraries(simulation PRIVATE gpe::solver)
```

Build with `cmake -S . -B build -DBUILD_TESTING=OFF` and
`cmake --build build --target simulation`.
The [build requirements](#build) must already be available. See
[C++ library](#c-library) below for all parameters, output files, continuation,
and the installed-package alternative. [`examples/library/`](examples/library)
contains a complete example using an installed package.

## Simulation example

Impurity velocity with a synchronized boson-density inset over $-15<x<15$
in the impurity frame. The animation uses all 51 saved configurations from
$t=0$ to $25$, with periodic system length $L=500$, $V_0=0.5$,
$m/M=0.333$, $`\tilde{G}=1`$, and $\gamma=0.1$.
The simulation uses $\Delta x=0.001$, $\Delta t=0.01$, Picard tolerance
$10^{-8}$, and at most 10 iterations per step. Playback is accelerated.

![Impurity velocity with a synchronized boson-density inset](docs/media/velocity_density_animation.gif)

[Download the MP4 video](docs/media/velocity_density_animation.mp4).

## Releases

Versioned source downloads are published on the [GitHub Releases page](https://github.com/saptarshi7889/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/releases).
See [RELEASE_NOTES.md](RELEASE_NOTES.md) for v0.2.0 requirements and numerical
limitations. Eigen is an external dependency; precompiled binaries are not included.

## Build

Requirements: a C++17 compiler, Eigen 3.3 or newer, zlib and libarchive
development headers/libraries, and either Make or CMake 3.17+.
Libarchive must provide built-in gzip and XZ support for compressed continuation.
Python 3 is needed only for tests. No Python packages are required. Dependencies
must already be installed; neither local build command downloads anything.

### Conda setup (Linux and macOS)

The supplied `environment.yml` declares the compiler, Eigen, zlib, libarchive, CMake, Make, and
Python needed to build and test the solver in a dedicated environment:

```sh
conda env create -f environment.yml
conda activate bose-polaron-gpe
cmake -S . -B build-conda -DCMAKE_BUILD_TYPE=Release
cmake --build build-conda --parallel 2
ctest --test-dir build-conda --output-on-failure
./build-conda/gpe --help
```

Environment creation downloads dependencies. The manifest specifies compatible
version ranges, not an exact lockfile. Eigen, zlib, and libarchive are the solver
dependencies; Python is used for the test driver. Plotting/animation scripts under
`examples/` have separate optional dependencies and are not part of this setup.

### Automatic tests

GitHub Actions runs on every pull request opening, reopening, and update, and
on pushes. It builds and tests with GCC and Clang and separately checks the
Conda environment above. CTest runs direct numerical function checks (the open
and periodic Thomas solvers against an independent dense solve, density, delta
placement, and uniform-field velocity) and the simulation integration tests
(restart consistency, uniform-solution conservation, output cadence, and invalid
inputs). Library tests additionally verify archive contents, in-memory continuation,
CLI equivalence, velocity records visible during evolution, and installation
into a separate consumer project. These tests do not establish conservation for
interacting simulations.

To prevent merging failed tests, configure a GitHub branch rule for `master`
requiring these status checks after the workflow has run successfully.

### Build with existing dependencies

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

## C++ library

Include `<gpe/solver.hpp>` and link to the CMake target `gpe::solver`:

```cpp
#include <gpe/solver.hpp>

gpe::SimulationParameters parameters;
parameters.t_final = final_time;
parameters.length = system_length;
parameters.V_0 = injection_velocity;
parameters.mass_ratio = boson_mass / impurity_mass; // m/M, not M/m
parameters.G = coupling;
parameters.gamma = gamma;
parameters.dx = grid_spacing;
parameters.dt = time_step;
parameters.picard_tolerance = tolerance;
parameters.max_picard_iterations = iteration_limit;
parameters.config_save_every = configuration_interval;
parameters.velocity_save_every = velocity_interval;
parameters.output_directory = "my-simulation";
auto result = gpe::GPE_solver(parameters);
```

The physical/numerical defaults are the small CLI defaults listed below.
Saving intervals count **time steps**, not Picard iterations. Both default to
10 independently. The initial state and final state are always recorded.
`length/dx` and `(t_final-t_initial)/dt` must be integers, with at least three
spatial points. The boundary conditions remain periodic.

Velocity rows are appended to a single open `velocity.dat` stream and flushed
immediately every `velocity_save_every` steps, plus the initial and final times.
The function retains no velocity history in memory. Each run requires a new
output directory; it never appends a new run to an existing trajectory file.

`result.configurations_archive` points to `configurations.tar.gz`, a standard
archive of `phi_STEP.dat` files with columns x, Re(phi), Im(phi). Snapshots are
compressed incrementally; memory holds only one snapshot's formatted text,
not the full trajectory. The ustar archive format limits a single snapshot to
8 GiB minus one byte. Extract with `tar -xzf configurations.tar.gz` before using
a snapshot with the CLI restart option. The CLI also writes plain snapshots.
`norm.dat`, `convergence.dat`, and `parameters.txt` accompany the archive.
`result.velocity_file` is the velocity file path. `result.final_field`,
`result.final_velocity` (complex), and `result.final_time` allow further evolution
without rereading files.

Optional initial-condition settings:

| Setting | Default / meaning |
| --- | --- |
| `initial_field` | Empty Eigen complex vector means uniform phi=1; otherwise one finite value per grid point |
| `t_initial` | Zero; starting time for evolution |
| `initial_velocity` | Unset means V_0; optionally supply a complex checkpoint velocity |
| `write_plain_snapshots` | False; also write plain checkpoint files when true |

The grid is `x_j = j*dx - length/2`. For continuation at positive `t_initial`,
provide both `initial_field` and `initial_velocity`, retain the original `V_0`,
and keep the physical/grid parameters consistent. For example:

```cpp
parameters.initial_field = result.final_field;
parameters.initial_velocity = result.final_velocity;
parameters.t_initial = result.final_time;
parameters.t_final = next_final_time;
parameters.output_directory = "continued-simulation";
auto continued = gpe::GPE_solver(parameters);
```

Invalid input, nonlinear nonconvergence, and I/O failures throw exceptions;
callers should catch `const std::exception&`. Partial files remain on failure.
An interrupted archive may be incomplete and require recovery; successful return
means it was finalized. Existing numerical limitations still apply.

For use inside another CMake project, `add_subdirectory(path/to/this/repository)`
and `target_link_libraries(your_target PRIVATE gpe::solver)` suffice. Alternatively,
install the library and use its exported package:

```sh
cmake -S . -B build-library -DCMAKE_BUILD_TYPE=Release
cmake --build build-library --parallel 2
ctest --test-dir build-library --output-on-failure
cmake --install build-library --prefix "$HOME/.local"
cmake -S examples/library -B build-library/example -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build-library/example --parallel 2
```

The consumer uses `find_package(gpe CONFIG REQUIRED)` followed by
`target_link_libraries(your_target PRIVATE gpe::solver)`. See
[`examples/library`](examples/library) for a complete example. Eigen, zlib, and libarchive
must also be discoverable when configuring the consumer.

### Continue directly from a compressed checkpoint

Include `<gpe/GPE_solver_contd.hpp>` and link to the same compiled target
`gpe::solver`. `GPE_solver_contd` loads the initial field from a compressed file
on the same machine, then uses the existing numerical evolution:

```cpp
#include <gpe/GPE_solver_contd.hpp>

// Retain the physical/grid parameters of the original run in parameters.
parameters.t_initial = checkpoint_time;       // Must be positive.
parameters.t_final = next_final_time;
parameters.initial_velocity = checkpoint_velocity;
parameters.output_directory = "continued-simulation";
// Leave parameters.initial_field empty; the checkpoint supplies it.
const std::filesystem::path checkpoint_file = "/path/to/checkpoints/phi_1000.dat.xz";
auto continued = gpe::GPE_solver_contd(parameters, checkpoint_file);
```

For the library's `configurations.tar.gz`, provide the exact snapshot member:

```cpp
auto continued = gpe::GPE_solver_contd(
    parameters, "/path/to/previous-run/configurations.tar.gz", "phi_1000.dat");
```

Supported inputs are single `.dat.gz` or `.dat.xz` fields and gzip/XZ-compressed
tar archives. The file contents determine the format. For tar archives,
`checkpoint_member` is required; for single compressed fields, omit it.
No files are extracted: libarchive decompresses sequentially through a 64 KiB
text buffer into one Eigen field vector. Memory also includes the compression
decoder's buffers/dictionary and the solver's usual working vectors. The full
expanded text and the archive's other fields are never retained in memory.
Reading a later member requires scanning earlier entries; the reader scans the
archive to its end to reject duplicate selected member names.

Supply the **instantaneous checkpoint velocity**, while retaining the original
injection velocity in `V_0`. The three-column field contains no time or physical
parameter metadata: matching `t_initial`, `V_0`, mass ratio, coupling, gamma,
dx, and dt to the original run remains the caller's responsibility. The reader
checks grid coordinates, row count, finite field values, and reports decompression
errors. Legacy low-precision snapshots retain their original rounding error.

Output uses a new directory, compressed configurations, and streamed velocity
records as for `GPE_solver`. Times start at `t_initial`; snapshot step numbers
start at zero for the continued segment. The example under `examples/library/`
demonstrates a complete fresh run followed by compressed-file continuation.

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
| `--save-every` | 10 | Steps between norm diagnostics and field snapshots |
| `--velocity-every` | save-every | Steps between velocity records (flushed immediately) |
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
| `configurations.tar.gz` | Compressed configuration snapshots |
| `velocity.dat` | time, Re(V), Im(V) |
| `norm.dat` | time, Euclidean norm of phi, dx × sum(abs(phi)^2) |
| `convergence.dat` | time, accepted iteration count, norm of last iterate difference |
| `phi_STEP.dat` | x, Re(phi), Im(phi), without a header |
| `parameters.txt` | Physical and numerical parameters and restart provenance |

Diagnostics have comment headers starting with `#`. Initial and final states
are always saved. Plain CLI configurations are published by atomic rename only after the
complete file is closed; diagnostics are flushed at each configuration save. The original norm output is `sqrt(sum(abs(phi)^2))`; the third
column of the new norm file gives the dimensionless particle-number integral.

The grid is `x_i = i*dx - length/2`. The supplied delta-potential routine places
its nonzero element at index `N/2-1`, which is `x=-dx` for an even grid. This
one-cell offset is retained for regression compatibility.

## Equation and units

The coefficient structure implements the dimensionless equations

```math
i\partial_t\phi = \left[-\frac{1+m/M}{2}\partial_x^2
 + |\phi|^2 + \widetilde G\delta(x) + i V(t)\partial_x\right]\phi,
```

```math
V(t)=V_0+\frac{i(m/M)}{\sqrt\gamma}\int\phi^*\partial_x\phi\,dx.
```

Here $x$ is measured in healing lengths, $t$ in inverse chemical-potential units,
$\phi$ in square roots of the background density, and $V$ in sound-speed units.
These conventions correspond to the coefficients of the supplied code.
The kinetic factor $1+m/M$ selects Eq. (5), rather than the alternative Eq. (7)
discussed in the paper. The discrete velocity routine uses symmetrization over both adjacent time slices,
as detailed below.


## Numerical algorithm

This section summarizes *GPE_simulation.pdf*, Chapter 1, printed pages 2–4,
Eqs. (1.3)–(1.7), before Section 1.1, and compares those notes with
[`src/solver.cpp`](src/solver.cpp) and [`src/numerics.hpp`](src/numerics.hpp).
The PDF is not bundled with this repository. Its Eqs. (1.6)/(1.7) describe the
updates referenced as Eqs. (4.6)/(4.7) in the inherited source comments.
Tildes on dimensionless quantities are omitted here; write $\mu=m/M$,
$h=\Delta x$, $\tau=\Delta t$, and $N=L/h$.

### Spatial and temporal discretization

The field is stored on $N$ periodic points, with indices interpreted modulo $N$.
Define the unscaled second difference and the velocity-selected first difference:

```math
(L_1\Phi)_j=\Phi_{j+1}-2\Phi_j+\Phi_{j-1},\qquad
(L_2(V)\Phi)_j=
\begin{cases}
\Phi_{j+1}-\Phi_j,&\mathrm{Re}(V)\geq0,\\
\Phi_j-\Phi_{j-1},&\mathrm{Re}(V)<0.
\end{cases}
```

For the ordered field vector $\Phi=(\Phi_0,\Phi_1,\ldots,\Phi_{N-1})^{\mathsf T}$,
the corresponding matrices are

```math
L_1=
\begin{pmatrix}
-2 & 1 & 0 & \cdots & 1 \\
1 & -2 & 1 & \ddots & 0 \\
0 & 1 & \ddots & \ddots & \vdots \\
\vdots & \ddots & \ddots & -2 & 1 \\
1 & 0 & \cdots & 1 & -2
\end{pmatrix}_{N\times N}.
```

For $\mathrm{Re}(V)\geq0$, the forward-difference matrix is

```math
L_2^{+}=
\begin{pmatrix}
-1 & 1 & 0 & \cdots & 0 \\
0 & -1 & 1 & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & 0 \\
0 & \cdots & 0 & -1 & 1 \\
1 & 0 & \cdots & 0 & -1
\end{pmatrix}_{N\times N}.
```

For $\mathrm{Re}(V)<0$, the backward-difference matrix is

```math
L_2^{-}=
\begin{pmatrix}
1 & 0 & \cdots & 0 & -1 \\
-1 & 1 & 0 & \ddots & 0 \\
0 & -1 & \ddots & \ddots & \vdots \\
\vdots & \ddots & \ddots & 1 & 0 \\
0 & \cdots & 0 & -1 & 1
\end{pmatrix}_{N\times N},\qquad
L_2(V)=
\begin{cases}
L_2^{+},&\mathrm{Re}(V)\geq0,\\
L_2^{-},&\mathrm{Re}(V)<0.
\end{cases}
```

The corner entries implement periodic boundaries. These are the unscaled
matrices shown in the notes: the derivative operators are $L_1/h^2$ and
$L_2(V)/h$. Their dimension is $N\times N$, where $N=L/h$ is the number of
grid points, rather than the physical system length $L$.

The forward stencil for positive velocity follows the sign of the drift in
$\partial_t\phi=V\partial_x\phi+\cdots$, whose transport speed is $-V$.
The delta potential is represented by $D_j=(\widetilde G/h)\delta_{j,j_*}$,
where $j_*=\lfloor N/2\rfloor-1$ is the inherited implementation location
(see the grid-offset note above).

The kinetic and drift terms are averaged between the old and new time slices.
The nonlinear term uses the average of the two densities multiplied by the
average field. With $W=V^n$ held fixed during the step, the implicit update is

```math
i\frac{\Phi^{n+1}-\Phi^n}{\tau}
=-\frac{1+\mu}{4h^2}L_1(\Phi^{n+1}+\Phi^n)
+\frac{iW}{2h}L_2(W)(\Phi^{n+1}+\Phi^n)
+\frac12\left[D+\frac{|\Phi^{n+1}|^2+|\Phi^n|^2}{2}\right]
\odot(\Phi^{n+1}+\Phi^n).
```

Here $\odot$ and the density operations act componentwise. This is a nonlinear
implicit scheme; the one-sided spatial drift stencil is first order in $h$.

### Picard iteration and cyclic solve

For each time step, set $\Psi_0=\Phi^n$ and form

```math
 a=\frac{i\tau(1+\mu)}{4h^2},\qquad b=\frac{\tau W}{2h},\qquad
 A_1=I-aL_1-bL_2(W),\qquad A_2=I+aL_1+bL_2(W).
```

Each Picard iteration solves the linear system

```math
A_1\Psi_{k+1}=A_2\Phi^n
-\frac{i\tau}{2}\left[D+\frac{|\Psi_k|^2+|\Phi^n|^2}{2}\right]
\odot(\Psi_k+\Phi^n).
```

`Temp` stores $A_2\Phi^n$; `Phi_k`, `Phi_tilde`, and `Phi_hat` store
$\Psi_k$, the right-hand side, and $\Psi_{k+1}$, respectively. Only vectors
and scalar stencil coefficients are stored, rather than dense matrices.
Periodicity makes $A_1$ cyclic tridiagonal. `Thomas_tridiagonal_solver_PBC`
uses two ordinary Thomas solves and a rank-one correction, giving $O(N)$
work and storage per Picard iteration.

Accept $\Phi^{n+1}=\Psi_{k+1}$ when
$\|\Psi_{k+1}-\Psi_k\|_2\leq\varepsilon$. This is an absolute, unweighted
Euclidean norm, not a relative error or an equation-residual test. The notes
specify $\varepsilon=10^{-9}$; the packaged CLI defaults to `--tolerance 1e-8`
and at most `--iterations 12` iterations. It stops with an error if convergence
is not reached within that limit.

### Coupled velocity update

After accepting the field, define $S_j=\Phi_j^n+\Phi_j^{n+1}$. The routine
`Velocity_calculation` evaluates the notes' symmetrized expression as

```math
V^{n+1}=V_0+\frac{i\mu}{8\sqrt\gamma}
\sum_{j=0}^{N-1}\left[
S_j^*(S_{j+1}-S_j)-S_j(S_{j+1}^*-S_j^*)\right].
```

The factors of $h$ from the spatial integral and first difference cancel.
The initial fresh-run field is $\Phi_j^0=1$, with $V^0=V_0$. The updated
velocity selects the stencil and coefficients for the **next** time step;
it is not iterated together with the field inside the current Picard loop.

### Conservation and scope

The notes base the nondrift discretization on the conservative finite-difference
scheme of V. A. Trofimov and N. V. Peskov, “Comparison of finite-difference schemes
for the Gross-Pitaevskii equation,” *Mathematical Modelling and Analysis*
**14**(1), 109–126 (2009), their reference [1]. They describe norm conservation
without drift and negligible norm decay with drift. These statements do not
establish exact conservation for the implemented one-sided drift scheme or a
finite Picard tolerance. Monitor $h\sum_j|\Phi_j|^2$ and check spatial, temporal,
and iteration convergence for the parameters of interest. The measured drift
and limits of the existing checks are documented in [VALIDATION.md](VALIDATION.md).

## License and citation

MIT license; see [LICENSE](LICENSE). Eigen is an external dependency and is not
vendored. Citation metadata is in [CITATION.cff](CITATION.cff).

Saptarshi Majumdar and Aleksandra Petković, “Relaxation dynamics of a mobile
impurity injected into a one-dimensional Bose gas,” Physical Review A **112**,
053326 (2025), [doi:10.1103/tgl5-rdn5](https://doi.org/10.1103/tgl5-rdn5).
