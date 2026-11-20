# Detailed build and API guide

[Quick start](../README.md) · [Numerical method](NUMERICS.md)

Run the commands below from the repository root.

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
[`examples/library`](../examples/library) for a complete example. Eigen, zlib, and libarchive
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


## Constant-force library headers

Include `<gpe/GPE_solver_driven.hpp>` for a fresh run or a field in memory;
include `<gpe/GPE_solver_driven_contd.hpp>` for compressed continuation.
Both link to the installed `gpe::solver` target.

```cpp
#include <gpe/GPE_solver_driven_contd.hpp>

// Inside your simulation program:
gpe::DrivenSimulationParameters parameters;
parameters.F = 0.5; // Dimensionless force, replacing the V_0 input.
parameters.output_directory = "driven-results";
auto result = gpe::GPE_solver_driven(parameters);
```

The remaining settings have the same names, defaults and meaning as in
`SimulationParameters`. There is no V_0 or separate `force` setting in this
driven interface: the initial drive velocity is zero and the drive contribution
is `(m/M)*F*t`. The default F is 0.5. The initial field defaults to uniform.
For in-memory continuation, set `initial_field`, `initial_velocity` and
`t_initial` as before, then call `GPE_solver_driven`.

For file continuation, leave `initial_field` empty, supply the absolute
`t_initial` and recorded `initial_velocity`, and call
`GPE_solver_driven_contd(parameters, checkpoint_file, checkpoint_member)`.
Retain F and the physical/grid settings. Compression formats, archive member
selection, output cadence and convergence checks match the undriven API.
The [complete driven example](../examples/driven/main.cpp) performs both calls.
The original standalone adiabatic programs now reside in `legacy/`.
