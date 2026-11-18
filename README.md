# GPE solver for a 1D Bose polaron

C++17 library for the coupled Gross–Pitaevskii equation and impurity velocity,
with periodic boundaries, finite differences, and Picard iteration.
Use it for a fresh simulation or to continue from a compressed checkpoint.

## 1. Download

Download **one complete source archive** from release **v0.2.0**:

- **[Download .tar.gz (Linux/macOS)](https://github.com/saptarshi7889/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/releases/download/v0.2.0/bose-polaron-gpe-0.2.0.tar.gz)**
- **[Download .zip](https://github.com/saptarshi7889/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/releases/download/v0.2.0/bose-polaron-gpe-0.2.0.zip)**
- [Release notes and checksums](https://github.com/saptarshi7889/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/releases/tag/v0.2.0)

Both archives contain the same source code. Download the complete archive:
these headers belong to a **compiled library** and require the implementation
in `src/`. The download contains source, not precompiled executables.

## 2. Build and run

Requirements: a C++17 compiler, Eigen 3.3+, zlib, libarchive with gzip/XZ support,
CMake 3.17+, and Python 3 for the tests. The supplied Conda environment keeps
these dependencies separate from your base environment.

For the `.tar.gz` download:

```sh
tar -xzf bose-polaron-gpe-0.2.0.tar.gz
cd bose-polaron-gpe-0.2.0
conda env create -f environment.yml
conda activate bose-polaron-gpe
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
./build/gpe --output results
```

For `.zip`, extract it and enter `bose-polaron-gpe-0.2.0`, then follow the same
commands starting at `conda env create`. If you already have the dependencies,
start at `cmake -S .`. Environment creation downloads packages; build commands
use existing dependencies.

The last command runs a small demonstration from a uniform field to `t=0.02`.
Its defaults are L=16, dx=0.125, dt=0.001, V0=2, m/M=0.125, G=-1.5,
and gamma=0.1. It is a functionality check, not a converged research simulation.
Each run needs a **new output directory**. List options with `./build/gpe --help`.

## 3. Run your own C++ simulation

Choose the public header for your starting state:

| Task | Header | Function |
| --- | --- | --- |
| Fresh run or supplied field in memory | `gpe/solver.hpp` | `gpe::GPE_solver` |
| Continue from a compressed checkpoint on disk | `gpe/GPE_solver_contd.hpp` | `gpe::GPE_solver_contd` |

The complete program in **[examples/library/main.cpp](examples/library/main.cpp)**
first runs a fresh simulation, then continues it directly from its compressed
archive. Build and run it from the repository root after step 2:

```sh
cmake --install build --prefix "$PWD/build/install"
cmake -S examples/library -B build/example -DCMAKE_PREFIX_PATH="$PWD/build/install"
cmake --build build/example --parallel 2
cmake -E chdir build/example ./simulation
```

Its outputs are `build/example/library-results/` and
`build/example/continued-library-results/`. Edit `gpe::SimulationParameters`
in the example to set your physical parameters, grid, time interval, save
intervals, and output directory. See the [full parameter guide](docs/USAGE.md#c-library).

In your own CMake project, use the installed package:

```cmake
cmake_minimum_required(VERSION 3.17)
project(my_simulation LANGUAGES CXX)
find_package(gpe CONFIG REQUIRED)
add_executable(simulation main.cpp)
target_link_libraries(simulation PRIVATE gpe::solver)
```

Pass the library installation directory through `-DCMAKE_PREFIX_PATH` as above.

## 4. Continue from a compressed checkpoint

Include **`<gpe/GPE_solver_contd.hpp>`** and retain the physical/grid parameters
of the original run. Set the checkpoint time and the velocity recorded at that
time; `V_0` remains the original injection velocity.

Inside your simulation program, with those values supplied:

```cpp
// parameters contains the same physical and numerical settings as the original run.
parameters.t_initial = checkpoint_time;       // Must be greater than zero.
parameters.t_final = next_final_time;
parameters.initial_velocity = checkpoint_velocity;
parameters.initial_field.resize(0);          // Load the field from the file below.
parameters.output_directory = "continued-results";
const std::filesystem::path checkpoint_file = "/path/to/previous-run/configurations.tar.gz";

// Select the snapshot matching checkpoint_time inside the archive.
auto continued = gpe::GPE_solver_contd(parameters, checkpoint_file, "phi_20.dat");
```

In the small example above, `phi_20.dat` is the checkpoint at `t=0.02`.
For a single compressed field, omit the archive member argument:

```cpp
const std::filesystem::path checkpoint_file = "/path/to/phi_20.dat.xz";
auto continued = gpe::GPE_solver_contd(parameters, checkpoint_file);
```

Supported inputs: `.dat.gz`, `.dat.xz`, `.tar.gz`, and `.tar.xz`.
The reader streams decompressed text through a 64 KiB buffer into the field
vector, without creating an expanded file. Memory also includes the decoder
buffers and the solver's working vectors.

Keep `length`, `dx`, `dt`, `mass_ratio`, `G`, `gamma`, and `V_0` consistent with
the checkpoint. The field file has no time/parameter metadata, so the caller
supplies the matching time and checkpoint velocity. Output times start at
`t_initial`; configuration step names start at zero for the new segment.
See the [continuation guide](docs/USAGE.md#continue-directly-from-a-compressed-checkpoint).

## Output files

| File | Contents |
| --- | --- |
| `velocity.dat` | Time, real velocity, imaginary velocity |
| `configurations.tar.gz` | Compressed `phi_STEP.dat` fields: x, Re(phi), Im(phi) |
| `parameters.txt` | Parameters and restart provenance |
| `norm.dat` | Field norm and particle-number diagnostic |
| `convergence.dat` | Picard iterations and last iterate difference |

Velocity and configuration save intervals count **time steps**, not Picard
iterations. Initial and final states are always saved. The command-line driver
also writes plain `phi_STEP.dat` snapshots.

## Repository guide

| Location | Purpose |
| --- | --- |
| `include/gpe/` | Public headers for your programs |
| `src/` | Compiled solver implementation and command-line driver |
| `examples/library/` | Complete fresh-run and compressed-continuation example |
| `examples/L500_v05/` | Advanced comparison and research reproduction scripts |
| `tests/` | Numerical, restart, archive, and interface tests |
| `legacy/` | Original research programs, preserved for reference |
| `docs/` | Detailed usage, equations, numerical method, and demo media |

[Detailed build/API guide](docs/USAGE.md) ·
[Equations and numerical method](docs/NUMERICS.md) ·
[Validation and limitations](VALIDATION.md) ·
[Release notes](RELEASE_NOTES.md) · [Maintainer publishing guide](PUBLISHING.md)

The original and packaged implementations matched all 251 saved velocities
through t=25 in the L=500, V0=0.5 comparison. Both continuation implementations
also matched all 51 saved velocities from t=20 to t=25. These checks establish
implementation agreement for those runs. Conservation and grid/time convergence
must be checked for each scientific study; see [VALIDATION.md](VALIDATION.md).

## Impurity velocity and density evolution

Impurity velocity with a synchronized boson-density inset in the impurity frame.

![Impurity velocity and density evolution](docs/media/velocity_density_animation.gif)

[Download the video](docs/media/velocity_density_animation.mp4)

## License and citation

MIT license; see [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff).

Saptarshi Majumdar and Aleksandra Petković, “Relaxation dynamics of a mobile
impurity injected into a one-dimensional Bose gas,” Physical Review A **112**,
053326 (2025), [doi:10.1103/tgl5-rdn5](https://doi.org/10.1103/tgl5-rdn5).
