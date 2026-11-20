# Constant-force impurity: fresh run and continuation

The reusable C++17 library supports the drive in both original adiabatic
research programs: the imposed momentum gives a velocity contribution
(m/M) F t, plus the coupled boson momentum contribution.
Include `<gpe/GPE_solver_driven.hpp>` for `gpe::GPE_solver_driven`, or
`<gpe/GPE_solver_driven_contd.hpp>` for `gpe::GPE_solver_driven_contd`.
Set `gpe::DrivenSimulationParameters`, with `F` replacing the previous `V_0`
input; the initial drive velocity is fixed at zero. Link to `gpe::solver`.
Positive G is repulsive; negative G is attractive. The mass parameter is m/M.

## Runnable library example

Build the library from current repository source, using existing dependencies
or the project Conda environment described in the root README. Then run:

```sh
cmake --install build --prefix "$PWD/build/install"
cmake -S examples/driven -B build/driven -DCMAKE_PREFIX_PATH="$PWD/build/install"
cmake --build build/driven --parallel 2
cmake -E chdir build/driven ./simulation
```

`main.cpp` sets F=0.5, G=0.5, gamma=1/400, M=3m and starts at rest.
It uses the small library grid L=16, dx=0.125, dt=0.001: a functionality
example, not a converged reproduction of the L=8500 animation.
It evolves to t=0.02, then continues from compressed `phi_20.dat` to t=0.04.
Outputs are `build/driven/driven-results` and
`build/driven/continued-driven-results`; rerunning requires new output paths.

Edit the named constants in `main.cpp` for another coupling or force.
For the strong-coupling animation, the physical parameters are G=15, F=4,
gamma=1/400 and M=3m. A suitable grid, time step and Picard limit must be
validated separately; changing these physical constants alone does not
reproduce that simulation.

The CLI also accepts the drive, for example:

```sh
./build/gpe --velocity 0 --force 0.5 --coupling 0.5 --gamma 0.0025 --mass-ratio 0.3333333333333333 --output driven-check
```

## Restart requirements

Keep the force and physical/grid parameters unchanged, set `t_initial` to the
absolute checkpoint time, and supply its recorded
complex velocity. The example passes the final velocity directly from the
fresh result and selects the corresponding compressed field. Do not reset
the force clock to zero at restart. See the root continuation guide for files
from an earlier process.

## Original standalone research programs

The two long-named `adiabatic` C++ files in `legacy/` are the original
fresh and continuation implementations. They use Eigen and GNU C++ headers
(including `bits/stdc++.h`), with simulation parameters embedded in `main`.
With GCC and an existing Eigen installation, compile from the repository root:

```sh
g++ -O2 -std=c++17 -I/path/to/eigen3 legacy/finite_difference_modified_impurity_Thomas_generalized_adiabatic_without_Matrix_hermitian_vel_massdep_Kamenev_new_cluster_upwind_scheme.cpp -o build/driven-original
g++ -O2 -std=c++17 -I/path/to/eigen3 legacy/finite_difference_modified_impurity_Thomas_generalized_adiabatic_without_Matrix_hermitian_vel_contd_massdep_Kamenev_new_cluster_upwind_scheme.cpp -o build/driven-original-contd
```

Create `build/` first if it does not exist. Before running, edit `m_tilde`,
`G_tilde`, `gamma`, F, grid/time settings, tolerances and time limit in each
program. The shipped settings use millions of grid points and long time limits.
Run in a dedicated output directory: these programs append velocity/norm logs.
The fresh program starts with phi=1 and zero velocity. The continuation program
requires a plain three-column field (x, Re(phi), Im(phi)), matching parameters,
`t_init`, the manually assigned checkpoint velocity V, and an appropriate `rep`.
Its input filename is constructed from those settings. Pass its directory with
a trailing slash, e.g. `./build/driven-original-contd folder=/path/to/checkpoints/`
from the repository root; use an absolute executable path when running elsewhere.
It does not read compressed checkpoints; the packaged library does.

The original continuation source has a repaired stream declaration and a
filename for its driven-velocity output. Its numerical update is unchanged.
Original-program output cadence, convergence handling and filenames differ
from the packaged library. No driven research-run equivalence claim is made.

## Recorded animations

- [Generalized case](../G05_F05_L8500/README.md): G=0.5, F=0.5,
  gamma=1/400, M=3m, L=8500.
- [Very strongly coupled impurity](../F4/README.md): G=15, F=4,
  gamma=1/400, M=3m.

The render scripts retain original local source-data paths. Watching the supplied
MP4/GIF requires no simulation data; rerendering requires those inputs (or an
edited path to your own matching data) and the optional tools described there.
