# Python interface to the C++ Bose-polaron solver

`bose_polaron_gpe` exposes the same compiled C++17 numerical implementation to
Python and Jupyter. All evolution, Picard iterations, checkpoint reading and
output remain in C++; there is no subprocess for a solver call.

This branch provides a source-buildable package. It has not been uploaded to
PyPI, and portable precompiled wheels have not yet been released. The Python
package version follows this branch's CMake project version.

## Install from this checkout

The separate environment specification includes Python, NumPy, pybind11,
scikit-build-core and the existing C++ dependencies. Creating it downloads
packages into a project environment:

```sh
conda env create -f environment-python.yml
conda activate bose-polaron-gpe-python
python -m pip install --no-build-isolation --no-deps .
```

Run the commands from the repository root. `--no-build-isolation` reuses the
build dependencies in this environment. `--no-deps` avoids downloading Python
runtime dependencies, which the environment already supplies. Do not run these
commands in a global/base environment. Existing compatible project environments
may be reused. For Jupyter, use a kernel running this environment's Python.

Without Conda, a source build requires a C++17 compiler, Eigen 3.3+, zlib,
libarchive with gzip/XZ support, CMake 3.17+, Python development headers,
pybind11, scikit-build-core, and NumPy. Python 3.10+ is required. Installing
from source requires compatible development libraries; a local wheel does not
automatically bundle their runtime libraries.

## Quenched impurity

```python
from pathlib import Path
import bose_polaron_gpe as gpe

parameters = gpe.SimulationParameters()
parameters.V_0 = 0.5
parameters.output_directory = Path("quenched-results")
result = gpe.GPE_solver(parameters)

phi = result.final_field       # Read-only NumPy complex128 vector.
density = abs(phi) ** 2
print(result.final_time, result.final_velocity)
print(result.velocity_file, result.configurations_archive)
```

The default grid and duration are a small functionality check, not a converged
research simulation. Each call needs a new output directory whose parent
already exists. Paths accept strings or `pathlib.Path`; returned paths are
`pathlib.Path` objects.

## Driven impurity

```python
parameters = gpe.DrivenSimulationParameters()
parameters.F = 0.5
parameters.G = 0.5
parameters.mass_ratio = 1.0 / 3.0
parameters.gamma = 1.0 / 400.0
parameters.output_directory = "driven-results"
result = gpe.GPE_solver_driven(parameters)
```

`F` replaces `V_0` for the driven interface. Its initial drive velocity is zero;
the drive is `(m/M)*F*t` at absolute time. Other settings retain their C++ names
and defaults: `length`, `dx`, `dt`, `t_initial`, `t_final`, `G`, `gamma`,
`mass_ratio`, `picard_tolerance`, `max_picard_iterations`,
`config_save_every`, `velocity_save_every`, `initial_velocity`,
`initial_field`, `output_directory`, and `write_plain_snapshots`.
`SimulationParameters` additionally retains its existing `force` field.
Save intervals count time steps. See [the C++ parameter guide](USAGE.md).

## Compressed continuation

After the small driven example above (t=0.02):

```python
parameters.t_initial = result.final_time
parameters.t_final = 0.04
parameters.initial_velocity = result.final_velocity
parameters.initial_field = None
parameters.output_directory = "continued-driven-results"
continued = gpe.GPE_solver_driven_contd(
    parameters, result.configurations_archive, "phi_20.dat"
)
```

For a quenched run, use `GPE_solver_contd` with `SimulationParameters` instead.
Retain all physical/grid settings and F (or the original V_0), and supply the
checkpoint's absolute time and recorded complex velocity. The field file does
not carry that metadata. For a tar archive, provide the exact member name;
for a single `.dat.gz` or `.dat.xz`, omit `checkpoint_member`.
Both gzip and XZ compressed tar archives are supported.

For continuation in memory, assign `parameters.initial_field = result.final_field`
and call the corresponding fresh-run function with `t_initial` and
`initial_velocity` set as above.

## NumPy arrays, ownership and threads

Assign a NumPy vector to `initial_field`; values are copied into the C++ parameter
object. Assign `None` or an empty complex vector to clear it. Reading
`initial_field` returns a separate array, so editing that array does not edit
the parameter object. A supplied field must match the grid and be finite.
The solver copies its settings before running, so later Python edits do not
change an ongoing simulation.

`result.final_field` is a read-only NumPy view of the C++ final field. The array
keeps the result alive, even if the original Python result variable is deleted.
Use `result.final_field.copy()` for an independently owned, writable array.
No velocity or configuration history is accumulated in Python; those outputs
retain the C++ streaming behavior.

The binding releases Python's GIL during C++ evolution. Other Python threads
can run, but this does not parallelize one simulation. Independent simulations
must use distinct output directories. A call blocks its calling Python thread
until completion or failure; live plotting and progress callbacks are not
provided. Prompt cancellation of a long C++ call with Ctrl+C is not guaranteed.

## Errors

```python
try:
    result = gpe.GPE_solver_driven(parameters)
except RuntimeError as error:
    print(f"Simulation stopped: {error}")
```

C++ validation, checkpoint and convergence exceptions become Python exceptions.
A failed Picard step is not accepted or saved. Convergence on the final permitted
iteration is accepted. Existing completed output remains on disk after failure;
use a new output directory to rerun.

## Build and verify distributions

In the activated environment:

```sh
python -m build --no-isolation --outdir dist/python
python tests/python_api.py
python examples/python/simulation.py
```

The example creates four run directories under `python-example-results/` and
requires that directory not to exist. The tests use temporary directories.
Tests can additionally compare Python outputs with an existing CLI executable:

```sh
python tests/python_api.py --cli /absolute/path/to/build/gpe
```

The build produces a source distribution and a local platform/Python-specific
wheel. Before public wheel distribution, bundle and audit native runtime
libraries using platform tooling and test on clean target machines. Do not
assume a wheel built in one Conda environment is portable to machines without
that environment. The C++ build remains available without Python bindings:
`GPE_BUILD_PYTHON` defaults to `OFF`.
