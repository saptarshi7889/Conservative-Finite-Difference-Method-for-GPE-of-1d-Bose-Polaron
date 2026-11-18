# Original research programs

These two C++ programs are preserved byte-for-byte from their original location
in the repository root. They provide the historical fresh-run and continuation
implementations used for numerical comparisons.

For new simulations, include `<gpe/solver.hpp>`, set
`gpe::SimulationParameters`, call `gpe::GPE_solver`, and link to the CMake target
`gpe::solver`. See the [main README](../README.md#3-run-your-own-c-simulation).
The public header is the interface to a compiled library; the implementation
is built from `src/`.

The historical programs contain their own parameter settings, output handling,
and stopping conditions. They are retained as reference sources and are not
built by the standard Make or CMake targets.
