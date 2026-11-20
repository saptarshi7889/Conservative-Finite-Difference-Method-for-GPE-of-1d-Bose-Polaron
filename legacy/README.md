# Original research programs

The two undriven C++ programs are preserved byte-for-byte from their original
location in the repository root. They provide the historical fresh-run and continuation
implementations used for numerical comparisons.

For new simulations, include `<gpe/solver.hpp>`, set
`gpe::SimulationParameters`, call `gpe::GPE_solver`, and link to the CMake target
`gpe::solver`. See the [main README](../README.md#3-run-your-own-c-simulation).
The public header is the interface to a compiled library; the implementation
is built from `src/`.

The historical programs contain their own parameter settings, output handling,
and stopping conditions. They are retained as reference sources and are not
built by the standard Make or CMake targets.

The two `adiabatic` programs provide the original constant-force fresh run and
continuation. They were moved from the repository root after the continuation
stream-declaration/filename fixes; this move does not alter their contents.
For driven simulations, use `<gpe/GPE_solver_driven.hpp>` and
`<gpe/GPE_solver_driven_contd.hpp>` with `gpe::DrivenSimulationParameters`:
F replaces V_0, and the imposed drive is (m/M) F t. See the
[driven example and standalone instructions](../examples/driven/README.md).
