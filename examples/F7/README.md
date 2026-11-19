# Constant force F=7

`main.cpp` uses F=7, gamma=1/400, dt=0.01, L=1500, dx=1/2000,
t_final=10, G=15, m/M=1/3, V_0=0, a uniform initial field, at most
12 Picard iterations, and absolute Euclidean iterate tolerance 1e-8.
For this requested diagnostic experiment, `accept_unconverged_steps=true`:
at the iteration cap, the last finite iterate is accepted even if its residual
exceeds tolerance. Early exit on convergence remains enabled. Nonfinite fields
and residuals still stop the solver. The default library behavior remains strict.
`convergence.dat` records the actual iteration count (at most 12) and residual
at every time step; `parameters.txt` records the explicit acceptance setting.

Velocity is saved every 10 time steps and full complex configurations every 50.
Initial and final states are also saved on a successful run: 101 velocity records
and 21 configurations. No spatial downsampling is applied.

Build the library using the root README, then from the repository root:

```sh
c++ -O3 -std=c++17 -Iinclude -I/path/to/eigen3 examples/F7/main.cpp \
    -Lbuild -lgpe_solver -lz -larchive -o build/F7
./build/F7 results_F7_t10_picard12
```

An optional argument selects a new output directory. Existing directories are
never overwritten. The library writes `configurations.tar.gz` directly,
without expanded configuration files. `phi_STEP.dat` members contain x,
Re(phi), Im(phi); their time is STEP*dt. Velocity columns are time, real and
imaginary velocity. These formats preserve the inputs needed for a synchronized
velocity/density animation.

For comparison, the following CLI command retains strict convergence checking
and stops if a step does not converge:

```sh
./build/gpe --force 7 --velocity 0 --gamma 0.0025 --dt 0.01 \
    --length 1500 --dx 0.0005 --end 10 --coupling 15 \
    --mass-ratio 0.3333333333333333 --iterations 12 --tolerance 1e-8 \
    --velocity-every 10 --save-every 50 --output results_F7_cli
```

The CLI additionally writes plain snapshots; use the library example to avoid
that duplication. To compress all outputs after a run:

```sh
tar -czf results_F7_t10_picard12.tar.gz results_F7_t10_picard12
```

The original 10-iteration run had a first-step residual of 0.35059867998455091.
The repeated run with 12 iterations had a residual of 0.14907334393382896
at t=0.01, still above the requested 1e-8 tolerance.
The earlier strict run failed at its first time step because Picard iteration did
not meet tolerance within 12 iterations. Its archive therefore contains only
the initial configuration, not a trajectory to t=30. The supplied legacy
program silently accepts the last iterate at its iteration cap; the packaged
solver rejects this condition. Only the iteration cap was changed from 10 to 12 with user approval; all
other parameters were retained. The failed outputs are compressed locally in
`results_F7_picard12_failed.tar.gz`, with a separate SHA-256 checksum.

The original force-driven source is preserved at the repository root under
`finite_difference_modified_impurity_Thomas_generalized_adiabatic_without_Matrix_hermitian_vel_massdep_Kamenev_new_cluster_upwind_scheme.cpp`.
Two output defects were repaired: references to an undeclared driven-velocity
stream were removed, and the equilibrium configuration stream is now opened.
Its numerical equations and original parameter defaults remain unchanged.

With an existing Python installation containing Matplotlib, run
`python examples/F7/run.py` from the repository root after building `build/F7`.
The runner updates `results_F7_t10_picard12/velocity.png` and `velocity.pdf`
while running, then compresses and verifies all outputs and writes a SHA-256
checksum. It refuses to overwrite an existing run directory.
