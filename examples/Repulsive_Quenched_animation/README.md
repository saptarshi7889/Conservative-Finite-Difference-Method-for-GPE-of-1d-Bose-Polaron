# Repulsive quenched animation

First version of the Fig. 8 (left) velocity animation for the single initial
velocity V0/v=1.5, with G tilde=1, M=3m and gamma=0.1, through t=40.
Paper: https://arxiv.org/abs/2507.10402.

Source: `Velocity_initial_vel_final.tar.gz`, `Actual/V1.5/dx3000`.
The grid spacing is 1/3000, dt=0.01, L=1000. All 16 saved density snapshots
at t=2.5, 5, ..., 40 are retained. Density is Re(psi)^2+Im(psi)^2.
Only -40<x<40 is used, sampled every 150th grid point (spacing 0.05).
Coordinates are reconstructed from grid indices and checked against the
rounded coordinates in the source files. Source data are never changed.

The layout and rendering are adapted from the existing driven F4 animation:
violet velocity history, amber tracer, teal density curve and boson beads,
and a shaded coral impurity at x=0. The illustrated uniform-bath opening
and transition to the first snapshot are labelled illustrative.

With the existing NumPy, Matplotlib, Pillow and GStreamer installation:

```bash
python3 load_archive.py
MPLCONFIGDIR=/tmp/gpe-animation-mpl python3 animate.py
```

Archive loading is needed only once. It caches the small sampled arrays in
`plot_data.npz`; later appearance changes require only `animate.py`.
Outputs are `velocity_density_animation.mp4`, `.gif`, intro/final PNGs,
`animation_manifest.json`, and the encoding log. Rerunning the renderer
replaces these generated outputs.
