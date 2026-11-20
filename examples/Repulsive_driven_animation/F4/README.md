# F4 animation

Run `python examples/F4/animate.py` from the repository root with the existing
NumPy, Matplotlib, Pillow and GStreamer installation used for the F7 animation.
The source directory is defined by `RESULTS` in the script; source data is read
in place and is never modified.

The original and continuation `impurity_velocity_*.dat` files are joined after
checking their shared t=22 record agrees. Only velocity records through t=50
are used. All 50 configurations at integer times t=1 through 50 are shown,
with the matching real impurity velocity. There is no supplied t=0 snapshot.

Density is Re(phi)^2 + Im(phi)^2 in the impurity frame, over -40 <= x < 40.
Every 100th spatial point is plotted (1600 points; spacing 0.05). Coordinates
are reconstructed from dx=0.0005 and L=2000 and checked against saved values.
The density inset uses a faint curve at every sampled point and small teal beads
at every tenth plotted sample (160 beads). A shaded coral impurity at x=0 is
9-point diameter (boson dots: 3.2 points); the velocity curve is deep violet and its tracer is an amber diamond.
The sphere uses diffuse lighting and a specular highlight.
A one-second illustrated t=0 opening shows immersion into a flat density of 1;
this is not a saved simulation configuration. Immediately on landing, a two-second smoothstep density interpolation connects the uniform bath to
the first saved snapshot; this transition is labelled illustrative. Then all
50 snapshots are shown at 5 fps, giving 13 seconds in total. The phase panel is removed.
The clock and main time axis read t g n_0 / hbar.

Outputs in this directory: `velocity_density_animation.mp4`, `.gif`,
`velocity_density_animation_preview.png`, `velocity_density_animation_intro.png`, `animation_manifest.json`, and
`animation_encoding.log`. Rerunning replaces these generated outputs.
