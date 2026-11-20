# G=0.5, F=0.5, L=8500 animation

`python examples/G05_F05_L8500/load_archive.py` streams F0.5.tar.gz with the
installed gzip command and Python tarfile; it never extracts full member files.
Only the L=8500 run in dx4000 and its impurity velocity files are selected.
Other runs in the mixed archive are skipped. Small sampled plotting arrays
and the archive member manifest are written locally.

`python examples/G05_F05_L8500/animate.py` reuses these arrays to render every
saved configuration through t=500. All archived velocity
continuations are combined in chronological order. Each restart supersedes
the earlier velocity branch from its initial time onward; overlapping older
records are preserved in the source cache but not mixed into that branch. No paper digitization or force rescaling is used.

Density is Re(phi)^2 + Im(phi)^2 in the impurity frame over -250 <= x < 250.
Every 200th spatial point is retained (spacing 0.05; 10,000 points at dx=0.00025).
The density curve uses every sampled point; decorative beads use every 40th
plotted sample. Violet velocity, amber tracer, teal bosons, shaded coral
impurity, black density-panel borders, and the t g n_0 / hbar clock follow the previous
animation. Density is shown in a wide panel below velocity so the complete multi-cycle
velocity curve remains unobscured and moving density features have more space.
Playback is 10.4 fps (1.3 times the previous 8 fps). The illustrated immersion opening is followed
immediately by a labelled smooth transition into the first saved snapshot.

If velocity and density records have unequal end times, available data through the selected endpoint are
shown. Missing velocity is not extrapolated; the marker disappears and the
missing interval is labelled. A terminal velocity record beyond the final
density uses the last saved density with a label.

Outputs are MP4, GIF, initial/final PNG previews, an encoding log and an
animation manifest. Rerunning replaces these generated outputs. The source
archive is not modified. Existing plotting and video tools are reused.

Current render ends at t=500, with 200 saved configurations and 5001
velocity records. The complete source caches remain available for rerendering.
