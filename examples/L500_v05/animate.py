"""Animate logged velocity with a synchronized density inset using every snapshot."""
import json
import lzma
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results_L500_v05'
OUTPUT = Path(__file__).resolve().parent
SPACING = 0.001
TIME_STEP = 0.01
LENGTH = 500.0
POINTS = 500000
WINDOW = 15.0
EXPECTED_FRAMES = 51
EXPECTED_VELOCITY_ROWS = 251
GRID_TOLERANCE = 1e-9
FPS = 5
FRAME_DURATION_MS = 1000 // FPS
FIGURE_SIZE = (12, 7)
DPI = 120
INSET_BOUNDS = (0.47, 0.49, 0.50, 0.44)
DENSITY_PADDING = 0.06
VELOCITY_PADDING = 0.015
LINE_WIDTH = 2.5
DENSITY_LINE_WIDTH = 1.6
MARKER_SIZE = 7
GRID_ALPHA = 0.18
PNG_COMPRESSION = 9
BACKGROUND = '#fafbfe'
VELOCITY_COLOR = '#1764a0'
DENSITY_COLOR = '#8a3e92'
CURRENT_COLOR = '#d47721'
CURSOR_ALPHA = 0.4
GUIDE_LINE_WIDTH = 1
REFERENCE_DENSITY = 1.0
MAIN_TITLE_SIZE = 19
CLOCK_SIZE = 16
SUBTITLE_SIZE = 10
MP4_QUANTIZER = 18
SAVE_INTERVAL = 50
FINAL_STEP = 2500
SUBPLOT_MARGINS = (0.10, 0.96, 0.12, 0.87)
TITLE_POSITION = (0.10, 0.98)
CLOCK_POSITION = (0.96, 0.925)
SUBTITLE_POSITION = (0.10, 0.90)


def step(path):
    return int(path.name.split('_')[1].split('.')[0])


def main():
    files = sorted(RESULTS.glob('phi_*.dat.xz'), key=step)
    assert len(files) == EXPECTED_FRAMES
    assert [step(path) for path in files] == list(range(0, FINAL_STEP+1, SAVE_INTERVAL))
    assert not list(RESULTS.glob('phi_*.dat'))
    velocity = np.loadtxt(RESULTS/'velocity.dat', ndmin=2)
    assert len(velocity) == EXPECTED_VELOCITY_ROWS and np.isfinite(velocity).all()
    assert np.all(np.diff(velocity[:, 0]) > 0)
    densities, times = [], []
    expected_grid = np.arange(POINTS)*SPACING - LENGTH/2
    selected = abs(expected_grid) < WINDOW
    positions = expected_grid[selected]
    for path in files:
        with lzma.open(path, 'rt') as stream:
            field = np.loadtxt(stream)
        assert field.shape == (POINTS, 3) and np.isfinite(field).all()
        assert np.max(abs(field[:, 0]-expected_grid)) < GRID_TOLERANCE
        density = field[selected, 1]**2 + field[selected, 2]**2
        densities.append(density)
        time = step(path)*TIME_STEP
        assert np.any(abs(velocity[:, 0]-time) < GRID_TOLERANCE)
        times.append(time)
        print(f'Loaded configuration at t={time:g}', flush=True)
    densities = np.asarray(densities)
    times = np.asarray(times)
    assert times[-1] == velocity[-1, 0]
    density_min, density_max = densities.min(), densities.max()
    density_padding = DENSITY_PADDING*(density_max-density_min)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI, facecolor=BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    fig.subplots_adjust(left=SUBPLOT_MARGINS[0], right=SUBPLOT_MARGINS[1],
                        bottom=SUBPLOT_MARGINS[2], top=SUBPLOT_MARGINS[3])
    ax.set(xlim=(times[0], times[-1]),
           ylim=(velocity[:, 1].min()-VELOCITY_PADDING, velocity[:, 1].max()+VELOCITY_PADDING),
           xlabel='Time', ylabel='Impurity velocity')
    ax.grid(alpha=GRID_ALPHA)
    ax.spines[['top', 'right']].set_visible(False)
    history, = ax.plot([], [], color=VELOCITY_COLOR, linewidth=LINE_WIDTH)
    marker, = ax.plot([], [], 'o', color=CURRENT_COLOR, markersize=MARKER_SIZE)
    cursor = ax.axvline(times[0], color=CURRENT_COLOR, alpha=CURSOR_ALPHA, linewidth=GUIDE_LINE_WIDTH)
    inset = ax.inset_axes(INSET_BOUNDS)
    inset.set_facecolor('white')
    inset.set(xlim=(-WINDOW, WINDOW),
              ylim=(density_min-density_padding, density_max+density_padding),
              xlabel='Position relative to impurity', ylabel=r'$|\psi(x)|^2$')
    inset.axhline(REFERENCE_DENSITY, color='0.75', linestyle=':', linewidth=GUIDE_LINE_WIDTH)
    inset.axvline(-SPACING, color='0.75', linestyle=':', linewidth=GUIDE_LINE_WIDTH)
    density_line, = inset.plot(positions, densities[0], color=DENSITY_COLOR,
                               linewidth=DENSITY_LINE_WIDTH)
    inset_title = inset.set_title('Boson density')
    fig.suptitle('Impurity velocity and boson density', x=TITLE_POSITION[0], y=TITLE_POSITION[1], ha='left', fontsize=MAIN_TITLE_SIZE)
    clock = fig.text(*CLOCK_POSITION, '', ha='right', fontsize=CLOCK_SIZE, color=VELOCITY_COLOR)
    fig.text(*SUBTITLE_POSITION, 'L = 500    |    $V_0 = 0.5$    |    Density shown in the impurity frame',
             fontsize=SUBTITLE_SIZE, color='0.4')
    gif_frames = []
    with tempfile.TemporaryDirectory(prefix='gpe-animation-') as directory:
        frames = Path(directory)
        for index, time in enumerate(times):
            mask = velocity[:, 0] <= time + GRID_TOLERANCE
            history.set_data(velocity[mask, 0], velocity[mask, 1])
            marker.set_data([time], [velocity[mask, 1][-1]])
            cursor.set_xdata([time, time])
            density_line.set_ydata(densities[index])
            clock.set_text(f't = {time:.1f}')
            inset_title.set_text(f'Boson density at t = {time:.1f}')
            path = frames / f'frame_{index:03d}.png'
            fig.savefig(path, dpi=DPI, pil_kwargs={'compress_level': PNG_COMPRESSION})
            with Image.open(path) as image:
                gif_frames.append(image.convert('RGB'))
            if index == len(times)-1:
                fig.savefig(OUTPUT/'velocity_density_animation_preview.png', dpi=DPI)
        plt.close(fig)
        video = OUTPUT/'velocity_density_animation.mp4'
        pipeline = ['gst-launch-1.0', '-e', 'multifilesrc', f'location={frames}/frame_%03d.png',
                    'index=0', f'stop-index={len(times)-1}',
                    f'caps=image/png,framerate=(fraction){FPS}/1', '!', 'pngdec', '!',
                    'videoconvert', '!', 'video/x-raw,format=I420', '!', 'x264enc',
                    'pass=qual', f'quantizer={MP4_QUANTIZER}', 'speed-preset=medium', '!',
                    'video/x-h264,profile=high', '!', 'mp4mux', 'faststart=true', '!',
                    'filesink', f'location={video}']
        result = subprocess.run(pipeline, capture_output=True, text=True,
                                env=dict(os.environ, GST_DEBUG='1'))
        (OUTPUT/'animation_encoding.log').write_text(result.stdout+result.stderr)
        if result.returncode:
            raise RuntimeError('MP4 encoding failed; see animation_encoding.log')
        gif_frames[0].save(OUTPUT/'velocity_density_animation.gif', save_all=True,
                           append_images=gif_frames[1:], duration=FRAME_DURATION_MS,
                           loop=0, optimize=False)
    manifest = {'configuration_count': len(files), 'files': [path.name for path in files],
                'times': times.tolist(), 'fps': FPS, 'duration_seconds': len(times)/FPS,
                'density_window': [-WINDOW, WINDOW], 'velocity_source': 'results_L500_v05/velocity.dat',
                'density_definition': 'Re(psi)^2 + Im(psi)^2',
                'description': 'Main velocity history with density inset at the same physical time; fixed axes; every saved configuration used.'}
    (OUTPUT/'animation_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Created MP4 and GIF: {len(times)} frames, {len(times)/FPS:g} seconds, t=0 to {times[-1]:g}', flush=True)


if __name__ == '__main__':
    main()
