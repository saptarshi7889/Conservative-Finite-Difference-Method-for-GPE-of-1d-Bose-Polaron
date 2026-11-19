"""Animate logged velocity with a particle-style density inset and an illustrated introduction."""
import json
from itertools import islice
import re
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.colors import to_rgb

ROOT = Path(__file__).resolve().parents[2]
RESULTS = Path('/run/media/sap/RIdata/Data/Repulsive_impurity_simulation_data/Repulsive_Driven/Kamenev_check/F4/dx2000/dt0.0005')
OUTPUT = Path(__file__).resolve().parent
SPACING = 0.0005
LENGTH = 2000.0
POINTS = 4000000
WINDOW = 40.0
SPATIAL_STRIDE = 100
EXPECTED_FRAMES = 50
EXPECTED_VELOCITY_ROWS = 501
GRID_TOLERANCE = 1e-9
POSITION_TOLERANCE = 0.000051
SNAPSHOT_INTERVAL = 1.0
FINAL_TIME = 50.0
FPS = 5
FRAME_DURATION_MS = 1000 // FPS
FIGURE_SIZE = (12, 7)
DPI = 120
INSET_BOUNDS = (0.47, 0.10, 0.50, 0.30)
INSET_FONT_SIZE = 9
INSET_BORDER_COLOR = 'black'
INSET_BORDER_WIDTH = 0.8
BOSON_DIAMETER = 3.2
IMPURITY_DIAMETER = 9.0
IMPURITY_SIZE_RATIO = IMPURITY_DIAMETER / BOSON_DIAMETER
BOSON_DISPLAY_STRIDE = 10
BOSON_COLOR = '#279aab'
IMPURITY_COLOR = '#e15b64'
IMPURITY_EDGE = '#873044'
IMPURITY_HIGHLIGHT = '#ffe0c9'
SPHERE_PIXELS = 96
SPHERE_LIGHT = (-0.4, 0.5, 0.75)
SPHERE_AMBIENT = 0.25
SPHERE_DIFFUSE = 0.75
SPHERE_SPECULAR = 0.65
SPHERE_SHININESS = 28
TRANSITION_FRAMES = 10
EASING_QUADRATIC = 3
EASING_CUBIC = 2
PARTICLE_ZORDER = 5
IMPURITY_ZORDER = 8
DENSITY_ALPHA = 0.35
IMMERSION_FRAMES = 5
INTRO_FRAMES = IMMERSION_FRAMES
IMMERSION_START_HEIGHT = 1.5
ANNOTATION_POSITION = (0.04, 0.86)
INTRO_LABEL_POSITION = (0.03, 0.07)
IMPURITY_LABEL_OFFSET = (16, 18)
DENSITY_PADDING = 0.06
VELOCITY_PADDING = 0.015
LINE_WIDTH = 2.5
DENSITY_LINE_WIDTH = 1.6
MARKER_SIZE = 7
GRID_ALPHA = 0.18
PNG_COMPRESSION = 9
BACKGROUND = '#fafbfe'
VELOCITY_COLOR = '#7025a8'
DENSITY_COLOR = '#279aab'
CURRENT_COLOR = '#d47721'
CURSOR_ALPHA = 0.4
GUIDE_LINE_WIDTH = 1
REFERENCE_DENSITY = 1.0
MAIN_TITLE_SIZE = 19
CLOCK_SIZE = 16
SUBTITLE_SIZE = 10
MP4_QUANTIZER = 18
SUBPLOT_MARGINS = (0.10, 0.96, 0.12, 0.87)
TITLE_POSITION = (0.10, 0.98)
CLOCK_POSITION = (0.96, 0.925)
SUBTITLE_POSITION = (0.10, 0.90)


def sphere_image():
    """Shade a sphere using surface normals and a specular highlight."""
    grid = np.linspace(-1, 1, SPHERE_PIXELS)
    x, y = np.meshgrid(grid, -grid)
    radius_squared = x*x + y*y
    inside = radius_squared <= 1
    z = np.sqrt(np.maximum(0, 1-radius_squared))
    light = np.asarray(SPHERE_LIGHT)
    light /= np.linalg.norm(light)
    diffuse = np.maximum(0, x*light[0] + y*light[1] + z*light[2])
    halfway = light + np.array([0, 0, 1])
    halfway /= np.linalg.norm(halfway)
    specular = np.maximum(0, x*halfway[0]+y*halfway[1]+z*halfway[2])**SPHERE_SHININESS
    rgba = np.zeros((SPHERE_PIXELS, SPHERE_PIXELS, 4))
    rgba[:, :, :3] = np.clip(
        np.asarray(to_rgb(IMPURITY_COLOR))*(SPHERE_AMBIENT+SPHERE_DIFFUSE*diffuse[:, :, None])
        + SPHERE_SPECULAR*specular[:, :, None], 0, 1)
    rgba[:, :, 3] = inside
    return rgba


def snapshot_time(path):
    return float(re.search(r'time([0-9.]+)dt', path.name).group(1))


def main():
    files = sorted((path for path in RESULTS.glob('Phi_solution_*.dat')
                    if snapshot_time(path) <= FINAL_TIME), key=snapshot_time)
    assert len(files) == EXPECTED_FRAMES
    assert np.allclose([snapshot_time(path) for path in files],
                       np.arange(SNAPSHOT_INTERVAL, FINAL_TIME + SNAPSHOT_INTERVAL,
                                 SNAPSHOT_INTERVAL), rtol=0, atol=GRID_TOLERANCE)
    velocity_path, = RESULTS.glob('impurity_velocity_*_massdep_finite_difference.dat')
    continuation_path, = RESULTS.glob('impurity_velocity_*_finite_difference_contd_rep1.dat')
    initial = np.loadtxt(velocity_path, ndmin=2)
    continuation = np.loadtxt(continuation_path, ndmin=2)
    assert np.allclose(initial[-1], continuation[0], rtol=0, atol=GRID_TOLERANCE)
    velocity = np.concatenate((initial, continuation[1:]))
    velocity = velocity[velocity[:, 0] <= FINAL_TIME]

    assert len(velocity) == EXPECTED_VELOCITY_ROWS and np.isfinite(velocity).all()
    assert np.all(np.diff(velocity[:, 0]) > 0)
    densities, times = [], []
    expected_grid = np.arange(POINTS)*SPACING - LENGTH/2
    selected = (expected_grid >= -WINDOW) & (expected_grid < WINDOW)
    positions = expected_grid[selected][::SPATIAL_STRIDE]
    assert positions[0] == -WINDOW
    selected_indices = np.flatnonzero(selected)
    first, stop = int(selected_indices[0]), int(selected_indices[-1]) + 1
    for path in files:
        # Read only the requested window; no expanded configuration files are written.
        # Reconstruct x from grid indices because legacy output rounds its coordinates.
        with path.open() as stream:
            field = np.loadtxt(islice(stream, first, stop, SPATIAL_STRIDE))
        assert field.shape == (len(positions), 3) and np.isfinite(field).all()
        assert np.max(abs(field[:, 0] - positions)) < POSITION_TOLERANCE
        density = field[:, 1]**2 + field[:, 2]**2
        densities.append(density)
        time = snapshot_time(path)
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
    ax.set(xlim=(velocity[0, 0], times[-1]),
           ylim=(velocity[:, 1].min()-VELOCITY_PADDING, velocity[:, 1].max()+VELOCITY_PADDING),
           xlabel=r'$tgn_0/\hbar$', ylabel='Impurity velocity')
    ax.grid(alpha=GRID_ALPHA)
    ax.spines[['top', 'right']].set_visible(False)
    history, = ax.plot([], [], color=VELOCITY_COLOR, linewidth=LINE_WIDTH)
    marker, = ax.plot([], [], 'D', color=CURRENT_COLOR, markersize=MARKER_SIZE)
    cursor = ax.axvline(times[0], color=CURRENT_COLOR, alpha=CURSOR_ALPHA, linewidth=GUIDE_LINE_WIDTH)
    inset = ax.inset_axes(INSET_BOUNDS)
    inset.set_facecolor('white')
    inset.set(xlim=(-WINDOW, WINDOW),
              ylim=(density_min-density_padding, density_max+density_padding),
              xlabel='Position relative to impurity', ylabel=r'$|\psi(x)|^2$')
    inset.tick_params(labelsize=INSET_FONT_SIZE)
    inset.yaxis.label.set_size(INSET_FONT_SIZE)
    inset.xaxis.label.set_size(INSET_FONT_SIZE)
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color(INSET_BORDER_COLOR)
        spine.set_linewidth(INSET_BORDER_WIDTH)
    inset.axhline(REFERENCE_DENSITY, color='0.75', linestyle=':', linewidth=GUIDE_LINE_WIDTH)
    inset.axvline(0.0, color='0.75', linestyle=':', linewidth=GUIDE_LINE_WIDTH)
    density_line, = inset.plot(positions, densities[0], color=DENSITY_COLOR,
                               linewidth=DENSITY_LINE_WIDTH, alpha=DENSITY_ALPHA)
    bead_positions = positions[::BOSON_DISPLAY_STRIDE]
    bosons, = inset.plot(bead_positions, densities[0][::BOSON_DISPLAY_STRIDE],
                          linestyle='None', marker='o', markersize=BOSON_DIAMETER,
                          markeredgewidth=0, color=BOSON_COLOR, zorder=PARTICLE_ZORDER)
    sphere = AnnotationBbox(OffsetImage(sphere_image(), zoom=IMPURITY_DIAMETER/SPHERE_PIXELS),
                            (0, 0), frameon=False, pad=0, zorder=IMPURITY_ZORDER)
    inset.add_artist(sphere)
    impurity_label = inset.annotate('impurity', (0, 0), xytext=IMPURITY_LABEL_OFFSET,
                                    textcoords='offset points', fontsize=INSET_FONT_SIZE,
                                    color=IMPURITY_EDGE,
                                    arrowprops=dict(arrowstyle='-', color=IMPURITY_EDGE))
    inset.text(*ANNOTATION_POSITION, 'BOSON DENSITY', transform=inset.transAxes,
               color=BOSON_COLOR, fontsize=INSET_FONT_SIZE)
    intro_label = inset.text(*INTRO_LABEL_POSITION, '', transform=inset.transAxes, ha='left',
                             color=IMPURITY_EDGE, fontsize=INSET_FONT_SIZE)
    fig.suptitle('Impurity velocity and boson density', x=TITLE_POSITION[0], y=TITLE_POSITION[1], ha='left', fontsize=MAIN_TITLE_SIZE)
    clock = fig.text(*CLOCK_POSITION, '', ha='right', fontsize=CLOCK_SIZE, color=VELOCITY_COLOR)
    fig.text(*SUBTITLE_POSITION, 'L = 2000    |    F = 4    |    $V_0 = 0$    |    Density in the impurity frame',
             fontsize=SUBTITLE_SIZE, color='0.4')
    transition_fractions = np.linspace(0, 1, TRANSITION_FRAMES + 2)[1:-1]
    transition_weights = (EASING_QUADRATIC*transition_fractions**2
                          - EASING_CUBIC*transition_fractions**3)
    assert np.all(np.diff(transition_weights) > 0)
    frame_times = np.concatenate((np.zeros(INTRO_FRAMES),
                                  transition_fractions*times[0], times))
    gif_frames = []
    with tempfile.TemporaryDirectory(prefix='gpe-animation-') as directory:
        frames = Path(directory)
        for index, time in enumerate(frame_times):
            mask = velocity[:, 0] <= time + GRID_TOLERANCE
            history.set_data(velocity[mask, 0], velocity[mask, 1])
            marker.set_data([time], [np.interp(time, velocity[:, 0], velocity[:, 1])])
            cursor.set_xdata([time, time])
            if index < INTRO_FRAMES:
                density = np.full_like(positions, REFERENCE_DENSITY)
                progress = min(index / (IMMERSION_FRAMES - 1), 1.0)
                impurity_height = IMMERSION_START_HEIGHT + progress*(REFERENCE_DENSITY-IMMERSION_START_HEIGHT)
                intro_label.set_text('Illustrated initial state · uniform bath')
            elif index < INTRO_FRAMES + TRANSITION_FRAMES:
                weight = transition_weights[index - INTRO_FRAMES]
                density = (1-weight)*REFERENCE_DENSITY + weight*densities[0]
                impurity_height = np.interp(0.0, positions, density)
                intro_label.set_text('Interpolated transition · illustrative')
            else:
                density = densities[index - INTRO_FRAMES - TRANSITION_FRAMES]
                impurity_height = np.interp(0.0, positions, density)
                intro_label.set_text('')
            density_line.set_ydata(density)
            bosons.set_ydata(density[::BOSON_DISPLAY_STRIDE])
            sphere.xy = (0, impurity_height)
            sphere.xybox = (0, impurity_height)
            impurity_label.xy = (0, impurity_height)
            clock.set_text(r'$tgn_0/\hbar = ' + f'{time:.1f}' + '$')
            path = frames / f'frame_{index:03d}.png'
            fig.savefig(path, dpi=DPI, pil_kwargs={'compress_level': PNG_COMPRESSION})
            with Image.open(path) as image:
                gif_frames.append(image.convert('RGB'))
            if index == 0:
                fig.savefig(OUTPUT/'velocity_density_animation_intro.png', dpi=DPI)
            if index == len(frame_times)-1:
                fig.savefig(OUTPUT/'velocity_density_animation_preview.png', dpi=DPI)
        plt.close(fig)
        video = OUTPUT/'velocity_density_animation.mp4'
        pipeline = ['gst-launch-1.0', '-e', 'multifilesrc', f'location={frames}/frame_%03d.png',
                    'index=0', f'stop-index={len(frame_times)-1}',
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
    manifest = {'configuration_directory': str(RESULTS), 'configuration_count': len(files), 'files': [path.name for path in files],
                'times': times.tolist(), 'frame_times': frame_times.tolist(),
                'intro_frames': INTRO_FRAMES, 'transition_frames': TRANSITION_FRAMES,
                'transition_description': 'Smoothstep blend from uniform bath to first saved density; illustrative interpolation',
                'boson_diameter_points': BOSON_DIAMETER, 'impurity_diameter_points': IMPURITY_DIAMETER, 'intro_description': 'Illustrated t=0 uniform bath; not a saved configuration',
                'fps': FPS, 'duration_seconds': len(frame_times)/FPS,
                'impurity_to_boson_diameter_ratio': IMPURITY_SIZE_RATIO,
                'boson_display_stride': BOSON_DISPLAY_STRIDE,
                'density_window': [-WINDOW, WINDOW],
                'spatial_stride': SPATIAL_STRIDE, 'density_points_per_frame': len(positions), 'velocity_source': [str(velocity_path), str(continuation_path)],
                'density_definition': 'Re(psi)^2 + Im(psi)^2',
                'description': 'Velocity history with particle-style density inset, shaded impurity at x=0, illustrated immersion opening; every saved configuration through t=50 used.'}
    (OUTPUT/'animation_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Created MP4 and GIF: {len(frame_times)} frames, {len(frame_times)/FPS:g} seconds, t=0 to {times[-1]:g}', flush=True)


if __name__ == '__main__':
    main()
