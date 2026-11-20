"""Animate logged velocity with a particle-style density inset and an illustrated introduction."""
import json
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

OUTPUT = Path(__file__).resolve().parent
SPACING = 1 / 3000
LENGTH = 1000.0
POINTS = 3000000
WINDOW = 40.0
SPATIAL_STRIDE = 150
EXPECTED_FRAMES = 16
GRID_TOLERANCE = 1e-9
POSITION_TOLERANCE = 0.000051
SNAPSHOT_INTERVAL = 2.5
FINAL_TIME = 40.0
BASE_FPS = 5
PLAYBACK_SPEED = 1.5
FPS = BASE_FPS * PLAYBACK_SPEED
FPS_NUMERATOR = 15
FPS_DENOMINATOR = 2
GIF_TIME_QUANTUM_MS = 10
VELOCITY_INITIALS = (0.5, 1.0, 1.5, 2.5, 3.5, 4.0)
VELOCITY_AXIS_MAX = 4.5
VELOCITY_COLORS = ("#e32929", "#c69200", "#2455cf", "#303030", "#7025a8", "#258130", "#29b9e8", "#ef2190")
LEGEND_COLUMNS = 4
LEGEND_FONT_SIZE = 9
FIGURE_SIZE = (12, 7)
DPI = 120
INSET_BOUNDS = (0.47, 0.45, 0.50, 0.30)
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


def main():
    data = np.load(OUTPUT/'plot_data.npz')
    positions, densities, times = data['positions'], data['densities'], data['times']
    velocity = data['velocity']
    velocity_data = np.load(OUTPUT/'all_velocities.npz')
    velocities = [velocity_data[f'v{initial:g}'] for initial in VELOCITY_INITIALS]
    for initial, values in zip(VELOCITY_INITIALS, velocities):
        assert values[0, 0] == 0 and values[-1, 0] == FINAL_TIME
        assert abs(values[0, 1]-initial) < GRID_TOLERANCE
        assert np.all(np.diff(values[:, 0]) > 0) and np.isfinite(values).all()
    metadata = json.loads((OUTPUT/'source_manifest.json').read_text())
    files = metadata['snapshot_members']
    assert len(times) == EXPECTED_FRAMES
    assert np.allclose(times, np.arange(SNAPSHOT_INTERVAL, FINAL_TIME + SNAPSHOT_INTERVAL, SNAPSHOT_INTERVAL))
    assert times[-1] == velocity[-1, 0] == FINAL_TIME
    assert np.all(np.diff(velocity[:, 0]) > 0)
    assert np.isfinite(densities).all() and np.isfinite(velocity).all()
    density_min, density_max = densities.min(), densities.max()
    density_padding = DENSITY_PADDING*(density_max-density_min)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI, facecolor=BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    fig.subplots_adjust(left=SUBPLOT_MARGINS[0], right=SUBPLOT_MARGINS[1],
                        bottom=SUBPLOT_MARGINS[2], top=SUBPLOT_MARGINS[3])
    ax.set(xlim=(velocity[0, 0], times[-1]),
           ylim=(min(values[:, 1].min() for values in velocities)-VELOCITY_PADDING,
                 VELOCITY_AXIS_MAX),
           xlabel=r'$tgn_0/\hbar$', ylabel=r'Impurity velocity $V/v$')
    ax.grid(alpha=GRID_ALPHA)
    ax.spines[['top', 'right']].set_visible(False)
    histories, markers = [], []
    for initial, color in zip(VELOCITY_INITIALS, VELOCITY_COLORS):
        history, = ax.plot([], [], color=color, linewidth=LINE_WIDTH,
                           label=rf'$V_0/v={initial:g}$')
        marker, = ax.plot([], [], 'D', color=CURRENT_COLOR, markersize=MARKER_SIZE)
        histories.append(history)
        markers.append(marker)
    ax.legend(loc='upper right', ncol=LEGEND_COLUMNS, fontsize=LEGEND_FONT_SIZE,
              frameon=False)
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
    inset.text(*ANNOTATION_POSITION, r'BOSON DENSITY · $V_0/v=1.5$ only', transform=inset.transAxes,
               color=BOSON_COLOR, fontsize=INSET_FONT_SIZE)
    intro_label = inset.text(*INTRO_LABEL_POSITION, '', transform=inset.transAxes, ha='left',
                             color=IMPURITY_EDGE, fontsize=INSET_FONT_SIZE)
    fig.suptitle('Impurity velocity and boson density', x=TITLE_POSITION[0], y=TITLE_POSITION[1], ha='left', fontsize=MAIN_TITLE_SIZE)
    clock = fig.text(*CLOCK_POSITION, '', ha='right', fontsize=CLOCK_SIZE, color=VELOCITY_COLOR)
    fig.text(*SUBTITLE_POSITION, r'$\widetilde G=1$    |    $M=3m$    |    $\gamma=0.1$    |    Impurity frame',
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
            for values, history, marker in zip(velocities, histories, markers):
                mask = values[:, 0] <= time + GRID_TOLERANCE
                history.set_data(values[mask, 0], values[mask, 1])
                marker.set_data([time], [np.interp(time, values[:, 0], values[:, 1])])
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
                    f'caps=image/png,framerate=(fraction){FPS_NUMERATOR}/{FPS_DENOMINATOR}', '!', 'pngdec', '!',
                    'videoconvert', '!', 'video/x-raw,format=I420', '!', 'x264enc',
                    'pass=qual', f'quantizer={MP4_QUANTIZER}', 'speed-preset=medium', '!',
                    'video/x-h264,profile=high', '!', 'mp4mux', 'faststart=true', '!',
                    'filesink', f'location={video}']
        result = subprocess.run(pipeline, capture_output=True, text=True,
                                env=dict(os.environ, GST_DEBUG='1'))
        (OUTPUT/'animation_encoding.log').write_text(result.stdout+result.stderr)
        if result.returncode:
            raise RuntimeError('MP4 encoding failed; see animation_encoding.log')
        gif_boundaries = np.rint(np.arange(len(frame_times)+1)*1000/FPS/GIF_TIME_QUANTUM_MS)
        gif_durations = (np.diff(gif_boundaries)*GIF_TIME_QUANTUM_MS).astype(int).tolist()
        gif_frames[0].save(OUTPUT/'velocity_density_animation.gif', save_all=True,
                           append_images=gif_frames[1:], duration=gif_durations,
                           loop=0, optimize=False)
    manifest = {'configuration_archive': metadata['archive'], 'configuration_count': len(files), 'files': files,
                'times': times.tolist(), 'frame_times': frame_times.tolist(),
                'intro_frames': INTRO_FRAMES, 'transition_frames': TRANSITION_FRAMES,
                'transition_description': 'Smoothstep blend from uniform bath to first saved density; illustrative interpolation',
                'boson_diameter_points': BOSON_DIAMETER, 'impurity_diameter_points': IMPURITY_DIAMETER, 'intro_description': 'Illustrated t=0 uniform bath; not a saved configuration',
                'fps': FPS, 'playback_speed_relative_to_first_version': PLAYBACK_SPEED,
                'velocity_initials': list(VELOCITY_INITIALS), 'density_initial_velocity': 1.5,
                'all_velocity_sources': json.loads((OUTPUT/'velocity_sources.json').read_text()), 'duration_seconds': len(frame_times)/FPS,
                'impurity_to_boson_diameter_ratio': IMPURITY_SIZE_RATIO,
                'boson_display_stride': BOSON_DISPLAY_STRIDE,
                'density_window': [-WINDOW, WINDOW],
                'spatial_stride': SPATIAL_STRIDE, 'density_points_per_frame': len(positions), 'velocity_source': metadata['velocity_members'],
                'density_definition': 'Re(psi)^2 + Im(psi)^2',
                'description': 'Velocity history with particle-style density inset, shaded impurity at x=0, illustrated immersion opening; every saved configuration through t=40 used.'}
    (OUTPUT/'animation_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Created MP4 and GIF: {len(frame_times)} frames, {len(frame_times)/FPS:g} seconds, t=0 to {times[-1]:g}', flush=True)


if __name__ == '__main__':
    main()
