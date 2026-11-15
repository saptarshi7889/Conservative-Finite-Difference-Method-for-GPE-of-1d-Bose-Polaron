"""Plot t<=15 velocity directly from velocity.dat and the t=15 density."""
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results_L500_v05'
OUTPUT = Path(__file__).resolve().parent
PLOT_TIME = 15.0
PLOT_STEP = 1500
POINTS = 500000
LENGTH = 500.0
SPACING = 0.001
WINDOW = 20.0
GRID_TOLERANCE = 1e-9
SIZE = (8, 4.5)
DPI = 160

path = RESULTS / f'phi_{PLOT_STEP}.dat'
try:
    stream = path.open()
except FileNotFoundError:
    stream = lzma.open(str(path)+'.xz', 'rt')
with stream:
    field = np.loadtxt(stream)
assert field.shape == (POINTS, 3) and np.isfinite(field).all()
assert np.max(abs(field[:, 0]-(np.arange(POINTS)*SPACING-LENGTH/2))) < GRID_TOLERANCE
velocity = np.loadtxt(RESULTS/'velocity.dat', ndmin=2)
velocity = velocity[velocity[:, 0] <= PLOT_TIME + GRID_TOLERANCE]
assert abs(velocity[-1, 0]-PLOT_TIME) < GRID_TOLERANCE
fig, ax = plt.subplots(figsize=SIZE)
ax.plot(velocity[:, 0], velocity[:, 1])
ax.set(xlabel='Time', ylabel='Impurity velocity', title='L=500: velocity from velocity.dat through t=15')
fig.tight_layout()
fig.savefig(OUTPUT/'velocity_t15.png', dpi=DPI)
plt.close(fig)
for suffix, selected, limit in [('t15', abs(field[:, 0])<WINDOW, WINDOW),
                                 ('t15_full_system', np.ones(len(field), dtype=bool), LENGTH/2)]:
    fig, ax = plt.subplots(figsize=SIZE)
    ax.plot(field[selected, 0], field[selected, 1]**2+field[selected, 2]**2)
    ax.set(xlabel='Position in impurity frame', ylabel=r'$|\psi(x)|^2$',
           title='L=500: density at t=15', xlim=(-limit, limit))
    fig.tight_layout()
    fig.savefig(OUTPUT/f'density_{suffix}.png', dpi=DPI)
    plt.close(fig)
print('t=15 plots ready:', OUTPUT, flush=True)
