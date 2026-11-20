"""Read Fig. 13/14 data in place; cache small arrays, never extract large fields."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

import numpy as np

from load_velocities import nested_archive

ROOT = Path('/run/media/sap/RIdata/Data/Repulsive_impurity_simulation_data/Repulsive_Quenched')
OUTPUT = Path(__file__).resolve().parent / 'quantum_flutter_data'
ARCHIVES = ('Different_coupling.tar.gz', 'Different_masses.tar.gz')
COUPLINGS = (3.0, 1.0, 0.85, 0.75, 0.65, 0.5)
DENSITY_COUPLING = 0.65
MASS_PARAMETER = 0.125
INITIAL_VELOCITY = 3.5
FINAL_TIME = 50.0
WINDOW = 40.0
DISPLAY_SPACING = 0.05
CHUNK_BYTES = 4 * 1024 * 1024
TAR_BUFFER_BYTES = 1024 * 1024
POSITION_TOLERANCE = 0.000051
PARAMETER_TOLERANCE = 1e-8
NUMBER = r'([0-9.eE+-]+)'
PARAMETERS = re.compile(r'G' + NUMBER + 'm' + NUMBER + 'V' + NUMBER)
GRID = re.compile(r'dx' + NUMBER + r'dt' + NUMBER + r'L' + NUMBER + '_')
FIELD_GRID = re.compile(r'time' + NUMBER + r'dt' + NUMBER + r'dx' + NUMBER + r'L' + NUMBER + '_')


def sample_field(stream, dx, length):
    """Reconstruct positions from the integer grid, checking printed coordinates."""
    spacing = 1 / round(1 / dx)
    first = round((length / 2 - WINDOW) / spacing) + 1
    stop = round((length / 2 + WINDOW) / spacing)
    stride = max(1, round(DISPLAY_SPACING / spacing))
    positions = np.arange(first, stop, stride) * spacing - length / 2
    row, pending, selected = 0, b'', []
    while row < stop:
        chunk = stream.read(CHUNK_BYTES)
        if not chunk:
            raise ValueError('Density field ends before requested window')
        data = pending + chunk
        end = data.rfind(b'\n') + 1
        complete, pending = data[:end], data[end:]
        count = complete.count(b'\n')
        if row + count > first:
            lines = complete.splitlines()
            start = max(first, row)
            offset = (-(start - first)) % stride
            selected.extend(lines[index - row] for index in range(start + offset, min(stop, row + count), stride))
        row += count
    field = np.loadtxt(selected)
    assert field.shape == (len(positions), 3) and np.isfinite(field).all()
    assert np.max(abs(field[:, 0] - positions)) < POSITION_TOLERANCE
    return positions, field[:, 1] ** 2 + field[:, 2] ** 2


def scan(filename):
    OUTPUT.mkdir(exist_ok=True)
    records = []
    catalog_path = OUTPUT / (Path(filename).stem + '_catalog.jsonl')

    def visit(archive, parent='', catalog=None):
        for member in archive:
            if not member.isfile():
                continue
            source = parent + member.name
            catalog.write(json.dumps({'name': source, 'bytes': member.size}) + '\n')
            catalog.flush()
            if member.name.endswith('.tar.gz'):
                print(f'Nested archive: {source}', flush=True)
                with archive.extractfile(member) as compressed:
                    with nested_archive(compressed) as inner:
                        visit(inner, source + '::', catalog)
                continue
            name = Path(member.name).name
            match = PARAMETERS.search(name)
            if not match or name.startswith('._'):
                continue
            coupling, mass, initial = map(float, match.groups())
            if coupling not in COUPLINGS or abs(mass - MASS_PARAMETER) > PARAMETER_TOLERANCE or initial != INITIAL_VELOCITY:
                continue
            kind = None
            if name.startswith('impurity_velocity_'):
                grid = GRID.search(name)
                if grid:
                    dx, dt, length = map(float, grid.groups())
                    kind = 'velocity'
            elif name.startswith('Phi_solution_') and coupling == DENSITY_COUPLING:
                grid = FIELD_GRID.search(name)
                if grid:
                    time, dt, dx, length = map(float, grid.groups())
                    if time <= FINAL_TIME:
                        kind = 'density'
            if kind is None:
                continue
            identity = hashlib.sha256((filename + '::' + source).encode()).hexdigest()[:16]
            cache = OUTPUT / (identity + '.npz')
            record = dict(archive=filename, member=source, kind=kind, coupling=coupling,
                          mass_parameter=mass, initial_velocity=initial, dx=dx, dt=dt,
                          length=length, cache=cache.name)
            if not cache.exists():
                with archive.extractfile(member) as stream:
                    if kind == 'velocity':
                        values = np.loadtxt(stream, ndmin=2)
                        assert values.shape[1] >= 2 and np.isfinite(values).all()
                        assert np.all(np.diff(values[:, 0]) > 0)
                        values = values[values[:, 0] <= FINAL_TIME + PARAMETER_TOLERANCE, :2]
                        np.savez_compressed(cache, values=values)
                    else:
                        positions, density = sample_field(stream, dx, length)
                        np.savez_compressed(cache, positions=positions, density=density)
            if kind == 'density':
                record['time'] = time
            else:
                values = np.load(cache)['values']
                record['time_range'] = values[[0, -1], 0].tolist() if len(values) else []
            records.append(record)
            (OUTPUT / (filename + '.json')).write_text(json.dumps(records, indent=2) + '\n')
            print(f'Cached {kind}: {source}', flush=True)

    print(f'Reading {filename}', flush=True)
    with catalog_path.open('w') as catalog:
        with subprocess.Popen(['gzip', '-dc', str(ROOT / filename)], stdout=subprocess.PIPE) as decoder:
            try:
                with tarfile.open(fileobj=decoder.stdout, mode='r|', bufsize=TAR_BUFFER_BYTES) as archive:
                    visit(archive, catalog=catalog)
                assert decoder.wait() == 0
            finally:
                decoder.stdout.close()
                if decoder.poll() is None:
                    decoder.terminate()
                    decoder.wait()
    (OUTPUT / (filename + '.complete')).write_text('Complete archive scan\n')
    print(f'Finished {filename}: {len(records)} useful records', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archives', nargs='*', choices=ARCHIVES)
    args = parser.parse_args()
    for filename in args.archives or ARCHIVES:
        scan(filename)
