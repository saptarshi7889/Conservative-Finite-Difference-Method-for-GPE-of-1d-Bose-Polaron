"""Stream L=8500 plotting arrays from a mixed-run archive without extraction."""
import json
from pathlib import Path
import re
import subprocess
import tarfile
import numpy as np

ARCHIVE = Path('/run/media/sap/RIdata/Data/Repulsive_impurity_simulation_data/Repulsive_Driven/Diff_G/G0.5/F0.5.tar.gz')
OUTPUT = Path(__file__).resolve().parent
LENGTH = 8500
WINDOW = 250.0
STRIDE = 200
POSITION_TOLERANCE = 0.000051
RECIPROCAL_GRID_TOLERANCE = 1e-5
READ_CHUNK_BYTES = 4 * 1024 * 1024
TAR_BUFFER_BYTES = 1024 * 1024


def sampled_rows(stream, first, stop):
    """Skip large row ranges by counting newlines in chunks, not Python loops."""
    row, pending, selected = 0, b'', []
    while row < stop:
        chunk = stream.read(READ_CHUNK_BYTES)
        if not chunk:
            raise ValueError('Configuration ends before the requested window')
        data = pending + chunk
        end = data.rfind(b'\n') + 1
        complete, pending = data[:end], data[end:]
        count = complete.count(b'\n')
        if row + count > first:
            lines = complete.splitlines()
            start = max(first, row)
            offset = (-(start-first)) % STRIDE
            selected.extend(lines[index-row] for index in range(start+offset, min(stop, row+count), STRIDE))
        row += count
    return selected


def main():
    snapshots, velocity_records, members, grids = [], [], [], set()
    with subprocess.Popen(['gzip', '-dc', str(ARCHIVE)], stdout=subprocess.PIPE) as decoder:
        with tarfile.open(fileobj=decoder.stdout, mode='r|', bufsize=TAR_BUFFER_BYTES) as archive:
            for index, member in enumerate(archive):
                if not member.isfile():
                    continue
                name = Path(member.name).name
                target = re.search(rf'L{LENGTH}(?:_|\.)', name) is not None
                match = re.search(r'Phi_solution_.*time([0-9.]+)dt([0-9.]+)dx([0-9.]+)L([0-9]+)_', name)
                if target and match:
                    time, dt, dx, length = map(float, match.groups())
                    reciprocal = round(1/dx)
                    if abs(dx*reciprocal-1) < RECIPROCAL_GRID_TOLERANCE:
                        dx = 1/reciprocal
                    grids.add((dt, dx, length))
                    assert len(grids) == 1, 'Multiple L=8500 discretizations found; select a run before combining'
                    first = round((length/2-WINDOW)/dx)
                    stop = round((length/2+WINDOW)/dx)
                    positions = np.arange(first, stop, STRIDE)*dx-length/2
                    with archive.extractfile(member) as stream:
                        field = np.loadtxt(sampled_rows(stream, first, stop))
                    assert field.shape == (len(positions), 3) and np.isfinite(field).all()
                    assert np.max(abs(field[:, 0]-positions)) < POSITION_TOLERANCE
                    snapshots.append((time, member.name, positions, field[:, 1]**2+field[:, 2]**2))
                    print(f'Loaded L={LENGTH}, t={time:g}, dx={dx:g}; {len(snapshots)} snapshots', flush=True)
                elif target and name.startswith('impurity_velocity') and name.endswith('.dat'):
                    with archive.extractfile(member) as stream:
                        values = np.loadtxt(stream, ndmin=2)
                    assert np.isfinite(values).all()
                    velocity_records.append((member.name, values))
                    print(f'Velocity: {member.name}; t={values[0,0]:g}..{values[-1,0]:g}, {len(values)} rows', flush=True)
                members.append({'name': member.name, 'size': member.size, 'selected_length': target})
                print(f'Scanned member {index+1}: {member.name}', flush=True)
        while decoder.stdout.read(READ_CHUNK_BYTES):
            pass
        assert decoder.wait() == 0, 'gzip decompression failed'
    snapshots.sort(key=lambda item: item[0])
    assert snapshots and velocity_records
    positions = snapshots[0][2]
    assert all(np.array_equal(item[2], positions) for item in snapshots)
    assert np.all(np.diff([item[0] for item in snapshots]) > 0)
    np.savez_compressed(OUTPUT/'plot_data_w250_s200.npz', times=[item[0] for item in snapshots],
                        positions=positions, densities=[item[3] for item in snapshots])
    for index, (_, values) in enumerate(velocity_records):
        np.save(OUTPUT/f'velocity_w250_s200_{index}.npy', values)
    metadata = {'archive': str(ARCHIVE), 'snapshot_members': [item[1] for item in snapshots],
                'velocity_members': [item[0] for item in velocity_records], 'members': members,
                'window': [-WINDOW, WINDOW], 'stride': STRIDE,
                'grid_dt_dx_length': list(next(iter(grids)))}
    (OUTPUT/'archive_manifest_w250_s200.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Finished: {len(snapshots)} snapshots, {len(velocity_records)} velocity files', flush=True)


if __name__ == '__main__':
    main()
