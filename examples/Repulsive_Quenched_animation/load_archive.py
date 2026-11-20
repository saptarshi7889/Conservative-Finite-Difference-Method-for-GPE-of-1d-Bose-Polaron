"""Cache only the sampled G=1, M=3m, V0=1.5 plotting data, without extraction."""
import json
from pathlib import Path
import re
import subprocess
import tarfile

import numpy as np

ARCHIVE = Path('/run/media/sap/RIdata/Data/Repulsive_impurity_simulation_data/Repulsive_Quenched/Velocity_initial_vel_final.tar.gz')
OUTPUT = Path(__file__).resolve().parent
PREFIX = 'Velocity_initial_vel_final/Actual/V1.5/dx3000/'
SPACING = 1 / 3000
LENGTH = 1000
WINDOW = 40
STRIDE = 150
FINAL_TIME = 40
SNAPSHOT_INTERVAL = 2.5
EXPECTED_COUNT = 16
CHUNK_BYTES = 4 * 1024 * 1024
TAR_BUFFER_BYTES = 1024 * 1024
POSITION_TOLERANCE = 0.000051
TIME_TOLERANCE = 1e-8


def sampled_rows(stream, first, stop):
    row, pending, selected = 0, b'', []
    while row < stop:
        chunk = stream.read(CHUNK_BYTES)
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
    first = round((LENGTH/2-WINDOW)/SPACING) + 1
    stop = round((LENGTH/2+WINDOW)/SPACING)
    positions = np.arange(first, stop, STRIDE)*SPACING-LENGTH/2
    snapshots, velocity_records = [], []
    with subprocess.Popen(['gzip', '-dc', str(ARCHIVE)], stdout=subprocess.PIPE) as decoder:
        try:
            with tarfile.open(fileobj=decoder.stdout, mode='r|', bufsize=TAR_BUFFER_BYTES) as archive:
                for member in archive:
                    if not member.isfile() or not member.name.startswith(PREFIX):
                        continue
                    name = Path(member.name).name
                    match = re.search(r'Phi_solution_herm_G1m0\.333333V1\.5time([0-9.]+)dt0\.01dx0\.000333333L1000_', name)
                    if match and float(match[1]) <= FINAL_TIME:
                        with archive.extractfile(member) as stream:
                            field = np.loadtxt(sampled_rows(stream, first, stop))
                        assert field.shape == (len(positions), 3) and np.isfinite(field).all()
                        assert np.max(abs(field[:, 0]-positions)) < POSITION_TOLERANCE
                        snapshots.append((float(match[1]), member.name, field[:, 1]**2+field[:, 2]**2))
                        print(f'Loaded t={match[1]} ({len(snapshots)}/{EXPECTED_COUNT})', flush=True)
                    elif name.startswith('impurity_velocity_') and name.endswith('.dat'):
                        with archive.extractfile(member) as stream:
                            velocity_records.append((member.name, np.loadtxt(stream, ndmin=2)))
                        print('Loaded velocity', flush=True)
                    if len(snapshots) == EXPECTED_COUNT and velocity_records:
                        break
        finally:
            decoder.stdout.close()
            decoder.terminate()
            decoder.wait()
    snapshots.sort()
    assert len(snapshots) == EXPECTED_COUNT and len(velocity_records) == 1
    times = np.array([item[0] for item in snapshots])
    assert np.allclose(times, np.arange(SNAPSHOT_INTERVAL, FINAL_TIME+SNAPSHOT_INTERVAL, SNAPSHOT_INTERVAL))
    velocity = velocity_records[0][1]
    velocity = velocity[velocity[:, 0] <= FINAL_TIME+TIME_TOLERANCE]
    assert abs(velocity[-1, 0]-FINAL_TIME) < TIME_TOLERANCE
    assert np.isfinite(velocity).all() and np.all(np.diff(velocity[:, 0]) > 0)
    np.savez_compressed(OUTPUT/'plot_data.npz', positions=positions, times=times,
                        densities=[item[2] for item in snapshots], velocity=velocity)
    metadata = dict(archive=str(ARCHIVE), snapshot_members=[item[1] for item in snapshots],
                    velocity_members=[item[0] for item in velocity_records], dx=SPACING,
                    dt=0.01, length=LENGTH, window=[-WINDOW, WINDOW], stride=STRIDE)
    (OUTPUT/'source_manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Cached {len(times)} snapshots, {len(positions)} points each', flush=True)


if __name__ == '__main__':
    main()
