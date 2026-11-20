"""Cache the eight Fig. 8 velocity histories from the original simulation archive."""
import json
from pathlib import Path
import re
import subprocess
import tarfile
import threading
from contextlib import contextmanager

import numpy as np

from load_archive import ARCHIVE, FINAL_TIME, OUTPUT, TAR_BUFFER_BYTES, TIME_TOLERANCE

INITIAL_VELOCITIES = (0.5, 1.0, 1.5, 2.5, 3.5, 4.0, 6.0, 6.5)
PREFIX = 'Velocity_initial_vel_final/Actual/'
NESTED = {
    'V2.5_dx30000.tar.gz': 2.5,
    'V3.5_dx25000.tar.gz': 3.5,
    'V4_dx25000.tar.gz': 4.0,
    'V6_dx60000.tar.gz': 6.0,
    'V6.5_dx60000.tar.gz': 6.5,
}
DIRECT = {
    'V0.5/dx900/dt0.0025_L1600/': 0.5,
    'V1_dx2500/': 1.0,
}
INITIAL_TOLERANCE = 1e-5


@contextmanager
def nested_archive(compressed):
    """Use the existing gzip executable to decode large nested archives."""
    decoder = subprocess.Popen(['gzip', '-dc'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    def feed():
        try:
            while chunk := compressed.read(TAR_BUFFER_BYTES):
                decoder.stdin.write(chunk)
        except BrokenPipeError:
            pass
        finally:
            try:
                decoder.stdin.close()
            except BrokenPipeError:
                pass
    worker = threading.Thread(target=feed)
    worker.start()
    try:
        with tarfile.open(fileobj=decoder.stdout, mode='r|', bufsize=TAR_BUFFER_BYTES) as inner:
            yield inner
    finally:
        decoder.stdout.close()
        decoder.terminate()
        decoder.wait()
        worker.join()


def main():
    records = {value: [] for value in INITIAL_VELOCITIES}
    records[1.5] = [('plot_data.npz (previously verified source)', np.load(OUTPUT/'plot_data.npz')['velocity'])]

    def record(name, stream, initial):
        values = np.loadtxt(stream, ndmin=2)
        assert values.shape[1] >= 2 and np.isfinite(values).all()
        assert np.all(np.diff(values[:, 0]) > 0)
        records[initial].append((name, values))
        print(f'V0={initial:g}: {name}, t={values[0,0]:g}..{values[-1,0]:g}', flush=True)
        np.save(OUTPUT/f'velocity_{initial:g}_{len(records[initial])}.npy', values)

    with subprocess.Popen(['gzip', '-dc', str(ARCHIVE)], stdout=subprocess.PIPE) as decoder:
        try:
            with tarfile.open(fileobj=decoder.stdout, mode='r|', bufsize=TAR_BUFFER_BYTES) as archive:
                for member in archive:
                    if not member.isfile() or not member.name.startswith(PREFIX):
                        continue
                    relative = member.name[len(PREFIX):]
                    if relative in NESTED:
                        initial = NESTED[relative]
                        print(f'Reading nested archive {relative}', flush=True)
                        with archive.extractfile(member) as compressed:
                            with nested_archive(compressed) as inner:
                                for child in inner:
                                    print(f'Nested member: {child.name} ({child.size} bytes)', flush=True)
                                    name = Path(child.name).name
                                    if child.isfile() and name.startswith('impurity_velocity_') and re.search(rf'G1m0\.333333V{initial:g}dx', name):
                                        with inner.extractfile(child) as stream:
                                            record(member.name+'::'+child.name, stream, initial)
                                        if any(v[0, 0] == 0 and v[-1, 0] >= FINAL_TIME for _, v in records[initial]):
                                            break
                    else:
                        for prefix, initial in DIRECT.items():
                            if relative.startswith(prefix) and Path(relative).name.startswith('impurity_velocity_'):
                                with archive.extractfile(member) as stream:
                                    record(member.name, stream, initial)
                    if all(any(v[0, 0] == 0 and v[-1, 0] >= FINAL_TIME for _, v in group) for group in records.values()):
                        break
        finally:
            decoder.stdout.close()
            decoder.terminate()
            decoder.wait()
    arrays, sources = {}, {}
    for initial in INITIAL_VELOCITIES:
        group = records[initial]
        assert group, f'Missing velocity V0={initial}'
        # Later continuation records replace earlier records from their restart time.
        group.sort(key=lambda item: item[1][0, 0])
        merged = group[0][1]
        for _, values in group[1:]:
            merged = np.concatenate((merged[merged[:, 0] < values[0, 0]-TIME_TOLERANCE], values))
        merged = merged[merged[:, 0] <= FINAL_TIME+TIME_TOLERANCE]
        assert abs(merged[0, 0]) < TIME_TOLERANCE and abs(merged[-1, 0]-FINAL_TIME) < TIME_TOLERANCE, f'Incomplete V0={initial}'
        assert abs(merged[0, 1]-initial) < INITIAL_TOLERANCE
        assert np.all(np.diff(merged[:, 0]) > 0)
        arrays[f'v{initial:g}'] = merged
        sources[str(initial)] = [name for name, _ in group]
    np.savez_compressed(OUTPUT/'all_velocities.npz', **arrays)
    (OUTPUT/'velocity_sources.json').write_text(json.dumps(sources, indent=2)+'\n')
    print('Cached all eight velocity histories through t=40', flush=True)


if __name__ == '__main__':
    main()
