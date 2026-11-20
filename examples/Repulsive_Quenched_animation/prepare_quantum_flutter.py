"""Select complete, finest-grid Fig. 13/14 series from the archive caches."""
import json
from collections import defaultdict

import numpy as np

from load_quantum_flutter import (ARCHIVES, COUPLINGS, DENSITY_COUPLING, FINAL_TIME,
                                 OUTPUT, PARAMETER_TOLERANCE, ROOT)

SNAPSHOT_INTERVAL = 2.5
INITIAL_TOLERANCE = 1e-5


def merge(records):
    parts = [(record, np.load(OUTPUT / record['cache'])['values']) for record in records]
    parts = [(record, values) for record, values in parts if len(values)]
    # A restart supersedes the preceding file from its starting time onward.
    parts.sort(key=lambda item: (item[1][0, 0], item[1][-1, 0], item[0]['member']))
    if not parts:
        return np.empty((0, 2))
    values = parts[0][1]
    for _, continuation in parts[1:]:
        values = np.concatenate((values[values[:, 0] < continuation[0, 0] - PARAMETER_TOLERANCE], continuation))
    return values


def main():
    records = []
    for archive in ARCHIVES:
        assert (OUTPUT / (archive + '.complete')).exists(), f'Archive scan incomplete: {archive}'
        records.extend(json.loads((OUTPUT / (archive + '.json')).read_text()))
    arrays, sources = {}, {}
    for coupling in COUPLINGS:
        groups = defaultdict(list)
        for record in records:
            if record['kind'] == 'velocity' and record['coupling'] == coupling:
                folder = record['member'].rsplit('/', 1)[0]
                groups[(record['archive'], folder, record['dx'], record['dt'], record['length'])].append(record)
        complete = []
        for key, group in groups.items():
            values = merge(group)
            if len(values) and abs(values[0, 0]) < PARAMETER_TOLERANCE and abs(values[-1, 0] - FINAL_TIME) < PARAMETER_TOLERANCE:
                assert abs(values[0, 1] - 3.5) < INITIAL_TOLERANCE
                assert np.all(np.diff(values[:, 0]) > 0) and np.isfinite(values).all()
                complete.append((key, group, values))
        assert complete, f'No complete velocity series for G={coupling}'
        # Highest dx-folder number means smallest spacing. Prefer smaller dt on ties.
        complete.sort(key=lambda item: (item[0][2], item[0][3], -item[0][4]))
        key, group, values = complete[0]
        arrays[f'g{coupling:g}'] = values
        sources[str(coupling)] = group
        print(f'G={coupling:g}: dx={key[2]:g}, dt={key[3]:g}, t={values[0,0]:g}..{values[-1,0]:g}; {key[1]}', flush=True)
    snapshots = [r for r in records if r['kind'] == 'density' and r['coupling'] == DENSITY_COUPLING]
    groups = defaultdict(list)
    for record in snapshots:
        groups[(record['archive'], record['member'].rsplit('/', 1)[0], record['dx'], record['dt'], record['length'])].append(record)
    expected_times = np.arange(SNAPSHOT_INTERVAL, FINAL_TIME + SNAPSHOT_INTERVAL, SNAPSHOT_INTERVAL)
    complete = []
    for key, group in groups.items():
        by_time = {record['time']: record for record in group}
        if all(time in by_time for time in expected_times):
            complete.append((key, [by_time[time] for time in expected_times]))
    assert complete, 'No single G=0.65 density folder contains every snapshot through t=50'
    complete.sort(key=lambda item: (item[0][2], item[0][3], -item[0][4]))
    key, snapshots = complete[0]
    fields = [np.load(OUTPUT / record['cache']) for record in snapshots]
    positions = fields[0]['positions']
    assert all(np.array_equal(positions, field['positions']) for field in fields)
    arrays.update(positions=positions, times=expected_times, densities=np.array([field['density'] for field in fields]))
    assert np.isfinite(arrays['densities']).all()
    np.savez_compressed(OUTPUT / 'plot_data.npz', **arrays)
    metadata = dict(archive_root=str(ROOT), paper='https://arxiv.org/pdf/2507.10402',
                    mass_ratio=8, initial_velocity=3.5, gamma=0.1, final_time=FINAL_TIME,
                    velocity_sources=sources, density_coupling=DENSITY_COUPLING,
                    density_sources=snapshots, density_definition='Re(psi)^2 + Im(psi)^2',
                    selection='Smallest spacing among folders with complete data through t=50; smallest dt on ties')
    (OUTPUT / 'source_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Density: {len(snapshots)} snapshots, dx={key[2]:g}, t={expected_times[0]:g}..{expected_times[-1]:g}', flush=True)


if __name__ == '__main__':
    main()
