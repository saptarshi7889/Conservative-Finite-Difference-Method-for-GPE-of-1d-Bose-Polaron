"""Compressed restart equivalence and failure checks; standard library only."""
import gzip
import io
import lzma
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

PROGRAM = str(Path(sys.argv[1]).resolve())
CHECKPOINT_MEMBER = 'phi_10.dat'
VELOCITIES = (2, -2)
FINAL_FIELD_TOLERANCE = 2e-13
GRID_POINTS = 128
TRUNCATION_DIVISOR = 2
PADDING_BYTES = 128 * 1024


def member_bytes(path, name):
    with tarfile.open(path) as archive:
        return archive.extractfile(name).read()


def run(root, velocity, checkpoint, output, member='-', invalid=None, success=True):
    command = [PROGRAM, str(root), str(velocity), str(checkpoint), output, member]
    if invalid:
        command.append(invalid)
    result = subprocess.run(command, capture_output=True, text=True)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    if not success:
        assert not (root / output).exists(), 'Invalid input created output'


with tempfile.TemporaryDirectory() as directory:
    for velocity in VELOCITIES:
        root = Path(directory) / str(velocity)
        root.mkdir()
        subprocess.run([PROGRAM, str(root), str(velocity)], check=True)
        archive = root / 'half/configurations.tar.gz'
        checkpoint = member_bytes(archive, CHECKPOINT_MEMBER)
        (root / 'field.dat.gz').write_bytes(gzip.compress(checkpoint))
        (root / 'field.dat.xz').write_bytes(lzma.compress(checkpoint))
        (root / 'padded.gz').write_bytes(gzip.compress(b' ' * PADDING_BYTES + checkpoint))
        with tarfile.open(root / 'checkpoints.tar.xz', 'w:xz') as output:
            for name, data in [('unselected.txt', b'not a field'), (CHECKPOINT_MEMBER, checkpoint)]:
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                output.addfile(entry, io.BytesIO(data))
        fixtures = [(archive, CHECKPOINT_MEMBER), (root / 'field.dat.gz', '-'),
                    (root / 'field.dat.xz', '-'), (root / 'checkpoints.tar.xz', CHECKPOINT_MEMBER),
                    (root / 'padded.gz', '-')]
        for index, (source, member) in enumerate(fixtures):
            before = source.read_bytes()
            name = f'restart_{index}'
            run(root, velocity, source, name, member)
            assert source.read_bytes() == before
            actual = member_bytes(root / name / 'configurations.tar.gz', CHECKPOINT_MEMBER)
            expected = member_bytes(root / 'full/configurations.tar.gz', 'phi_20.dat')
            actual_rows = [list(map(float, line.split())) for line in actual.splitlines()]
            expected_rows = [list(map(float, line.split())) for line in expected.splitlines()]
            assert len(actual_rows) == len(expected_rows) == GRID_POINTS
            assert max(abs(x-y) for a, b in zip(actual_rows, expected_rows)
                       for x, y in zip(a, b)) < FINAL_FIELD_TOLERANCE
            assert (root / name / 'velocity.dat').read_bytes().splitlines()[1:] == (
                root / 'full/velocity.dat').read_bytes().splitlines()[2:]
            assert not list((root / name).glob('phi_*.dat'))
            metadata = (root / name / 'parameters.txt').read_text()
            assert f'restart {source}\n' in metadata
            if member != '-':
                assert f'restart_member {member}\n' in metadata
        assert not list(root.glob('phi_*.dat')), 'Checkpoint was extracted to disk'

        run(root, velocity, archive, 'missing-member', 'absent.dat', success=False)
        run(root, velocity, archive, 'unspecified-member', success=False)
        run(root, velocity, root / 'field.dat.xz', 'raw-member', CHECKPOINT_MEMBER, success=False)
        run(root, velocity, root / 'absent.gz', 'missing-file', success=False)
        for invalid in ('zero-time', 'missing-velocity', 'nonfinite-velocity', 'invalid-grid', 'initial-field'):
            run(root, velocity, archive, invalid, CHECKPOINT_MEMBER, invalid, success=False)

        rows = checkpoint.splitlines(keepends=True)
        malformed = {
            'short': b''.join(rows[:-1]),
            'extra': checkpoint + rows[-1],
            'grid': b'99 1 0\n' + b''.join(rows[1:]),
            'nonfinite': b'-8 nan 0\n' + b''.join(rows[1:]),
            'trailing': checkpoint + b'invalid',
        }
        for name, data in malformed.items():
            source = root / f'{name}.gz'
            source.write_bytes(gzip.compress(data))
            run(root, velocity, source, f'rejected-{name}', success=False)
        for suffix in ('gz', 'xz'):
            source = root / f'truncated.{suffix}'
            data = (root / f'field.dat.{suffix}').read_bytes()
            source.write_bytes(data[:len(data)//TRUNCATION_DIVISOR])
            run(root, velocity, source, f'rejected-truncated-{suffix}', success=False)
        plain = root / 'plain.txt'
        plain.write_bytes(checkpoint)
        run(root, velocity, plain, 'rejected-plain', success=False)
        with tarfile.open(root / 'duplicate.tar.gz', 'w:gz') as output:
            for _ in range(2):
                entry = tarfile.TarInfo(CHECKPOINT_MEMBER)
                entry.size = len(checkpoint)
                output.addfile(entry, io.BytesIO(checkpoint))
        run(root, velocity, root / 'duplicate.tar.gz', 'rejected-duplicate', CHECKPOINT_MEMBER, success=False)
    print('Compressed continuation matches continuous evolution for both velocity signs; invalid inputs rejected.')
