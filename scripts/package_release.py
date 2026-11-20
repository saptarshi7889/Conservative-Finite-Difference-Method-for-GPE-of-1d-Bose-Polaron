"""Create source-only release archives and SHA-256 checksums (standard library)."""
import argparse
import gzip
import hashlib
from pathlib import Path
import re
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_NAME = 'bose-polaron-gpe'
ARCHIVE_TIMESTAMP = 0
ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
ZIP_COMPRESSION = 9
FILE_MODE = 0o644
PACKAGE_PATHS = (
    'environment.yml', 'examples/library', 'examples/driven', 'examples/F4', 'examples/G05_F05_L8500',
    'finite_difference_modified_impurity_Thomas_generalized_adiabatic_without_Matrix_hermitian_vel_massdep_Kamenev_new_cluster_upwind_scheme.cpp',
    'finite_difference_modified_impurity_Thomas_generalized_adiabatic_without_Matrix_hermitian_vel_contd_massdep_Kamenev_new_cluster_upwind_scheme.cpp',
    'examples/L500_v05/compare_library.py', 'examples/L500_v05/compare_continuation.py', 'examples/L500_v05/COMPARISON.md',
    'CMakeLists.txt', 'Makefile', 'LICENSE', 'CITATION.cff', 'README.md',
    'VALIDATION.md', 'PUBLISHING.md', 'RELEASE_NOTES.md', 'src', 'tests', 'include', 'cmake',
    'scripts/package_release.py', '.github/workflows', 'docs',
    'legacy',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'dist')
    args = parser.parse_args()
    version = re.search(r'project\(bose_polaron_gpe VERSION ([0-9.]+)',
                        (ROOT/'CMakeLists.txt').read_text()).group(1)
    prefix = f'{PACKAGE_NAME}-{version}'
    files = []
    for name in PACKAGE_PATHS:
        path = ROOT/name
        if not path.exists():
            raise SystemExit(f'Missing release input: {name}')
        candidates = path.rglob('*') if path.is_dir() else [path]
        files.extend(item for item in candidates if item.is_file()
                     and '__pycache__' not in item.parts and item.suffix != '.pyc')
    files = sorted(set(files))
    archives = [args.output/(prefix+'.tar.gz'), args.output/(prefix+'.zip')]
    checksums = args.output/'SHA256SUMS'
    if any(path.exists() for path in [*archives, checksums]):
        raise SystemExit('Release output already exists; choose a new --output directory.')
    args.output.mkdir(parents=True, exist_ok=True)
    with archives[0].open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=ARCHIVE_TIMESTAMP) as compressed:
        with tarfile.open(fileobj=compressed, mode='w') as archive:
            for path in files:
                info = archive.gettarinfo(str(path), arcname=f'{prefix}/{path.relative_to(ROOT)}')
                info.uid = info.gid = ARCHIVE_TIMESTAMP
                info.uname = info.gname = ''
                info.mtime = ARCHIVE_TIMESTAMP
                info.mode = FILE_MODE
                with path.open('rb') as stream:
                    archive.addfile(info, stream)
    with zipfile.ZipFile(archives[1], 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=ZIP_COMPRESSION) as archive:
        for path in files:
            info = zipfile.ZipInfo(f'{prefix}/{path.relative_to(ROOT)}', date_time=ZIP_TIMESTAMP)
            info.external_attr = FILE_MODE << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, path.read_bytes(), compresslevel=ZIP_COMPRESSION)
    checksums.write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n' for path in archives))
    print(f'Packaged {len(files)} files for v{version} in {args.output}')


if __name__ == '__main__':
    main()
