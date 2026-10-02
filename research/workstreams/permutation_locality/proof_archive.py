"""Preserve proof sources and selected witnesses in a verified local snapshot.

This is an integrity archive, not a numerical proof replay. It never replaces
an existing archive, changes a witness, or extracts files over the workspace.
Large experiment caches and binaries are intentionally excluded.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import platform
import subprocess
import sys
from datetime import datetime, timezone
import zipfile

ROOT = Path(__file__).resolve().parents[3]
RECEIPTS = (
    'four-bit-dense-adaptive-half-final.json',
    'two-bit-affine-dense-401-d0925.json',
    'gf16-birth-classes-d09-smallgrid-complete.json',
    'gf16-r2-d09-independent-p384-complete.json',
    'gf16-birth-classes-d09-tilt3over16-q97.json',
    'gf16-r4-d099-assembly-p384-complete.json',
    'gf16-parallel-r4-d099-q49.json',
    'gf16-r4-d10-assembly-p384-complete.json',
    'gf16-lowthreads-r4-d10-q49.json',
    'gf16-r4-d10-q1-48-prefix.json',
    'gf16-r4-d10-q1-48-prefix-p384.json',
    'shared-gf16-q1-screen.json',
    'shared-gf16-q1-p384.json',
    'shared-gf16-sparse-screen.json',
    'shared-gf16-q16-outward-40.json',
    'shared-gf16-r4-d10-q2-16-fill.json',
    'shared-gf16-performance.log',
    'shared-gf16-profile.log',
)
WITNESS_DIRS = (
    'two-bit-covers', 'two-bit-density-covers', 'two-bit-bridge-covers',
    'two-bit-bridge-subset-split', 'two-bit-bridge-complete',
    'two-bit-bridge-retuned', 'two-bit-bridge-322-400',
    'shared-goal', 'pairwise-goal', 'shared-relaxed',
)


def digest_stream(stream):
    result = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1 << 20), b''):
        result.update(chunk)
    return result.hexdigest()


def digest(path):
    with path.open('rb') as stream:
        return digest_stream(stream)


def selected_files():
    command = ['rg', '--files', 'research', 'spin']
    for pattern in ('*.py', '*.md', '*.tex', '*.bib', '*.cpp', '*.h', '*.hpp',
                    '*.sh', '*.ps1', '*.cmake', 'CMakeLists.txt', 'requirements*.txt'):
        command.extend(('-g', pattern))
    paths = {ROOT / p for p in subprocess.check_output(command, cwd=ROOT, text=True).splitlines()}
    paths.update(ROOT / 'tmp' / name for name in RECEIPTS)
    for name in WITNESS_DIRS:
        directory = ROOT / 'tmp' / name
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        for pattern in ('*.json','*.log','*.txt','*.csv'):
            paths.update(directory.rglob(pattern))
    paths.update((ROOT/'tmp').glob('pairwise-goal-*.log'))

    # Authenticate the existing BCH premises, without invoking the driver's
    # main() or rewriting its historical PROGRESS_VERIFIED receipt.
    sys.path.insert(0, str(ROOT / 'research/workstreams/inner_design/asymmetric/bch256'))
    import verify_progress
    dependencies = {**verify_progress.outer_caps(), **verify_progress.model.sources()}
    for name, expected in dependencies.items():
        path = ROOT / 'research' / name
        if digest(path) != expected:
            raise ValueError(f'Changed authenticated dependency: {path}')
        paths.add(path)
    paths.add(ROOT / 'research/bch_spectrum_work/bch_spectrum_codex_bundle/generated/bch256_shift_rank_q30_refined.json')
    for path in paths:
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise ValueError(f'Missing or nonlocal snapshot input: {path}')
    return sorted(paths)


def verify(path):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError('Duplicate archive members')
        manifest = json.loads(archive.read('MANIFEST.json'))
        if manifest.get('schema') != 'spin-proof-archive-1':
            raise ValueError('Unknown archive schema')
        records = manifest['files']
        if set(names) != {'MANIFEST.json', *records}:
            raise ValueError('Archive membership differs from manifest')
        for name, record in records.items():
            relative = PurePosixPath(name)
            if relative.is_absolute() or '..' in relative.parts or '\\' in name or ':' in name:
                raise ValueError('Unsafe member path')
            if archive.getinfo(name).file_size != record['bytes']:
                raise ValueError(f'Size mismatch: {name}')
            with archive.open(name) as stream:
                if digest_stream(stream) != record['sha256']:
                    raise ValueError(f'Hash mismatch: {name}')
    return manifest


def create(path):
    if path.exists():
        raise FileExistsError(f'Will not overwrite snapshot: {path}')
    paths = selected_files()
    records = {p.relative_to(ROOT).as_posix(): dict(bytes=p.stat().st_size, sha256=digest(p)) for p in paths}
    manifest = dict(
        schema='spin-proof-archive-1', created_utc=datetime.now(timezone.utc).isoformat(),
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        includes_uncommitted_sources=True, python=platform.python_version(),
        packages={p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'python-flint')},
        scope='Research/library sources, authenticated BCH premises, selected completed and partial proof witnesses. No numerical replay implied.',
        files=records)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for file in paths:
            archive.write(file, file.relative_to(ROOT).as_posix())
        archive.writestr('MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    verify(path)  # Also detects source changes while the archive was written.
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--create', type=Path)
    mode.add_argument('--verify', type=Path)
    args = parser.parse_args()
    path = args.create or args.verify
    manifest = create(path) if args.create else verify(path)
    print('Verified', len(manifest['files']), 'files;', path.stat().st_size, 'archive bytes')
    print('SHA256', digest(path))
    print(path.resolve())


if __name__ == '__main__':
    main()
