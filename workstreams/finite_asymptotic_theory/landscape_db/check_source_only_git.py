"""Reject generated landscape outputs in the tip or unpublished history."""
import argparse
from pathlib import PurePosixPath
import subprocess

PREFIX = 'workstreams/finite_asymptotic_theory/landscape_db/'
SOURCE_SUFFIXES = {'.py', '.cpp', '.ps1', '.sql', '.md'}


def source_path(name):
    path = PurePosixPath(name)
    return (path.suffix in SOURCE_SUFFIXES
            or path.name in ('.gitignore', '.gitattributes', 'catalog.json'))


def git(*args, input=None):
    return subprocess.check_output(['git', *args], input=input)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', default='origin/main')
    parser.add_argument('--head', default='HEAD')
    args = parser.parse_args()
    failures = []
    for name in git('ls-tree', '--full-tree', '-r', '--name-only', args.head, '--', PREFIX).decode().splitlines():
        if not source_path(name):
            failures.append('generated file in tree: '+name)
    objects = git('rev-list', '--objects', args.head, '--not', args.base)
    descriptions = git('cat-file', '--batch-check=%(objectname) %(objecttype) %(objectsize) %(rest)', input=objects)
    blob_count = total = maximum = 0
    for line in descriptions.decode().splitlines():
        fields = line.split(' ', 3)
        if len(fields) < 3 or fields[1] != 'blob':
            continue
        size = int(fields[2]); name = fields[3] if len(fields) > 3 else ''
        blob_count += 1; total += size; maximum = max(maximum, size)
        if size > 1024*1024:
            failures.append('new blob exceeds 1 MiB: '+name)
        if name.startswith(PREFIX) and not source_path(name):
            failures.append('generated file in unpublished history: '+name)
    if failures:
        raise SystemExit('\n'.join(failures))
    print(f'Source-only check passed: {blob_count} new blobs, {total} bytes total, largest {maximum} bytes.')


if __name__ == '__main__':
    main()
