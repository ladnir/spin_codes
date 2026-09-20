"""Copy frozen BCH sources and local audit inputs without importing Git history.

The manifest is compact and tracked; bulk runtime inputs remain ignored.
Copies are byte-preserving and never overwrite different existing content.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.set_int_max_str_digits(0)
BRIDGE = Path('workstreams/bch_rm2sub_bridge')
BCH = Path('bch_spectrum_work/bch_spectrum_codex_bundle')
MANIFEST = BRIDGE / 'MIGRATION_MANIFEST.json'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def contained(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f'Path leaves workspace: {relative}')
    return path


def discover(root):
    bridge, bch = root / BRIDGE, root / BCH
    seeds = [bridge / 'generated/larger_t64_s20_full_closure_audit.json',
             bridge / 'generated/christoffel_oa29_caps.json']
    seeds += list((bridge / 'generated').glob('joint_shell_*/cap.json'))
    seeds += list((bridge / 'inputs').glob('*.json'))
    seeds += list((bridge / 'generated/larger_state_inputs_v1').glob('*.json'))
    seen, hashes, external = set(), {}, []
    aliases = {'diagnostic': bch / 'generated/random_inner_oa15_majorant_diagnostic.json',
               'johnson_w52': bch / 'generated/johnson_n256_w52_d38.json',
               'johnson_w54': bch / 'generated/johnson_n256_w54_d38.json'}
    while seeds:
        path = seeds.pop().resolve()
        if path in seen:
            continue
        if not path.is_file() or not path.is_relative_to(root):
            raise ValueError(f'Missing or external input: {path}')
        seen.add(path)
        if path.suffix != '.json':
            continue
        obj = json.loads(path.read_text(encoding='utf-8-sig'))
        scope = bridge if path.is_relative_to(bridge) else bch

        def visit(value):
            if isinstance(value, list):
                for item in value:
                    visit(item)
            if not isinstance(value, dict):
                return
            for field, contents in value.items():
                if field.endswith('sha256') and isinstance(contents, dict):
                    for name, expected in contents.items():
                        if not isinstance(expected, str) or len(expected) != 64:
                            continue
                        relative = Path(name.replace('\\', '/'))
                        # Original source-worktree hashes record provenance;
                        # snapshot_sha256 authenticates the copied runtime inputs.
                        if field == 'source_sha256' and (
                                'source_worktree' in value or 'source_directory' in value):
                            external.append({'receipt': path.relative_to(root).as_posix(),
                                             'path': name, 'sha256': expected})
                            continue
                        if relative.is_absolute():
                            raise ValueError(f'Unexpected absolute dependency: {name}')
                        bases = [path.parent] if field == 'snapshot_sha256' else (
                            [bch] if field == 'outer_sha256' else
                            [scope, scope / 'code', path.parent, root])
                        candidates = [base / relative for base in bases]
                        if path.name == 'random_inner_threshold_outward.json' and name in aliases:
                            candidates.insert(0, aliases[name])
                        source = next((p.resolve() for p in candidates if p.is_file()), None)
                        if source is None or not source.is_relative_to(root):
                            raise ValueError(f'Unresolved dependency: {path}: {field}: {name}')
                        actual = hashes.get(source)
                        if actual is None:
                            actual = hashes[source] = digest(source)
                        if actual != expected:
                            raise ValueError(f'Original receipt hash mismatch: {source}')
                        seeds.append(source)
                elif isinstance(contents, (dict, list)):
                    visit(contents)
        visit(obj)
    return seen, hashes, external


def copy_checked(source, target, expected):
    if digest(source) != expected:
        raise ValueError(f'Source changed: {source}')
    if target.exists():
        if digest(target) != expected:
            raise ValueError(f'Refusing to overwrite different content: {target}')
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    if digest(target) != expected:
        raise ValueError(f'Copy verification failed: {target}')


def migrate(source, target):
    manifest_path = target / MANIFEST
    if manifest_path.exists():
        raise ValueError('Migration already recorded; use --restore or --verify')
    inputs, hashes, external = discover(source)
    sources = {p for p in (source / BRIDGE).iterdir()
               if p.is_file() and (p.suffix in ('.py', '.md') or p.name == '.gitignore')}
    sources |= {p for p in (source / BCH / 'code').iterdir()
                if p.is_file() and p.suffix in ('.py', '.cpp')}
    # Existing BCH prose is already on main. Refuse semantic differences.
    for p in (source / BCH).glob('*.md'):
        dest = target / p.relative_to(source)
        if dest.exists():
            if dest.read_bytes().replace(b'\r\n', b'\n') != p.read_bytes().replace(b'\r\n', b'\n'):
                raise ValueError(f'BCH prose differs from main; review first: {dest}')
        else:
            sources.add(p)
    files = []
    for p in sorted(inputs | sources):
        relative = p.relative_to(source)
        expected = hashes.get(p) or digest(p)
        copy_checked(p, target / relative, expected)
        files.append({'path': relative.as_posix(), 'bytes': p.stat().st_size,
                      'sha256': expected, 'kind': 'source' if p in sources else 'local_audit_input'})
    payload = {'schema': 'bch-workspace-migration-v1', 'github_base': 'cf78d822',
               'legacy_checkpoint': 'c3d8225',
               'legacy_source_note': 'Bridge sources were untracked; per-file hashes identify the migrated content.',
               'files': files, 'external_snapshot_provenance': external}
    with manifest_path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'files': len(files), 'source_files': sum(f['kind'] == 'source' for f in files),
                      'source_bytes': sum(f['bytes'] for f in files if f['kind'] == 'source'),
                      'local_input_bytes': sum(f['bytes'] for f in files if f['kind'] == 'local_audit_input')}))


def restore_or_verify(target, source=None):
    payload = json.loads((target / MANIFEST).read_text(encoding='utf-8'))
    for row in payload['files']:
        dest = contained(target, row['path'])
        if source is not None:
            copy_checked(contained(source, row['path']), dest, row['sha256'])
        if not dest.is_file() or digest(dest) != row['sha256']:
            raise ValueError(f'Missing or changed migration file: {dest}')
    print(f"Verified {len(payload['files'])} migrated files byte-for-byte")


def extend_manifest(source, target, relatives):
    path = target / MANIFEST
    payload = json.loads(path.read_text(encoding='utf-8'))
    rows = {row['path']: row for row in payload['files']}
    for relative in relatives:
        src = contained(source, relative)
        dest = contained(target, relative)
        expected = digest(src)
        copy_checked(src, dest, expected)
        key = src.relative_to(source).as_posix()
        row = dict(path=key, bytes=src.stat().st_size, sha256=expected,
                   kind='source' if src.suffix in ('.py', '.cpp') else 'local_audit_input')
        if key in rows and rows[key] != row:
            raise ValueError(f'Manifest entry changed: {key}')
        rows[key] = row
    payload['files'] = [rows[key] for key in sorted(rows)]
    with path.open('w', encoding='utf-8', newline='\n') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')
    print(f"Recorded {len(relatives)} additional runtime dependencies")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--target-root', type=Path, default=Path(__file__).resolve().parents[2])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--migrate', action='store_true')
    mode.add_argument('--restore', action='store_true')
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--extend', nargs='+', metavar='RELATIVE_PATH')
    args = parser.parse_args()
    target = args.target_root.resolve()
    source = args.source_root.resolve() if args.source_root else None
    if not (target / '.git').exists() or (source and not (source / '.git').exists()):
        parser.error('Expected source/target Git workspaces')
    if args.verify:
        restore_or_verify(target)
    elif source is None or source == target:
        parser.error('Provide a distinct source workspace')
    elif args.restore:
        restore_or_verify(target, source)
    elif args.extend:
        extend_manifest(source, target, args.extend)
    else:
        migrate(source, target)


if __name__ == '__main__':
    main()
