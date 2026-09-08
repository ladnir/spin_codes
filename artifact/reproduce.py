"""Reader-facing SPIN artifact commands. No search or benchmark is implicit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = Path('workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub')
MANIFEST_SHA = 'cd5afa44c245bc46b6c55376d5f86908c5c0d069d17b1b5669c64a28ed466786'
AUDIT = Path('workstreams/bch_rm2sub_bridge/generated/t128_s19_ladder_completion_audit_v1.json')
AUDIT_SHA = 'dc06c1e9770a33ca37c8ea8824ae554d3b2ce84011477e55b210d39a38446aa9'
OLD_SELECTION = ('constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_'
                 'rm2sub_t128_s20/receipts/min_state/s19_rm2sub_selection.json')


def checked_path(name, root=ROOT):
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes artifact: {name}')
    return path


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def authenticated_json(path, expected):
    if sha(path) != expected:
        raise ValueError(f'Manifest hash mismatch: {path}')
    return json.loads(path.read_text(encoding='utf-8'))


def asymptotic_entries():
    manifest = authenticated_json(ROOT / BUNDLE / 'SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json', MANIFEST_SHA)
    entries = [(BUNDLE / Path(e['path']).name, e['sha256']) for e in manifest['local_files']]
    entries += [(BUNDLE / 'linear_time_audit' / Path(e['path']).name, e['sha256'])
                for e in manifest['companion_linear_time_files']]
    entries += [(Path('artifact/data/asymptotic/s19_rm2sub_selection.json')
                 if e['path'] == OLD_SELECTION else Path(e['path']), e['sha256'])
                for e in manifest['frozen_dependencies']]
    if len(entries) != 31:
        raise ValueError('Expected 31 asymptotic manifest entries')
    return entries


def check_asymptotic():
    entries = asymptotic_entries()
    for name, expected in entries:
        if sha(checked_path(name)) != expected:
            raise ValueError(f'Asymptotic dependency changed: {name}')
    print('Asymptotic manifest: 31 entries authenticated; no interval replay.', flush=True)


def dependency_inventory(root=ROOT):
    audit = authenticated_json(root / AUDIT, AUDIT_SHA)
    pins = dict(audit['source_sha256'])
    pins[AUDIT.as_posix()] = AUDIT_SHA
    perf = json.loads((root / 'workstreams/bare_bch_rm2sub/PERFORMANCE.json').read_text())
    for name, expected in perf['source_sha256'].items():
        if name in pins and pins[name] != expected:
            raise ValueError(f'Conflicting hash: {name}')
        pins[name] = expected
    tracked = set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0'))
    files = []
    for name, expected in sorted(pins.items()):
        path = checked_path(name, root)
        exists = path.is_file()
        files.append(dict(path=name, sha256=expected, tracked=name in tracked,
                          bytes=path.stat().st_size if exists else None,
                          status=('missing' if not exists else 'ok' if sha(path) == expected else 'mismatch')))
    return dict(scope='BCH ladder audit and performance pins; not all research dependencies',
                numerical_replay=False, entries=len(files),
                missing=sum(e['status'] == 'missing' for e in files),
                mismatched=sum(e['status'] == 'mismatch' for e in files),
                untracked=sum(not e['tracked'] for e in files),
                present_bytes=sum(e['bytes'] or 0 for e in files),
                untracked_present_bytes=sum(e['bytes'] or 0 for e in files if not e['tracked']),
                files=files)


def run(*args, cwd=ROOT):
    print('+ ' + ' '.join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=cwd, check=True)


def quick():
    run(sys.executable, '-B', ROOT / 'paper/build_parameter_figures.py', '--check')
    run(sys.executable, '-B', ROOT / 'paper/check_finite_integration.py')
    check_asymptotic()
    print('QUICK CHECK PASSED (retained evidence, not full numerical replay)', flush=True)


def pack_evidence(output):
    result = dependency_inventory()
    if result['missing'] or result['mismatched']:
        raise ValueError('Complete matching evidence is required before packaging')
    output = output.resolve()
    if output.suffix != '.zip':
        raise ValueError('Evidence output must end in .zip')
    # Exclusive creation: never overwrite a prior release or a source file.
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for entry in result['files']:
            data = checked_path(entry['path']).read_bytes()
            if hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError(f'Input changed while packaging: {entry["path"]}')
            archive.writestr(entry['path'], data)
        archive.writestr('evidence-manifest.json', json.dumps(result, indent=2))
    # Check every archived byte against the same pins, not just ZIP CRCs.
    with zipfile.ZipFile(output) as archive:
        for entry in result['files']:
            if hashlib.sha256(archive.read(entry['path'])).hexdigest() != entry['sha256']:
                raise ValueError(f'Archive verification failed: {entry["path"]}')
    print(json.dumps(dict(status='EVIDENCE_ARCHIVE_VERIFIED', path=str(output),
                          bytes=output.stat().st_size, sha256=sha(output), entries=result['entries'],
                          numerical_replay=False), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['quick', 'figures', 'paper', 'inventory', 'evidence', 'pack-evidence'])
    parser.add_argument('--output', type=Path, help='New ZIP path; only for pack-evidence')
    args = parser.parse_args()
    if (args.command == 'pack-evidence') != (args.output is not None):
        parser.error('--output is required only with pack-evidence')
    if args.command == 'quick':
        quick()
    elif args.command == 'figures':
        run(sys.executable, '-B', ROOT / 'paper/build_parameter_figures.py')
        run(sys.executable, '-B', ROOT / 'paper/build_parameter_figures.py', '--check')
    elif args.command == 'paper':
        quick()
        (ROOT / 'output/pdf').mkdir(parents=True, exist_ok=True)
        run('latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error',
            '-outdir=../output/pdf', '-jobname=spin_codes_draft', 'main.tex', cwd=ROOT / 'paper')
    elif args.command == 'inventory':
        result = dependency_inventory()
        print(json.dumps(result, indent=2))
        if result['mismatched']:
            return 1
    elif args.command == 'evidence':
        run(sys.executable, '-B', ROOT / 'workstreams/bch_rm2sub_bridge/verify_paper_milestone.py',
            '--require-local-evidence')
    elif args.command == 'pack-evidence':
        pack_evidence(args.output)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f'Artifact check failed: {error}', file=sys.stderr)
        raise SystemExit(1)
