"""Current finite IMT evidence: check, figures, inventory, or package.

No search, interval replay, benchmark, download, or source commit is implicit.
The default bundle covers selected-certificate/timing pins. --include-q1 adds
the engineering grid. Neither scope bundles installed external dependencies.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'paper'))
import imt_results


def accepted_roots(include_q1=False):
    entries = list(imt_results.HALF_PINS.values()) + [imt_results.TIMING_PIN,
        imt_results.QUARTER_PIN, imt_results.QUARTER_PROOF_PIN]
    if include_q1:
        import build_imt_parameter_figures as figures
        import build_imt_length_figure as length
        entries += [figures.GRID_PIN, figures.REPLAY_PIN]
        entries += [length.GRID_PIN, length.REPLAY_PIN]
    return {'workstreams/inner_design/' + name: digest for name,digest in entries}


def checked_path(root, name):
    if Path(name).is_absolute() or '\\' in name:
        raise ValueError(f'Expected repository-relative POSIX path: {name}')
    path = (root/name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes repository: {name}')
    return path


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def add_pin(pins,name,digest):
    if not re.fullmatch('[0-9a-f]{64}',digest):
        raise ValueError(f'Invalid SHA-256: {name}')
    if name in pins and pins[name] != digest:
        raise ValueError(f'Conflicting pin: {name}')
    pins[name] = digest


def inventory(root=ROOT, roots=None):
    root = root.resolve()
    roots = accepted_roots() if roots is None else roots
    pins, unreadable = {},[]
    for name,digest in roots.items():
        add_pin(pins,name,digest)
        path = checked_path(root,name)
        if not path.is_file() or sha(path) != digest:
            unreadable.append(name)
            continue
        record = json.loads(path.read_text())
        for dependency,expected in record['source_sha256'].items():
            checked_path(root,dependency)
            add_pin(pins,dependency,expected)
    try:
        tracked = set(subprocess.check_output(['git','ls-files','-z'],cwd=root,
            stderr=subprocess.DEVNULL).decode().split('\0'))
    except (subprocess.CalledProcessError,FileNotFoundError):
        tracked = None
    files = []
    for name,digest in sorted(pins.items()):
        path = checked_path(root,name)
        exists = path.is_file()
        files.append(dict(path=name,sha256=digest,bytes=path.stat().st_size if exists else None,
            tracked=None if tracked is None else name in tracked,
            status='missing' if not exists else 'ok' if sha(path)==digest else 'mismatch'))
    q1_included = any(name.endswith('/PARAMETER_NO_CONSTANT_Q1_v1.json') for name in roots)
    return dict(scope='accepted selected finite IMT certificate and timing pins'
                + (' plus diagnostic Q1 grid' if q1_included else ''),
        root_receipts=len(roots),unreadable_roots=unreadable,inventory_complete=not unreadable,
        entries=len(files),missing=sum(r['status']=='missing' for r in files),
        mismatched=sum(r['status']=='mismatch' for r in files),
        present_bytes=sum(r['bytes'] or 0 for r in files),
        untracked_bytes=sum(r['bytes'] or 0 for r in files if r['tracked'] is False),
        q1_grid_included=q1_included,
        full_parameter_grid_certified=False,numerical_replay=False,
        self_contained_runtime_bundle=False,files=files)


def verify_archive(path, report):
    with zipfile.ZipFile(path) as archive:
        expected = {r['path'] for r in report['files']} | {'imt-evidence-manifest.json'}
        if len(archive.namelist()) != len(expected) or set(archive.namelist()) != expected:
            raise ValueError('Archive members differ from inventory')
        for row in report['files']:
            with archive.open(row['path']) as stream:
                digest = hashlib.file_digest(stream,'sha256').hexdigest()
            if digest != row['sha256']:
                raise ValueError(f"Archived bytes differ: {row['path']}")
        if json.loads(archive.read('imt-evidence-manifest.json')) != report:
            raise ValueError('Archive manifest differs')


def pack(output, root=ROOT, roots=None):
    report = inventory(root,roots)
    if not report['inventory_complete'] or report['missing'] or report['mismatched']:
        raise ValueError('Complete matching selected IMT evidence is required')
    output = output.resolve()
    if output.suffix != '.zip':
        raise ValueError('Evidence output must end in .zip')
    output.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive creation. A failure leaves an explicitly incomplete file;
    # there is no destructive overwrite or automatic cleanup of user data.
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for row in report['files']:
            source = checked_path(root,row['path'])
            with source.open('rb') as reader, archive.open(row['path'],'w',force_zip64=True) as writer:
                shutil.copyfileobj(reader,writer,length=1 << 20)
        archive.writestr('imt-evidence-manifest.json',json.dumps(report,indent=2))
    verify_archive(output,report)
    return dict(status='SELECTED_IMT_EVIDENCE_ARCHIVE_VERIFIED',path=str(output),
                sha256=sha(output),bytes=output.stat().st_size,entries=report['entries'],
                numerical_replay=False,q1_grid_included=report['q1_grid_included'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('check','figures','inventory','pack'))
    p.add_argument('--output',type=Path)
    p.add_argument('--include-q1',action='store_true',help='Add Q1 grid pins to inventory or pack')
    a = p.parse_args()
    if a.command == 'pack' and a.output is None:
        p.error('pack requires --output with a new ZIP path')
    if a.include_q1 and a.command not in ('inventory','pack'):
        p.error('--include-q1 applies only to inventory or pack; check and figures already include Q1')
    if a.command == 'check':
        if a.output is not None: p.error('check does not produce a data file')
        subprocess.run([sys.executable,'-B',str(ROOT/'paper/check_finite_integration.py')],
                       cwd=ROOT,check=True)
    elif a.command == 'figures':
        if a.output is not None: p.error('figures writes the fixed paper inputs')
        subprocess.run([sys.executable,'-B',str(ROOT/'paper/build_imt_parameter_figures.py')],
                       cwd=ROOT,check=True)
    elif a.command == 'pack':
        print(json.dumps(pack(a.output,roots=accepted_roots(a.include_q1)),indent=2))
    else:
        report = inventory(roots=accepted_roots(a.include_q1))
        if a.output:
            a.output.parent.mkdir(parents=True,exist_ok=True)
            with a.output.open('x',encoding='utf-8',newline='\n') as stream:
                json.dump(report,stream,indent=2)
                stream.write('\n')
        print(json.dumps({k:v for k,v in report.items() if k != 'files'},indent=2))
        if not report['inventory_complete'] or report['missing'] or report['mismatched']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
