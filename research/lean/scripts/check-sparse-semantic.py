"""Kernel-check the matrix bridge in dependency order, recording source hashes."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/sparse_data'
MODULES = ['SparseCancellation', 'SparseColumnBounds', 'SparseNonneg',
           'SparseContraction', 'SparseIteration', 'SparseBridge/Pin']


def main():
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT / '.lake/build/lib/lean'] +
        sorted((ROOT / '.lake/packages').glob('*/.lake/build/lib/lean')))
    lean = shutil.which('lean')
    if lean is None:
        raise RuntimeError('lean is not on PATH')
    report = {'status': 'RUNNING', 'modules': []}
    destination = DATA / 'semantic_verification.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    for module in MODULES:
        source = Path('SpinCodes/Structured') / (module + '.lean')
        digest = hashlib.sha256((ROOT/source).read_bytes()).hexdigest()
        output = ROOT / '.lake/build/lib/lean' / source.with_suffix('.olean')
        log = DATA / ('semantic_' + module.replace('/', '_') + '.log')
        start = time.monotonic()
        with log.open('w', encoding='utf-8') as stream:
            result = subprocess.run([lean, '-o', str(output), str(source)],
                                    cwd=ROOT, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=900)
        assert digest == hashlib.sha256((ROOT/source).read_bytes()).hexdigest(), source
        report['modules'].append(dict(path=source.as_posix(), source_sha256=digest,
            exit_code=result.returncode, elapsed_seconds=round(time.monotonic()-start, 2),
            log=log.name))
        if result.returncode:
            report['status'] = 'FAIL'
        destination.write_text(json.dumps(report, indent=2) + '\n')
        print(f'{module}: {result.returncode}', flush=True)
        if result.returncode:
            print(log.read_text(encoding='utf-8'))
            raise SystemExit(result.returncode)
    report['status'] = 'PASS'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: matrix bridge semantic replay.', flush=True)


if __name__ == '__main__':
    main()
