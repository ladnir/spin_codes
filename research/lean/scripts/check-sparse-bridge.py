"""Kernel-check all fixed-weight polynomial identities with bounded workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--weight', type=int)
    parser.add_argument('--min-weight', type=int, default=0)
    parser.add_argument('--max-weight', type=int, default=128)
    parser.add_argument('--phase', choices=['weights', 'dominance', 'contributions'], default='weights')
    args = parser.parse_args()
    manifest = {'weights': 'manifest.json', 'dominance': 'dominance_manifest.json',
                'contributions': 'contribution_manifest.json'}[args.phase]
    records = json.loads((ROOT/'scripts/sparse_data'/manifest).read_text())
    if args.weight is not None:
        records = [r for r in records if r['weight'] == args.weight]
    else:
        records = [r for r in records if args.min_weight <= r['weight'] <= args.max_weight]
    assert records and 1 <= args.jobs <= 6
    logdir = ROOT/'scripts/sparse_data/check_logs'
    logdir.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean', *[p/'.lake/build/lib/lean' for p in
         (ROOT/'.lake/packages').iterdir() if p.is_dir()]])
    dependencies = ['PolyIdentityDefs', 'SparsePowers', 'SparsePolynomialDefs', 'SparseModelData']
    if args.phase == 'dominance':
        dependencies += ['PolyCertDefs', 'PolySignIdentityDefs']
    if args.phase == 'contributions':
        dependencies += ['PolyPackedDefs'] + [f'SparseBridge/Weight{j}' for j in range(129)]
    def dependency_hashes():
        return {name: hashlib.sha256((ROOT/'.lake/build/lib/lean/SpinCodes/Structured'/
                f'{name}.olean').read_bytes()).hexdigest() for name in dependencies}
    initial_dependencies = dependency_hashes()

    def check(record):
        source = ROOT / record['path']
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        assert digest == record['sha256']
        output = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
        output.parent.mkdir(exist_ok=True, parents=True)
        start = time.monotonic()
        proc = subprocess.Popen(['lean', '-o', str(output), str(source)], cwd=ROOT, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')
        try:
            log, _ = proc.communicate(timeout=900)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            if os.name == 'nt':
                subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                proc.kill()
            proc.communicate()
            code, log = 124, 'Kernel check timed out after 900 seconds.\n'
        (logdir/f'{args.phase}{record["weight"]}.log').write_text(log, encoding='utf-8')
        assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
        return dict(weight=record['weight'], exit_code=code,
                    seconds=time.monotonic()-start, source_sha256=digest)

    results = []
    suffix = '' if (args.min_weight, args.max_weight) == (0,128) else f'_{args.min_weight}_{args.max_weight}'
    report = ROOT / (f'scripts/sparse_data/kernel_{args.phase}{suffix}.json' if args.weight is None
                     else f'scripts/sparse_data/kernel_{args.phase}{args.weight}.json')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = [pool.submit(check, r) for r in records]
        for future in as_completed(pending):
            result = future.result()
            results.append(result)
            print(f'{len(results)}/{len(records)} weight={result["weight"]} '
                  f'exit={result["exit_code"]} seconds={result["seconds"]:.2f}', flush=True)
            report.write_text(json.dumps(dict(status='RUNNING', weights=results), indent=2)+'\n')
    stable = dependency_hashes() == initial_dependencies
    status = 'PASS' if stable and all(r['exit_code'] == 0 for r in results) else 'FAIL'
    report.write_text(json.dumps(dict(status=status, weights=results,
                                     stable_dependencies=stable,
                                     dependency_olean_sha256=initial_dependencies), indent=2)+'\n')
    if status != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
