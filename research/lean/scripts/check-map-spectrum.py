"""Kernel replay of concrete map spectrum blocks with bounded workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--basis', action='store_true')
    parser.add_argument('--resume', action='store_true',
                        help='Reuse successful checks only with matching source, output, and dependency hashes.')
    parser.add_argument('--min-block', type=int, default=0)
    parser.add_argument('--max-block', type=int, default=511)
    args = parser.parse_args()
    assert 1 <= args.jobs <= 4
    if args.basis:
        records = [dict(block=-1, **json.loads((DATA/'basis_manifest.json').read_text()))]
        report_path = DATA/'kernel_basis.json'
    else:
        records = [r for r in json.loads((DATA/'manifest.json').read_text())
                   if args.min_block <= r['block'] <= args.max_block]
        suffix = '' if (args.min_block,args.max_block) == (0,511) else f'_{args.min_block}_{args.max_block}'
        report_path = DATA/f'kernel_blocks{suffix}.json'
    assert records
    logdir = DATA/'kernel_logs'
    logdir.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
    dependencies = ['PackedMapDefs', 'ConcreteMapData', 'MapSpectrumDefs']
    if not args.basis:
        dependencies += ['MapSpectrumData/Basis']
    def dependency_hashes():
        return {name: sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')
                for name in dependencies}
    initial = dependency_hashes()
    results = []
    if args.resume and report_path.exists():
        previous = json.loads(report_path.read_text())
        assert previous['dependency_olean_sha256'] == initial, 'Dependencies changed; start a fresh replay.'
        by_block = {r['block']: r for r in records}
        for result in previous['blocks']:
            record = by_block.get(result['block'])
            if record is None or result['exit_code'] != 0:
                continue
            source = ROOT/record['path']
            output = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
            if (result['source_sha256'] == record['sha256'] == sha(source)
                    and output.exists() and result.get('output_sha256') == sha(output)):
                results.append(result)
    completed = {r['block'] for r in results}
    def save(status):
        report_path.write_text(json.dumps(dict(status=status, blocks=results,
            dependency_olean_sha256=initial),indent=2)+'\n')
    save('RUNNING')
    def check(record):
        source = ROOT/record['path']
        digest = sha(source)
        assert digest == record['sha256']
        output = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
        output.parent.mkdir(parents=True,exist_ok=True)
        start = time.monotonic()
        proc = subprocess.Popen(['lean','-o',str(output),str(source)],cwd=ROOT,env=env,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding='utf-8')
        try:
            log,_ = proc.communicate(timeout=900)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            if os.name == 'nt':
                subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],
                    stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            else:
                proc.kill()
            proc.communicate()
            code,log = 124,'Kernel check timed out after 900 seconds.\n'
        logfile = f'block{record["block"]}.log'
        (logdir/logfile).write_text(log,encoding='utf-8')
        assert sha(source) == digest
        return dict(block=record['block'],path=record['path'],source_sha256=digest,
            output_sha256=sha(output) if code == 0 else None,
            exit_code=code,seconds=round(time.monotonic()-start,2),log=logfile)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = [pool.submit(check,r) for r in records if r['block'] not in completed]
        for future in as_completed(pending):
            result = future.result()
            results.append(result)
            save('RUNNING')
            print(f'{len(results)}/{len(records)} block={result["block"]} '
                  f'exit={result["exit_code"]} seconds={result["seconds"]}',flush=True)
    assert dependency_hashes() == initial, 'A dependency changed during replay.'
    status = 'PASS' if all(r['exit_code'] == 0 for r in results) else 'FAIL'
    save(status)
    if status != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
