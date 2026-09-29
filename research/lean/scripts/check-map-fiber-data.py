"""Sequential kernel replay of the 129 Fourier rows and numerical cap checks."""
import argparse
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
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    records = json.loads((DATA/'fiber_manifest.json').read_text())
    assert len(records) == 129 and {r['weight'] for r in records} == set(range(129))
    target = DATA/'kernel_fiber_data.json'
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
    names = ['PolyIdentityDefs','FiberNumericsDefs','FiberNumericsData/Tables']
    def dependencies():
        return {name: sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean') for name in names}
    initial = dependencies()
    results = []
    if args.resume and target.exists():
        previous = json.loads(target.read_text())
        assert previous['dependency_olean_sha256'] == initial, 'Dependency changed.'
        by_weight = {r['weight']:r for r in records}
        for result in previous['weights']:
            record = by_weight[result['weight']]
            output = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
            if (result['exit_code'] == 0 and output.exists()
                    and result['source_sha256'] == record['sha256'] == sha(ROOT/record['path'])
                    and result['output_sha256'] == sha(output)):
                results.append(result)
    done = {r['weight'] for r in results}
    def save(status):
        target.write_text(json.dumps(dict(status=status,
            scope='Numerical Fourier rows and cap inequalities only. Concrete spectrum identification remains separate.',
            weights=results,dependency_olean_sha256=initial),indent=2)+'\n')
    save('RUNNING')
    for record in records:
        if record['weight'] in done:
            continue
        source = ROOT/record['path']
        digest = sha(source)
        assert digest == record['sha256']
        output = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
        output.parent.mkdir(parents=True,exist_ok=True)
        logfile = DATA/f'fiber_numeric_{record["weight"]}.log'
        start = time.monotonic()
        with logfile.open('w',encoding='utf-8') as stream:
            proc = subprocess.Popen(['lean','-o',str(output),str(source)],cwd=ROOT,env=env,
                stdout=stream,stderr=subprocess.STDOUT)
            try:
                code = proc.wait(timeout=900)
            except subprocess.TimeoutExpired:
                if os.name == 'nt':
                    subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],
                        stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                else:
                    proc.kill()
                proc.wait()
                code = 124
        assert digest == sha(source)
        results.append(dict(weight=record['weight'],path=record['path'],method=record['method'],
            source_sha256=digest,output_sha256=sha(output) if code == 0 else None,
            exit_code=code,seconds=round(time.monotonic()-start,2),log=logfile.name))
        save('RUNNING' if code == 0 else 'FAIL')
        print(f'{len(results)}/129 weight={record["weight"]} exit={code} '
              f'seconds={results[-1]["seconds"]}',flush=True)
        if code:
            print(logfile.read_text(encoding='utf-8'))
            raise SystemExit(code)
    assert dependencies() == initial, 'Dependency changed during replay.'
    save('PASS')
    print('PASS: 129 Fourier rows and numerical cap checks.',flush=True)


if __name__ == '__main__':
    main()
