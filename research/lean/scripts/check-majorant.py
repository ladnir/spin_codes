"""Kernel-check generated majorant modules with bounded build concurrency.

Every requested part is rebuilt from its Lean source. No result is accepted on
the basis of the generator's status or a previously cached output file.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import hashlib
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',type=int,default=6)
    parser.add_argument('--segment',type=int)
    parser.add_argument('--phase',choices=['parts','assembly'],default='parts')
    args=parser.parse_args()
    records=json.loads((ROOT/'scripts/majorant_data/manifest.json').read_text())
    if args.segment is not None:
        records=[r for r in records if r['segment']==args.segment]
    modules=[p['module'] for r in records for p in r['parts']]
    if args.phase == 'assembly':
        modules=[f'SpinCodes.Majorant.Segment{r["segment"]}' for r in records]
    logdir=ROOT/'scripts/majorant_data/check_logs'
    logdir.mkdir(exist_ok=True)
    env=os.environ.copy()
    env['LEAN_PATH']=os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean',*[p/'.lake/build/lib/lean' for p in (ROOT/'.lake/packages').iterdir() if p.is_dir()]])

    def check(module):
        relative=Path(*module.split('.'))
        source=ROOT/relative.with_suffix('.lean')
        output=(ROOT/'.lake/build/lib/lean'/relative).with_suffix('.olean')
        output.parent.mkdir(parents=True,exist_ok=True)
        start=time.monotonic()
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        run=subprocess.run(['lean','-o',str(output),str(source)],cwd=ROOT,env=env,
                           stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding='utf-8')
        (logdir/(module+'.log')).write_text(run.stdout,encoding='utf-8')
        if hashlib.sha256(source.read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'Source changed during kernel check: {source}')
        return dict(module=module,exit_code=run.returncode,seconds=time.monotonic()-start,
                    source_sha256=digest)

    results=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending={pool.submit(check,m):m for m in modules}
        for future in as_completed(pending):
            result=future.result()
            results.append(result)
            print(len(results),'/',len(modules),json.dumps(result),flush=True)
    report=dict(status='PASS' if all(r['exit_code']==0 for r in results) else 'FAIL',
                parts=results)
    (ROOT/f'scripts/majorant_data/kernel_{args.phase}.json').write_text(json.dumps(report,indent=2)+'\n')
    if report['status']!='PASS':
        raise SystemExit(1)


if __name__=='__main__':
    main()
