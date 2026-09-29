"""Check scalar local boxes independently, keeping a resumable hash manifest.

These are proof replays, not benchmarks. Local box validity does not certify
the outer segment interpretation, coverage, or the full first-moment sum.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
REPORT = DATA/'dense_scalar_exact_boxes_verification.json'
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(index):
    result = []
    for suffix in ['Data', '']:
        name = f'DenseScalarExactB{index:03d}{suffix}'
        src = ROOT/f'SpinCodes/Structured/{name}.lean'
        obj = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
        log = DATA/f'peach_{name}.log'
        assert src.exists() and obj.exists() and log.exists(), name
        assert ': error:' not in log.read_text(encoding='utf-8'), name
        result.append({'name':name,'source_sha256':sha(src),'object_sha256':sha(obj)})
    audit = (DATA/f'peach_DenseScalarExactB{index:03d}.log').read_text(encoding='utf-8')
    rows = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", audit)
    assert len(rows) == 2, index
    assert all(set(a.strip() for a in ax.split(',')) <= ALLOWED for _, ax in rows), index
    return result


def check(index):
    started = time.monotonic()
    tag = f'B{index:03d}'
    log = DATA/f'dense_scalar_exact_batch_{tag}.log'
    with log.open('w',encoding='utf-8') as output:
        steps = [[sys.executable,'scripts/dense_scalar_exact.py',str(index)]]
        steps += [[sys.executable,'scripts/peach-lean.py',f'SpinCodes/Structured/DenseScalarExact{tag}{s}.lean']
                  for s in ['Data','']]
        for command in steps:
            result = subprocess.run(command,cwd=ROOT,stdout=output,stderr=subprocess.STDOUT)
            output.flush()
            if result.returncode:
                return tag, {'status':'FAILED','command':command,'log':log.name}
    return tag, {'status':'PASS','elapsed_seconds':time.monotonic()-started,
                 'modules':records(index),'axiom_audits':2,'log':log.name}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start',type=int,default=1)
    parser.add_argument('--stop',type=int,default=283)
    parser.add_argument('--workers',type=int,default=2)
    args = parser.parse_args()
    report = json.loads(REPORT.read_text()) if REPORT.exists() else {
        'status':'RUNNING','scope':'Exact scalar transfer and local four-vertex exponent boxes only. '
        'Full-domain coverage, outer affine semantics, and summed first moment are separate.',
        'concurrency':args.workers,'boxes':{},
    }
    if 'B000' not in report['boxes']:
        report['boxes']['B000']={'status':'PASS','modules':records(0),'axiom_audits':2,'source':'Initial individually checked prototype'}
    todo = []
    for index in range(args.start,args.stop):
        old = report['boxes'].get(f'B{index:03d}',{})
        if old.get('status') == 'PASS' and old['modules'] == records(index):
            continue
        todo.append(index)
    def save():
        tmp = REPORT.with_suffix('.tmp')
        tmp.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        tmp.replace(REPORT)
    save()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(check,index) for index in todo]
        for future in as_completed(futures):
            tag, result = future.result()
            report['boxes'][tag] = result
            save()
            print(tag,result['status'],flush=True)
    report['status'] = 'PASS' if all(x['status']=='PASS' for x in report['boxes'].values()) else 'FAILED'
    report['checked_boxes'] = sum(x['status']=='PASS' for x in report['boxes'].values())
    save()
    print(report['status'],report['checked_boxes'],'local scalar boxes',flush=True)


if __name__ == '__main__':
    main()
