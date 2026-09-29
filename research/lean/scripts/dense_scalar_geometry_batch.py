"""Associate checked scalar boxes with global geometry while numeric replay proceeds."""
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import subprocess
import sys
import time
from dense_scalar_geometry import generate, ROOT

DATA=ROOT/'scripts/map_data'
REPORT=DATA/'dense_scalar_geometry_verification.json'
ALLOWED={'propext','Classical.choice','Quot.sound'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check(i):
    tag=f'B{i:03d}'
    source=ROOT/f'SpinCodes/Structured/DenseScalarExact{tag}.lean'
    obj=ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/DenseScalarExact{tag}.olean'
    log=DATA/f'peach_DenseScalarExact{tag}.log'
    # Objects are fetched only after successful remote compilation.
    while not (source.exists() and obj.exists() and log.exists() and
               obj.stat().st_mtime>=source.stat().st_mtime and
               ': error:' not in log.read_text(encoding='utf-8')):
        time.sleep(5)
    path=generate(i)
    batchlog=DATA/f'dense_scalar_geometry_batch_{tag}.log'
    with batchlog.open('w',encoding='utf-8') as output:
        result=subprocess.run([sys.executable,'scripts/peach-lean.py',path.relative_to(ROOT).as_posix()],
                              cwd=ROOT,stdout=output,stderr=subprocess.STDOUT)
    if result.returncode:
        return tag,{'status':'FAILED','log':batchlog.name}
    obj=ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/DenseScalarGeometry{tag}.olean'
    log=DATA/f'peach_DenseScalarGeometry{tag}.log'
    audits=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",log.read_text(encoding='utf-8'))
    assert len(audits)==1 and set(a.strip() for a in audits[0][1].split(','))<=ALLOWED, tag
    return tag,{'status':'PASS','source_sha256':sha(path),'object_sha256':sha(obj),
               'axiom_audits':1,'log':log.name}


def main():
    report=json.loads(REPORT.read_text()) if REPORT.exists() else {
        'scope':'283 scalar boxes associated with exact indexed global geometry and actual outer supports. Other numerical families remain separate.',
        'boxes':{}}
    report['status']='RUNNING'
    report['concurrency']=3
    todo=[]
    for i in range(283):
        tag=f'B{i:03d}'
        old=report['boxes'].get(tag,{})
        path=ROOT/f'SpinCodes/Structured/DenseScalarGeometry{tag}.lean'
        obj=ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/DenseScalarGeometry{tag}.olean'
        if (old.get('status')=='PASS' and path.exists() and obj.exists()
            and sha(path)==old['source_sha256'] and sha(obj)==old['object_sha256']):
            continue
        todo.append(i)
    def save():
        temp=REPORT.with_suffix('.tmp')
        temp.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        temp.replace(REPORT)
    save()
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs=[pool.submit(check,i) for i in todo]
        for job in as_completed(jobs):
            tag,result=job.result()
            report['boxes'][tag]=result
            save()
            print(tag,'geometry',result['status'],flush=True)
    report['checked_boxes']=sum(x['status']=='PASS' for x in report['boxes'].values())
    report['status']='PASS' if report['checked_boxes']==283 else 'FAILED'
    save()

if __name__=='__main__':
    main()
