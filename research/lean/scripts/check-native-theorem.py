"""Compile and audit the final native theorem, reusing checked dependencies.

--wait waits for the independently owned dense aggregate and fixed-range objects.
This is an assembly check, not a fresh replay of every numerical proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'scripts/map_data'
REPORT=DATA/'native_theorem_verification.json'
PIN_HASH='447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca'
ALLOWED={'propext','Classical.choice','Quot.sound'}
TARGETS=['ConcreteNativeTheorem','ConcreteNativeTheoremPin']
REQUIRED=['ConcreteNativeFixedAll','DenseOccupationMixedCertified']
INPUT_REPORTS={'ConcreteNativeFixedAll':'native_fixed_all_verification.json',
               'DenseOccupationMixedCertified':'dense_mixed_aggregate_verification.json'}
PROVENANCE_REPORT='native_dependency_timestamp_verification.json'
PROVENANCE_SOURCES=({'SpinCodes/Structured/LowCancellationData/Weight1Block0.lean'} |
    {f'SpinCodes/Structured/LowCancellationData/Weight2Block{i}.lean' for i in range(62)} |
    {'SpinCodes/Structured/MapSpectrumData/Basis.lean', 'SpinCodes/Structured/SparseProgramDefs.lean',
     'SpinCodes/Structured/SparseMaxima.lean', 'SpinCodes/Structured/SparseContributionSound.lean',
     'SpinCodes/Majorant/RefinedData.lean'})
PROVENANCE_INPUTS={'scripts/map_data/low_cancellation_verification.json',
                   'scripts/map_data/encoder_mtime_gap_peach_verification.json'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source(stem):
    return ROOT/f'SpinCodes/Structured/{stem}.lean'

def obj(stem):
    return ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'

def ready(stem):
    log=DATA/f'peach_{stem}.log'
    try:
        report=json.loads((DATA/INPUT_REPORTS[stem]).read_text(encoding='utf-8'))
        if report.get('status')!='PASS':
            return False
        if stem=='ConcreteNativeFixedAll':
            entry=report['modules'][source(stem).relative_to(ROOT).as_posix()]
        else:
            entry=report['groups'][stem]
            if entry.get('status')!='PASS':
                return False
        return (log.exists() and ': error:' not in log.read_text(encoding='utf-8')
                and sha(source(stem))==entry['source_sha256']
                and sha(obj(stem))==entry['object_sha256'])
    except (FileNotFoundError,json.JSONDecodeError,KeyError):
        return False

def upstream_failures():
    failed={}
    for stem,name in INPUT_REPORTS.items():
        try:
            report=json.loads((DATA/name).read_text(encoding='utf-8'))
        except (FileNotFoundError,json.JSONDecodeError):
            continue
        if report.get('status') in {'FAIL','FAILED','BLOCKED'}:
            failed[stem]={'report':name,'status':report['status'],
                          'detail':report.get('error',report.get('reason',report.get('tail')))}
    path=DATA/PROVENANCE_REPORT
    try:
        report=json.loads(path.read_bytes())
        if report.get('status') in {'FAIL','FAILED','BLOCKED'}:
            failed['dependency provenance']={'report':path.name,'status':report['status'],
                                             'detail':report.get('error')}
    except (OSError,json.JSONDecodeError):
        pass
    return failed

def provenance_ready(report=None):
    try:
        if report is None:
            report=json.loads((DATA/PROVENANCE_REPORT).read_bytes())
        if report.get('status')!='PASS' or report.get('checked_modules')!=68:
            return False
        if set(report['modules'])!=PROVENANCE_SOURCES or set(report['input_report_sha256'])!=PROVENANCE_INPUTS:
            return False
        for rel,entry in report['modules'].items():
            if sha(ROOT/rel)!=entry['source_sha256']:
                return False
            if sha((ROOT/'.lake/build/lib/lean'/rel).with_suffix('.olean'))!=entry['object_sha256']:
                return False
        return (len(report['modules'])==68 and
                all(sha(ROOT/rel)==digest for rel,digest in report['input_report_sha256'].items()))
    except (OSError,json.JSONDecodeError,KeyError):
        return False

def closure(module, result):
    if module in result:
        return
    src=ROOT/(module.replace('.','/')+'.lean')
    for dep in re.findall(r'^import (SpinCodes\S*)',src.read_text(encoding='utf-8'),re.M):
        closure(dep,result)
    result[module]=src

def save(report):
    tmp=REPORT.with_suffix('.tmp')
    tmp.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    tmp.replace(REPORT)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--wait',action='store_true')
    args=parser.parse_args()
    sys.stdout.reconfigure(encoding='utf-8')
    report={'status':'WAITING' if args.wait else 'RUNNING',
            'scope':'Final unconditional native distance/rate statement compiled on Peach; '
                    'previously checked project and mathlib dependencies reused. '
                    'Complete project import graph recorded with source/object hashes. '
                    'This does not assert a fresh full numerical source replay.',
            'checks':[]}
    save(report)
    previous=None
    while args.wait and not (all(ready(s) for s in REQUIRED) and provenance_ready()):
        failed=upstream_failures()
        if failed:
            report.update(status='FAIL',upstream_failures=failed)
            save(report)
            raise SystemExit('Upstream verification failed; inspect native_theorem_verification.json')
        missing=[s for s in REQUIRED if not ready(s)]
        if not provenance_ready():
            missing.append('dependency source/object provenance')
        if missing!=previous:
            print('WAITING:',', '.join(missing),flush=True)
            report['waiting_for']=missing
            save(report)
            previous=missing
        time.sleep(20)
    try:
        assert all(ready(s) for s in REQUIRED), 'Required independently checked aggregates missing'
        provenance_bytes=(DATA/PROVENANCE_REPORT).read_bytes()
        provenance=json.loads(provenance_bytes)
        assert provenance_ready(provenance), 'Dependency content provenance missing or changed'
        assert sha(ROOT/'SpinCodes/Pin.lean')==PIN_HASH
        papers={str(p.relative_to(ROOT.parent)):sha(p) for p in sorted((ROOT.parent/'paper').glob('*.tex'))}
        baseline=json.loads((DATA/'parallel_finite_verification.json').read_text())
        assert papers==baseline['paper_sha256'],'Paper changed since prior verified checkpoint'
        graph={}
        closure('SpinCodes.Structured.ConcreteNativeTheoremPin',graph)
        snapshot={m:sha(s) for m,s in graph.items()}
        report.update(status='RUNNING',project_module_count=len(graph))
        report['input_report_sha256']={name:sha(DATA/name) for name in INPUT_REPORTS.values()}
        report['dependency_provenance_report_sha256']=hashlib.sha256(provenance_bytes).hexdigest()
        report.pop('waiting_for',None)
        save(report)
        audits=[]
        for stem in TARGETS:
            start=time.monotonic()
            run=subprocess.run([sys.executable,'scripts/peach-lean.py',str(source(stem).relative_to(ROOT))],
                               cwd=ROOT,capture_output=True,encoding='utf-8',timeout=1200)
            output=run.stdout+run.stderr
            log=DATA/f'native_theorem_check_{stem}.log'
            log.write_text(output,encoding='utf-8')
            report['checks'].append({'module':stem,'exit_code':run.returncode,
                                     'seconds':round(time.monotonic()-start,2),'log':log.name})
            save(report)
            assert run.returncode==0,output[-8000:]
            rows=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",output)
            assert len(rows)==2,(stem,len(rows))
            for name,axioms in rows:
                assert set(a.strip() for a in axioms.split(','))<=ALLOWED,(name,axioms)
                audits.append(name)
            print('PASS',stem,flush=True)
        records={}
        for module,src in graph.items():
            assert sha(src)==snapshot[module],module
            object_path=(ROOT/'.lake/build/lib/lean'/src.relative_to(ROOT)).with_suffix('.olean')
            assert object_path.exists(),module
            if object_path.stat().st_mtime<src.stat().st_mtime:
                entry=provenance['modules'].get(src.relative_to(ROOT).as_posix())
                assert entry is not None,('unverified timestamp inversion',module)
                assert (sha(src),sha(object_path))==(entry['source_sha256'],entry['object_sha256']),module
            records[src.relative_to(ROOT).as_posix()]={
                'source_sha256':sha(src),'object_sha256':sha(object_path)}
        assert sha(ROOT/'SpinCodes/Pin.lean')==PIN_HASH
        assert provenance_ready(), 'Dependency content provenance changed during compilation'
        assert sha(DATA/PROVENANCE_REPORT)==report['dependency_provenance_report_sha256']
        assert papers=={str(p.relative_to(ROOT.parent)):sha(p) for p in sorted((ROOT.parent/'paper').glob('*.tex'))}
        report.update(status='PASS',modules=records,axiom_audits=len(audits),axiom_theorems=audits,
                      original_pin_sha256=PIN_HASH,paper_sha256=papers)
        print('PASS: unconditional actual native minimum-distance limit and exact rate1/2.',flush=True)
    except Exception as exc:
        report.update(status='FAIL',error=str(exc))
        raise
    finally:
        save(report)

if __name__=='__main__':
    main()
