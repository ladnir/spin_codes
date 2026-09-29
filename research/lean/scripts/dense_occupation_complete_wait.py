"""Read-only final occupation evidence audit; performs no Lean compilation."""
from pathlib import Path
import hashlib, json, re, time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
OUT = DATA / 'dense_occupation_complete_verification.json'
FILES = ['dense_occupation_fixed_verification.json',
         'dense_occupation_fixed_batch_verification.json',
         'dense_occupation_fixed_remaining_verification.json',
         'dense_occupation_fixed_all_boxes_verification.json',
         'dense_occupation_fixed_rates_verification.json',
         'dense_occupation_aggregate_verification.json']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(record):
    tmp = OUT.with_suffix('.audit.tmp')
    tmp.write_text(json.dumps(record, indent=2) + '\n')
    tmp.replace(OUT)

def main():
    script_sha256=sha(Path(__file__))
    save({'status': 'WAITING', 'required_witnesses': 133, 'required_local_boxes': 433,
          'required_indexed_boxes': 433})
    while True:
        try:
            producer_bytes = {file: (DATA / file).read_bytes() for file in FILES}
            producer_hashes = {file: hashlib.sha256(content).hexdigest() for file, content in producer_bytes.items()}
            reports = [json.loads(producer_bytes[file]) for file in FILES]
        except (FileNotFoundError, json.JSONDecodeError):
            time.sleep(5)
            continue
        if any(report.get('status') in ('FAIL', 'FAILED', 'BLOCKED') for report in reports):
            save({'status': 'BLOCKED', 'reason': 'An occupation producer did not pass.'})
            raise SystemExit(1)
        if all(report.get('status') == 'PASS' for report in reports):
            break
        time.sleep(20)

    first, batch, remaining, boxes, rates, aggregate = reports
    assert set(batch['witnesses']) == {f'W{i:03d}' for i in range(1, 9)}
    assert set(remaining['witnesses']) == {f'W{i:03d}' for i in range(9, 133)}
    assert set(boxes['boxes']) == set(rates['boxes']) == {f'B{i:03d}' for i in range(433)}
    modules = {}
    audits = {}

    def add_module(rec):
        name = rec['name']
        assert name not in modules, f'Duplicate module {name}'
        source = ROOT / f'SpinCodes/Structured/{name}.lean'
        obj = ROOT / f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
        assert sha(source) == rec['source_sha256'], f'Source hash mismatch {name}'
        assert sha(obj) == rec['olean_sha256'], f'Object hash mismatch {name}'
        modules[name] = {'source_sha256': rec['source_sha256'], 'olean_sha256': rec['olean_sha256']}

    def add_audits(mapping):
        for name, axioms in mapping.items():
            assert name not in audits
            assert set(x.strip() for x in axioms.split(',')) <= {'propext', 'Classical.choice', 'Quot.sound'}
            audits[name] = axioms

    for module in first['modules']:
        if module['name'] in ('DenseOccupationFixedFirstData', 'DenseOccupationFixedFirst'):
            add_module(module)
    add_audits({name: axioms for name, axioms in first['axioms'].items() if '.First.' in name})
    for report, key in [(batch, 'witnesses'), (remaining, 'witnesses'), (boxes, 'boxes')]:
        for rec in report[key].values():
            assert rec['status'] == 'PASS' and len(rec['modules']) == 2
            for module in rec['modules']:
                add_module(module)
            add_audits(rec['axioms'])
    for tag, rec in rates['boxes'].items():
        assert rec['status'] == 'PASS'
        add_module(dict(rec, name='DenseOccupationFixed' + tag + 'Rate'))
        add_audits(rec['axioms'])
    assert len(modules) == 1565 and len(audits) == 2131
    add_module(dict(aggregate, name='DenseOccupationAllCertified'))
    add_audits(aggregate['axioms'])
    expected_modules = ({'DenseOccupationFixed'+tag+suffix
        for tag in ['First']+[f'W{i:03d}' for i in range(1,133)] for suffix in ['Data','']} |
        {'DenseOccupationFixed'+f'B{i:03d}'+suffix for i in range(433) for suffix in ['Data','','Rate']} |
        {'DenseOccupationAllCertified'})
    expected_audits = ({'Spin.Structured.DenseOccupationFixed.'+tag+'.'+theorem
        for tag in ['First']+[f'W{i:03d}' for i in range(1,133)] for theorem in ['collatz','iterate','routed_probability']} |
        {'Spin.Structured.DenseOccupationFixed.'+f'B{i:03d}'+'.'+theorem
        for i in range(433) for theorem in ['exponent_bound','collatz','certified','certified_uniform']} |
        {'Spin.Structured.DenseGeometry.occupation_certified'})
    assert set(modules)==expected_modules and len(modules)==1566
    assert set(audits)==expected_audits and len(audits)==2132
    # Recheck the immutable evidence snapshot immediately before publishing PASS.
    for file, digest in producer_hashes.items():
        assert sha(DATA/file)==digest, f'Producer changed during audit: {file}'
    for name, rec in modules.items():
        assert sha(ROOT/f'SpinCodes/Structured/{name}.lean')==rec['source_sha256']
        assert sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')==rec['olean_sha256']
    assert sha(Path(__file__))==script_sha256, 'Audit script changed while running'
    save({'status': 'PASS', 'witnesses': 133, 'local_boxes': 433, 'indexed_boxes': 433,
          'module_count': len(modules), 'audit_count': len(audits),
          'scope': 'Current source/object hashes and standard-only audits for every occupation numerical module and its aggregate. Checked cached dependencies reused; this is not a fresh full dependency replay.',
          'modules': modules, 'axioms': audits,
          'producer_manifests': producer_hashes, 'script_sha256': script_sha256,
          'verified_at_utc': datetime.now(timezone.utc).isoformat()})
    print('PASS occupation complete: 133 witnesses, 433 local, 433 indexed, aggregate', flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        save({'status':'FAIL','error':type(exc).__name__+': '+str(exc),'script_sha256':sha(Path(__file__)),'failed_at_utc':datetime.now(timezone.utc).isoformat()})
        raise
