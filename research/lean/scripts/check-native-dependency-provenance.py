"""Verify compilation evidence for the 68 imports with inverted file timestamps.

Regenerating identical source text changes mtime but not the checked content.
Only paired successful compile records are accepted here, never cache inventory.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
OUT = DATA / 'native_dependency_timestamp_verification.json'
LOW = DATA / 'low_cancellation_verification.json'
FRESH = DATA / 'encoder_mtime_gap_peach_verification.json'
LOW_MODULES = ['SpinCodes.Structured.LowCancellationData.Weight1Block0'] + [
    f'SpinCodes.Structured.LowCancellationData.Weight2Block{i}' for i in range(62)]
FRESH_MODULES = {'SpinCodes.Structured.MapSpectrumData.Basis',
                 'SpinCodes.Structured.SparseProgramDefs',
                 'SpinCodes.Structured.SparseMaxima',
                 'SpinCodes.Structured.SparseContributionSound',
                 'SpinCodes.Majorant.RefinedData'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(record):
    tmp = OUT.with_suffix('.tmp')
    tmp.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    tmp.replace(OUT)


def main():
    record = {'status': 'RUNNING', 'scope': 'Current source/object pairs backed by successful '
              'compilation records for timestamp-inverted imports. 63 existing records and '
              '5 fresh isolated source rechecks; no live files modified and no full replay claim.',
              'started_utc': datetime.now(timezone.utc).isoformat(), 'modules': {}}
    save(record)
    try:
        raw = {p: p.read_bytes() for p in (LOW, FRESH)}
        low, fresh = (json.loads(raw[p]) for p in (LOW, FRESH))
        assert low['status'] == fresh['status'] == 'PASS'
        inputs = {p.relative_to(ROOT).as_posix(): hashlib.sha256(b).hexdigest() for p, b in raw.items()}

        def add(module, entry, producer, isolated=False):
            relative = module.replace('.', '/') + '.lean'
            src = ROOT / relative
            obj = (ROOT / '.lake/build/lib/lean' / relative).with_suffix('.olean')
            assert entry['exit_code'] == 0, module
            source_hash, object_hash = sha(src), sha(obj)
            assert source_hash == entry['source_sha256'], module
            assert object_hash == entry['object_sha256'], module
            if isolated:
                assert entry['status'] == 'PASS' and entry['fresh_object_bytes_match_original'], module
                assert sha(ROOT / entry['fresh_object']) == object_hash, module
                assert entry['fresh_object_sha256'] == entry['original_object_sha256'] == object_hash, module
            record['modules'][relative] = {'source_sha256': source_hash, 'object_sha256': object_hash,
                                            'verification_report': producer.relative_to(ROOT).as_posix()}

        for module in LOW_MODULES:
            add(module, low['modules'][module.replace('.', '/') + '.lean'], LOW)
        entries = {entry['module']: entry for entry in fresh['modules']}
        assert set(entries) == FRESH_MODULES and len(fresh['modules']) == 5
        for module, entry in entries.items():
            add(module, entry, FRESH, isolated=True)
        assert len(record['modules']) == 68
        for relative, entry in record['modules'].items():
            assert sha(ROOT / relative) == entry['source_sha256']
            assert sha((ROOT / '.lake/build/lib/lean' / relative).with_suffix('.olean')) == entry['object_sha256']
        assert all(sha(ROOT / path) == digest for path, digest in inputs.items())
        record.update(status='PASS', checked_modules=68, input_report_sha256=inputs,
                      script_sha256=sha(Path(__file__)))
        print('PASS: 68 timestamp inversions resolved by verified current compilation hashes.')
    except BaseException as exc:
        record.update(status='FAIL', error=repr(exc))
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save(record)


if __name__ == '__main__':
    main()
