"""Run the final default invariants only after the actual theorem/corollaries pass.

This complements their dedicated checks; the default build alone does not
import the actual native theorem. Existing dependency objects are reused.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
REPORT = DATA / 'native_final_invariants_verification.json'
PIN_HASH = '447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca'
COMMAND = [r'C:\Program Files\Git\bin\bash.exe', 'scripts/check.sh']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(record):
    tmp = REPORT.with_suffix('.tmp')
    tmp.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    for attempt in range(40):
        try:
            tmp.replace(REPORT)
            return
        except PermissionError:
            if attempt == 39:
                raise
            time.sleep(.05)


def verify_pairs(base, consequences):
    records = dict(base['modules'])
    for entry in consequences['modules']:
        records[f"SpinCodes/Structured/{entry['module']}.lean"] = entry
    for relative, entry in records.items():
        obj = (ROOT / '.lake/build/lib/lean' / relative).with_suffix('.olean')
        assert sha(ROOT / relative) == entry['source_sha256'], relative
        assert sha(obj) == entry.get('object_sha256', entry.get('olean_sha256')), relative
    return len(records)


def main():
    record = {'status': 'RUNNING', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Final default build, original pin, prohibited proof constructs, and default '
                       'axiom checks after dedicated actual native theorem/corollary checks. '
                       'Cached dependencies reused; not a fresh full source replay.',
              'command': COMMAND}
    save(record)
    try:
        names = ['native_theorem_verification.json',
                 'encoder_native_distance_consequences_verification.json']
        reports = [json.loads((DATA / name).read_bytes()) for name in names]
        assert all(r['status'] == 'PASS' for r in reports), 'Actual final theorem/corollaries not PASS'
        inputs = {name: sha(DATA / name) for name in names}
        assert reports[1]['base_report_sha256'] == inputs[names[0]], 'Corollaries base report differs'
        count = verify_pairs(*reports)
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        papers = {str(p.relative_to(ROOT.parent)): sha(p)
                  for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        baseline = json.loads((DATA / 'parallel_finite_verification.json').read_bytes())
        assert papers == baseline['paper_sha256'], 'Paper changed'
        checker_hash = sha(ROOT / 'scripts/check.sh')
        result = subprocess.run(COMMAND, cwd=ROOT, capture_output=True,
                                encoding='utf-8', errors='replace', timeout=1200)
        output = result.stdout + result.stderr
        log = DATA / 'native_final_invariants.log'
        log.write_text(output, encoding='utf-8')
        record.update(exit_code=result.returncode, log=log.name, log_sha256=sha(log))
        assert result.returncode == 0 and 'ALL CHECKS PASS' in output, output[-8000:]
        assert sha(ROOT / 'scripts/check.sh') == checker_hash, 'Invariant checker changed'
        assert verify_pairs(*reports) == count
        assert all(sha(DATA / name) == value for name, value in inputs.items()), 'Final report changed'
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        assert papers == {str(p.relative_to(ROOT.parent)): sha(p)
                          for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        record.update(status='PASS', original_pin_sha256=PIN_HASH, paper_sha256=papers,
                      checker_sha256=checker_hash, current_pairs_verified=count,
                      base_report_sha256=inputs[names[0]], corollaries_report_sha256=inputs[names[1]])
        print('PASS: final default invariants and actual theorem/corollary artifact preservation.')
    except BaseException as exc:
        record.update(status='FAIL', error=repr(exc))
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save(record)


if __name__ == '__main__':
    main()
