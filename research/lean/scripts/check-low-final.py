"""Verify fetched low-layer objects locally and run the project invariants."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    report = {'status':'RUNNING', 'checks':[],
              'scope':'Local statement/axiom audit of the fetched objects and default project invariant check. The full distance theorem remains open.'}
    target = DATA/'low_final_verification.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    try:
        semantic = json.loads((DATA/'low_cancellation_verification.json').read_text())
        assert semantic['status'] == 'PASS' and semantic['axiom_audits'] == 16
        assert len(semantic['modules']) == 76
        for rel, entry in semantic['modules'].items():
            assert entry['source_sha256'] == sha(ROOT/rel), rel
            assert entry['object_sha256'] == sha((ROOT/'.lake/build/lib/lean'/rel).with_suffix('.olean')), rel
        original_pin = sha(ROOT/'SpinCodes/Pin.lean')
        assert original_pin == '447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca'
        env = os.environ.copy()
        env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
            [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
        commands = [('low_local_pin.log', ['lean', '-j1', '-M', '10000', 'SpinCodes/Structured/ConcreteLowPin.lean']),
                    ('low_invariant_check.log', [r'C:\Program Files\Git\bin\bash.exe', 'scripts/check.sh'])]
        for log, cmd in commands:
            start = time.monotonic()
            out = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, encoding='utf-8', timeout=1200)
            (DATA/log).write_text(out.stdout+out.stderr, encoding='utf-8')
            report['checks'].append({'log':log, 'exit_code':out.returncode, 'seconds':round(time.monotonic()-start,2)})
            target.write_text(json.dumps(report, indent=2)+'\n')
            assert out.returncode == 0, out.stdout[-6000:]+out.stderr[-1000:]
            if log == 'low_local_pin.log':
                entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", out.stdout)
                assert len(entries) == 16
                for _, closure, _ in entries:
                    assert {s.strip() for s in closure.split(',') if s.strip()} <= {'propext','Classical.choice','Quot.sound'}
            else:
                assert 'ALL CHECKS PASS' in out.stdout
            print('PASS', log, flush=True)
        assert original_pin == sha(ROOT/'SpinCodes/Pin.lean')
        report.update(status='PASS', modules=76, axiom_audits=16, original_pin_sha256=original_pin,
                      semantic_report_sha256=sha(DATA/'low_cancellation_verification.json'))
    except Exception as err:
        report.update(status='FAIL', error=str(err))
        raise
    finally:
        target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
