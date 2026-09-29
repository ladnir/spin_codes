"""Replay the concrete route/encoder bridge, then verify pins locally.

Earlier numerical certificates are hash-checked and reused. All proof replay
and invariant commands run sequentially; this is not a full numerical replay.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
PIN_HASH = '447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca'
NAMES = [
    'FiniteLaw', 'FiniteProductLaw',
    'Structured/ConcreteShuffle', 'Structured/ConcreteRowLaw',
    'Structured/ConcreteShuffleMixture', 'Structured/ConcreteShufflePermutation',
    'Structured/ConcreteShufflePin', 'Structured/ConcreteRouteTranspose',
    'Structured/ConcreteRoute', 'Structured/ConcreteRouteDomination',
    'Structured/ConcreteRoutePermutation', 'Structured/ConcreteRoutePin',
    'Structured/ConcreteEncoder', 'Structured/ConcreteEncoderMoment',
    'Structured/ConcreteEncoderPin', 'Structured/ConcreteReshape',
    'Structured/ConcreteSerialization', 'Structured/ConcreteRoutedMoment',
    'Structured/ConcreteRoutedPin',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(text, expected):
    entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", text)
    assert len(entries) == expected, (len(entries), expected)
    for name, closure, _ in entries:
        assert {s.strip() for s in closure.split(',') if s.strip()} <= {
            'propext', 'Classical.choice', 'Quot.sound'}, (name, closure)
    return [name for name, _, _ in entries]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    target = DATA / 'routed_bridge_verification.json'
    report = {'status': 'RUNNING', 'modules': {}, 'checks': [],
              'scope': 'Nineteen bridge modules replayed on Peach; earlier numerical dependencies reused. Windows consolidated statement/axiom check and default project invariants. Full distance theorem remains open.'}

    def save():
        target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    save()
    try:
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        paper = {str(p.relative_to(ROOT.parent)): sha(p)
                 for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        prior = json.loads((DATA / 'low_cancellation_verification.json').read_text())
        assert prior['status'] == 'PASS' and len(prior['modules']) == 76
        for rel, entry in prior['modules'].items():
            assert sha(ROOT / rel) == entry['source_sha256'], rel
            assert sha((ROOT / '.lake/build/lib/lean' / rel).with_suffix('.olean')) == entry['object_sha256'], rel
        report['prior_modules_hash_checked'] = len(prior['modules'])
        for name in NAMES:
            source = f'SpinCodes/{name}.lean'
            start = time.monotonic()
            result = subprocess.run([sys.executable, 'scripts/peach-lean.py', source],
                                    cwd=ROOT, capture_output=True, encoding='utf-8', timeout=1200)
            assert result.returncode == 0, result.stdout + result.stderr
            report['modules'][source] = {
                'source_sha256': sha(ROOT / source),
                'object_sha256': sha((ROOT / '.lake/build/lib/lean' / source).with_suffix('.olean')),
                'seconds': round(time.monotonic() - start, 2),
                'log': f'peach_{Path(name).name}.log',
            }
            save()
            print('PASS', source, flush=True)
        report['axiom_theorems'] = audit((DATA / 'peach_ConcreteRoutedPin.log').read_text(), 26)
        env = os.environ.copy()
        env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
            [ROOT / '.lake/build/lib/lean'] + sorted((ROOT / '.lake/packages').glob('*/.lake/build/lib/lean')))
        commands = [
            ('routed_local_pin.log', ['lean', '-j1', '-M', '10000', 'SpinCodes/Structured/ConcreteRoutedPin.lean']),
            ('routed_invariant_check.log', [r'C:\Program Files\Git\bin\bash.exe', 'scripts/check.sh']),
        ]
        for log, command in commands:
            print('START', log, flush=True)
            start = time.monotonic()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                    encoding='utf-8', timeout=1200)
            output = result.stdout + result.stderr
            (DATA / log).write_text(output, encoding='utf-8')
            report['checks'].append({'log': log, 'exit_code': result.returncode,
                                     'seconds': round(time.monotonic() - start, 2)})
            save()
            assert result.returncode == 0, output[-8000:]
            if log == 'routed_local_pin.log':
                assert audit(output, 26) == report['axiom_theorems']
            else:
                assert 'ALL CHECKS PASS' in output
            print('PASS', log, flush=True)
        for source, entry in report['modules'].items():
            assert sha(ROOT / source) == entry['source_sha256'], source
            assert sha((ROOT / '.lake/build/lib/lean' / source).with_suffix('.olean')) == entry['object_sha256'], source
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        assert paper == {str(p.relative_to(ROOT.parent)): sha(p)
                         for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        report.update(status='PASS', axiom_audits=26, original_pin_sha256=PIN_HASH,
                      paper_sha256=paper)
    except Exception as err:
        report.update(status='FAIL', error=str(err))
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
