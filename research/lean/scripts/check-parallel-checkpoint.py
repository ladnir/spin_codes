"""Replay the stable parallel finite-proof checkpoint, not the full certificate tree.

All 44 modules added since the routed checkpoint are replayed in dependency order.
Earlier checked numerical dependencies are reused after hash verification. Workers
may continue in new modules that are outside this fixed checkpoint's import graph.
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def closure(module, ordered):
    if module in ordered:
        return
    src = ROOT / (module.replace('.', '/') + '.lean')
    for dep in re.findall(r'^import (SpinCodes\S*)', src.read_text(encoding='utf-8'), re.M):
        closure(dep, ordered)
    ordered[module] = src


def audit(output, expected):
    entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", output)
    assert len(entries) == expected, (len(entries), expected)
    for name, axioms, _ in entries:
        assert {s.strip() for s in axioms.split(',') if s.strip()} <= {
            'propext', 'Classical.choice', 'Quot.sound'}, (name, axioms)
    return [name for name, _, _ in entries]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    target = DATA / 'parallel_finite_verification.json'
    report = {'status': 'RUNNING', 'modules': {}, 'checks': [],
              'scope': 'Stable finite concrete parallel checkpoint replayed on Peach; earlier numerical dependencies reused. Windows combined pin and default invariant check. Unconditional distance theorem remains open.'}

    def save():
        target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    save()
    try:
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        paper = {str(p.relative_to(ROOT.parent)): sha(p)
                 for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        prior_count = 0
        for prior_name in ['low_cancellation_verification.json', 'routed_bridge_verification.json']:
            prior = json.loads((DATA / prior_name).read_text())
            assert prior['status'] == 'PASS'
            for rel, entry in prior['modules'].items():
                assert sha(ROOT / rel) == entry['source_sha256'], rel
                assert sha((ROOT / '.lake/build/lib/lean' / rel).with_suffix('.olean')) == entry['object_sha256'], rel
                prior_count += 1
        report['prior_module_records_hash_checked'] = prior_count
        old, new = {}, {}
        closure('SpinCodes', old)
        closure('SpinCodes.Structured.ConcreteRoutedPin', old)
        closure('SpinCodes.Structured.ParallelClosurePin', new)
        selected = {m: s for m, s in new.items() if m not in old}
        assert len(selected) == 44, len(selected)
        snapshot = {m: sha(s) for m, s in selected.items()}
        report['axiom_theorems'] = []
        for module, src in selected.items():
            rel = src.relative_to(ROOT).as_posix()
            start = time.monotonic()
            run = subprocess.run([sys.executable, 'scripts/peach-lean.py', rel], cwd=ROOT,
                                 capture_output=True, encoding='utf-8', timeout=1200)
            assert run.returncode == 0, run.stdout + run.stderr
            assert sha(src) == snapshot[module], rel
            report['modules'][rel] = {
                'source_sha256': sha(src),
                'object_sha256': sha((ROOT / '.lake/build/lib/lean' / rel).with_suffix('.olean')),
                'seconds': round(time.monotonic() - start, 2),
                'log': f'peach_{src.stem}.log',
            }
            count = len(re.findall(r'^#print axioms ', src.read_text(encoding='utf-8'), re.M))
            if count:
                report['axiom_theorems'].extend(audit(run.stdout, count))
            save()
            print('PASS', rel, flush=True)
        env = os.environ.copy()
        env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
            [ROOT / '.lake/build/lib/lean'] + sorted((ROOT / '.lake/packages').glob('*/.lake/build/lib/lean')))
        for log, command in [
            ('parallel_local_pin.log', ['lean', '-j1', '-M', '10000', 'SpinCodes/Structured/ParallelClosurePin.lean']),
            ('parallel_invariant_check.log', [r'C:\Program Files\Git\bin\bash.exe', 'scripts/check.sh']),
        ]:
            print('START', log, flush=True)
            start = time.monotonic()
            run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                 encoding='utf-8', timeout=1200)
            output = run.stdout + run.stderr
            (DATA / log).write_text(output, encoding='utf-8')
            report['checks'].append({'log': log, 'exit_code': run.returncode,
                                     'seconds': round(time.monotonic() - start, 2)})
            save()
            assert run.returncode == 0, output[-8000:]
            if log == 'parallel_local_pin.log':
                audit(output, 4)
            else:
                assert 'ALL CHECKS PASS' in output
            print('PASS', log, flush=True)
        for module, src in selected.items():
            assert sha(src) == snapshot[module], module
        assert sha(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH
        assert paper == {str(p.relative_to(ROOT.parent)): sha(p)
                         for p in sorted((ROOT.parent / 'paper').glob('*.tex'))}
        report.update(status='PASS', axiom_audits=len(report['axiom_theorems']),
                      original_pin_sha256=PIN_HASH, paper_sha256=paper)
    except Exception as err:
        report.update(status='FAIL', error=str(err))
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
