"""Sequentially replay the low-cancellation modules on Peach, recording hashes.

Use --blocks for the finite witnesses or --assembly for the semantic assembly.
Already successful blocks with unchanged source/object hashes may be resumed.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT/'scripts/map_data/low_cancellation_verification.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--blocks', action='store_true')
    ap.add_argument('--assembly', action='store_true')
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    report = json.loads(REPORT.read_text()) if args.resume and REPORT.exists() else {'modules': {}}
    paths = []
    if args.blocks:
        for j in (1,2):
            paths += sorted((ROOT/'SpinCodes/Structured/LowCancellationData').glob(f'Weight{j}Block*.lean'),
                            key=lambda p: int(p.stem.split('Block')[1]))
    if args.assembly:
        paths += [ROOT/f'SpinCodes/Structured/{s}.lean' for s in
                  ['LowCancellationDefs', 'LowCancellationSort', 'LowCancellationBridge',
                   'LowCancellationMoments', 'LowCancellationShells', 'LowCancellationRows',
                   'LowCancellationData/Weight1', 'LowCancellationData/Weight2',
                   'ConcreteLowCancellation', 'ConcreteTransferAll', 'ConcreteBernoulli',
                   'ConcreteOccupation', 'ConcreteLowPin']]
    report['status'] = 'RUNNING'
    for p in paths:
        rel = p.relative_to(ROOT).as_posix()
        obj = ROOT/'.lake/build/lib/lean'/p.relative_to(ROOT).with_suffix('.olean')
        previous = report['modules'].get(rel, {})
        if args.resume and previous.get('exit_code') == 0 and previous.get('source_sha256') == sha(p) and previous.get('object_sha256') == sha(obj):
            print('RESUME', rel, flush=True)
            continue
        start = time.monotonic()
        out = subprocess.run([sys.executable, str(ROOT/'scripts/peach-lean.py'), rel],
                             cwd=ROOT, capture_output=True, encoding='utf-8')
        entry = {'exit_code':out.returncode, 'seconds':round(time.monotonic()-start, 2),
                 'source_sha256':sha(p), 'object_sha256':sha(obj) if out.returncode == 0 else None}
        report['modules'][rel] = entry
        report['status'] = 'RUNNING' if out.returncode == 0 else 'FAIL'
        REPORT.write_text(json.dumps(report, indent=2)+'\n')
        print(('PASS' if out.returncode == 0 else 'FAIL'), rel, entry['seconds'], flush=True)
        if out.returncode:
            print(out.stdout, out.stderr, flush=True)
            raise SystemExit(out.returncode)
    if args.assembly:
        log = (ROOT/'scripts/map_data/peach_ConcreteLowPin.log').read_text(encoding='utf-8')
        entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", log)
        expected = re.findall(r'#print axioms (\S+)', (ROOT/'SpinCodes/Structured/ConcreteLowPin.lean').read_text(encoding='utf-8'))
        assert len(entries) == len(expected) == 16 and {e[0] for e in entries} == set(expected)
        for _, closure, _ in entries:
            assert {s.strip() for s in closure.split(',') if s.strip()} <= {'propext', 'Classical.choice', 'Quot.sound'}, closure
        report['axiom_audits'] = len(entries)
    report['status'] = 'PASS'
    report['scope'] = 'Only the modules listed here; inherited project dependencies are reused.'
    REPORT.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
