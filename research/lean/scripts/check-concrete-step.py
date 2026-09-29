"""Sequential Peach replay of the concrete one-step semantic modules and pins.

Numerical and existing mathematical dependencies are reused, hash-pinned by
peach-lean.py. The final theorem excludes input weights 1 and 2 explicitly.
"""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
MODULES = ['ConcreteMoments', 'ConcreteCancellation', 'ConcreteStep', 'ConcreteRows',
           'ConcreteCancellationRows', 'ConcreteZeroRow', 'ConcreteShellRows',
           'LiveKernel', 'ConcreteClosedShells', 'ConcreteTransfer', 'ConcreteStepPin']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    report = dict(status='RUNNING', scope='Fresh sequential semantic replay and statement/axiom audit for the concrete fixed-weight one-step law, excluding weights 1 and 2. Numerical dependencies are reused. The full distance theorem remains open.', checks=[])
    target = DATA/'concrete_step_verification.json'
    try:
        for name in MODULES:
            source = ROOT/f'SpinCodes/Structured/{name}.lean'
            start = time.monotonic()
            result = subprocess.run([sys.executable, str(ROOT/'scripts/peach-lean.py'), str(source)],
                cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8', timeout=1200)
            obj = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
            entry = dict(source=source.relative_to(ROOT).as_posix(), source_sha256=sha(source),
                         exit_code=result.returncode, seconds=round(time.monotonic()-start, 2), log=f'peach_{name}.log')
            if result.returncode == 0:
                entry['output_sha256'] = sha(obj)
            report['checks'].append(entry)
            target.write_text(json.dumps(report, indent=2)+'\n')
            assert result.returncode == 0, result.stdout[-8000:]
            print(f'PASS {name} ({entry["seconds"]} s including transfer)', flush=True)
        log = (DATA/'peach_ConcreteStepPin.log').read_text(encoding='utf-8')
        entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", log)
        expected = re.findall(r'#print axioms (\S+)', (ROOT/'SpinCodes/Structured/ConcreteStepPin.lean').read_text(encoding='utf-8'))
        assert len(entries) == len(expected) == 17 and {e[0] for e in entries} == set(expected)
        for _, closure, _ in entries:
            assert {s.strip() for s in closure.split(',') if s.strip()} <= {'propext', 'Classical.choice', 'Quot.sound'}, closure
        report.update(status='PASS', axiom_audits=len(entries))
        print('PASS: concrete step for 127 weights and 17 standard-axiom audits.')
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
