"""Check the final concrete-count assembly using already replayed dependencies."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    report = dict(status='RUNNING', scope='Final concrete spectrum, kernel counts, fiber caps, pair counts and shell assembly. Imported numerical certificates are reused; the encoder one-step law remains open.', checks=[])
    target = DATA / 'concrete_counts_verification.json'
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT / '.lake/build/lib/lean'] + sorted((ROOT / '.lake/packages').glob('*/.lake/build/lib/lean')))
    try:
        for module in ['FiberNumericsAll', 'ConcreteCounts', None]:
            if module == 'ConcreteCounts':
                import concrete_counts_assemble
                concrete_counts_assemble.main()
            source = (f'SpinCodes/Structured/{module}.lean' if module else 'scripts/concrete_counts_axioms.lean')
            output = (ROOT / '.lake/build/lib/lean' / source).with_suffix('.olean') if module else None
            command = ['lean', '-j1', '-M', '10000'] + (['-o', str(output)] if output else []) + [source]
            start = time.monotonic()
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, encoding='utf-8', timeout=900)
            log = (module or 'concrete_counts_axioms') + '.log'
            (DATA / log).write_text(result.stdout, encoding='utf-8')
            entry = dict(source=source, source_sha256=sha(ROOT/source), exit_code=result.returncode,
                         seconds=round(time.monotonic()-start, 2), log=log)
            if output and result.returncode == 0:
                entry['output_sha256'] = sha(output)
            report['checks'].append(entry)
            target.write_text(json.dumps(report, indent=2)+'\n')
            assert result.returncode == 0, result.stdout[-5000:]
            print(f'PASS {source} ({entry["seconds"]} s)', flush=True)
        entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", result.stdout)
        expected = re.findall(r'#print axioms (\S+)', (ROOT/'scripts/concrete_counts_axioms.lean').read_text())
        assert len(entries) == len(expected) and {e[0] for e in entries} == set(expected)
        for _, closure, _ in entries:
            assert {s.strip() for s in closure.split(',') if s.strip()} <= {'propext', 'Classical.choice', 'Quot.sound'}, closure
        report.update(status='PASS', axiom_audits=len(entries))
        print(f'PASS: concrete-count connection and {len(entries)} standard-axiom audits.')
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
