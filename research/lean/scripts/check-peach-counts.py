"""Linux check of the final count assembly against hash-pinned local dependencies."""
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
    report = dict(status='RUNNING', scope='Linux compilation of the final count assembly, reusing hash-checked local project dependencies and Linux mathlib cache.', checks=[])
    target = DATA / 'peach_counts_verification.json'
    try:
        manifest = json.loads((DATA / 'peach_dependencies.json').read_text())
        for record in manifest['files']:
            assert sha(ROOT / record['path']) == record['sha256'], record['path']
        report['dependency_modules'] = manifest['modules']
        env = os.environ.copy()
        env['LEAN_PATH'] = ':'.join(str(p) for p in [ROOT / '.lake/build/lib/lean'] +
                                  sorted((ROOT / '.lake/packages').glob('*/.lake/build/lib/lean')))
        for module in ['FiberNumericsAll', 'ConcreteCounts', None]:
            source = f'SpinCodes/Structured/{module}.lean' if module else 'scripts/concrete_counts_axioms.lean'
            output = (ROOT / '.lake/build/lib/lean' / source).with_suffix('.olean') if module else None
            tag = module or 'concrete_counts_axioms'
            command = ['lean', '-j1', '-M', '20000']
            command += (['-o', str(output)] if output else []) + [source]
            start = time.monotonic()
            log = f'peach_{tag}.log'
            with (DATA / log).open('w') as stream:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
                while True:
                    pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                    if pid:
                        process.returncode = os.waitstatus_to_exitcode(status)
                        break
                    if time.monotonic() - start > 900:
                        process.kill()
                        process.wait()
                        raise TimeoutError(source)
                    time.sleep(0.25)
            text = (DATA / log).read_text()
            entry = dict(source=source, source_sha256=sha(ROOT/source), exit_code=process.returncode,
                         seconds=round(time.monotonic()-start, 2), peak_rss_kib=usage.ru_maxrss, log=log)
            if output and process.returncode == 0:
                entry['output_sha256'] = sha(output)
            report['checks'].append(entry)
            target.write_text(json.dumps(report, indent=2)+'\n')
            assert process.returncode == 0, text[-5000:]
            print(f'PASS {tag}: {entry["seconds"]} s, peak {entry["peak_rss_kib"] / 1024**2:.2f} GiB', flush=True)
        entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))", text)
        expected = re.findall(r'#print axioms (\S+)', (ROOT/'scripts/concrete_counts_axioms.lean').read_text())
        assert len(entries) == len(expected) and {e[0] for e in entries} == set(expected)
        for _, closure, _ in entries:
            assert {s.strip() for s in closure.split(',') if s.strip()} <= {'propext', 'Classical.choice', 'Quot.sound'}, closure
        report.update(status='PASS', axiom_audits=len(entries))
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
