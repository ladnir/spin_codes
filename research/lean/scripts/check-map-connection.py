"""Compile connection statement pins and audit the imported proof closures.

The semantic modules are built separately with lean-local.ps1. This audit
checks their current source/object timestamps and records hashes; it does
not rerun the numerical worker dependencies or claim a full spectrum proof.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
MODULES = ['FiberNumericsDefs','FiberNumericsData/Tables','FiberNumerics','FiberNumericsBridge',
           'FiberFrozenData','FiberFrozen','FiberNumericsSample','ConcreteShells',
           'NonnegativeInduction','LiveInduction','SparseMomentBridge','ConcreteConnectionPin']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
    report = dict(status='RUNNING',scope='Statement-pin compilation and axiom audit for the numerical semantic connection, shell construction, and live nonnegative induction. The concrete spectrum and one-step law remain open.')
    target = DATA/'connection_audit.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    pin = 'SpinCodes/Structured/ConcreteConnectionPin.lean'
    for source,output,log in [(pin,str((ROOT/'.lake/build/lib/lean'/pin).with_suffix('.olean')),'connection_pin.log'),
                              ('scripts/concrete_connection_axioms.lean',None,'connection_axioms.log')]:
        command = ['lean'] + (['-o',output] if output else []) + [source]
        result = subprocess.run(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,encoding='utf-8',timeout=900)
        (DATA/log).write_text(result.stdout,encoding='utf-8')
        if result.returncode:
            report['status'] = 'FAIL'
            target.write_text(json.dumps(report,indent=2)+'\n')
            print(result.stdout)
            raise SystemExit(result.returncode)
    entries = re.findall(r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|(does not depend on any axioms))",result.stdout)
    expected = re.findall(r'#print axioms (\S+)',(ROOT/'scripts/concrete_connection_axioms.lean').read_text())
    assert len(entries) == len(expected) == 16 and {e[0] for e in entries} == set(expected), result.stdout
    for _,closure,_ in entries:
        assert {part.strip() for part in closure.split(',') if part.strip()} <= {'propext','Classical.choice','Quot.sound'}, closure
    snapshots = []
    for name in MODULES:
        source = ROOT/f'SpinCodes/Structured/{name}.lean'
        output = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
        assert output.stat().st_mtime_ns >= source.stat().st_mtime_ns, f'Rebuild {name} first.'
        snapshots.append(dict(path=source.relative_to(ROOT).as_posix(),source_sha256=sha(source),output_sha256=sha(output)))
    report.update(status='PASS',modules_snapshot=snapshots,axiom_audits=len(entries),
        logs=['connection_pin.log','connection_axioms.log'])
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: connection pins and 16 standard-axiom audits. Spectrum/one-step premises remain explicit.')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        target = DATA/'connection_audit.json'
        report = json.loads(target.read_text()) if target.exists() else {}
        report.update(status='FAIL',error=str(error))
        target.write_text(json.dumps(report,indent=2)+'\n')
        raise
