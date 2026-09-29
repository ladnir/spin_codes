"""Check the concrete Fourier/fiber connection without rebuilding spectrum dependencies."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
MODULES = ['ConcretePairing', 'ConcreteFourier', 'ConcreteVariance', 'ConcreteKernel',
           'ConcreteWeightSums', 'FiberCertificate', 'ConcreteFiberBounds', 'ConcreteFiberPin']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
    report = dict(status='RUNNING', scope='Concrete character pairing, Fourier/Parseval identities, '
        'kernel distance and four fiber bounds. Full spectra and the distance theorem remain open.', modules=[])
    target = DATA/'fiber_verification.json'
    def save():
        target.write_text(json.dumps(report,indent=2)+'\n')
    save()
    for name in MODULES:
        source = Path(f'SpinCodes/Structured/{name}.lean')
        digest = sha(ROOT/source)
        output = (ROOT/'.lake/build/lib/lean'/source).with_suffix('.olean')
        logfile = DATA/f'fiber_{name}.log'
        start = time.monotonic()
        with logfile.open('w',encoding='utf-8') as stream:
            result = subprocess.run(['lean','-o',str(output),str(source)],cwd=ROOT,env=env,
                stdout=stream,stderr=subprocess.STDOUT,timeout=900)
        assert digest == sha(ROOT/source)
        report['modules'].append(dict(path=source.as_posix(),source_sha256=digest,
            output_sha256=sha(output) if result.returncode == 0 else None,
            exit_code=result.returncode,seconds=round(time.monotonic()-start,2),log=logfile.name))
        if result.returncode:
            report['status'] = 'FAIL'
        save()
        print(f'{name}: {result.returncode}',flush=True)
        if result.returncode:
            print(logfile.read_text(encoding='utf-8'))
            raise SystemExit(result.returncode)
    source = ROOT/'scripts/concrete_fibers_axioms.lean'
    audit = subprocess.run(['lean',str(source)],cwd=ROOT,env=env,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding='utf-8',timeout=900)
    (DATA/'fiber_axioms.log').write_text(audit.stdout,encoding='utf-8')
    closures = re.findall(r'depends on axioms: \[([^\]]*)\]',audit.stdout)
    assert audit.returncode == 0 and len(closures) == 16, audit.stdout
    for closure in closures:
        assert set(filter(None,closure.split(', '))) <= {'propext','Classical.choice','Quot.sound'}, closure
    forbidden = re.compile(r'\b(?:sorry|admit|native_decide)\b|^\s*axiom\s',re.M)
    for module in (ROOT/'SpinCodes').rglob('*.lean'):
        assert not forbidden.search(module.read_text(encoding='utf-8')), module
    report['audit'] = dict(path='scripts/concrete_fibers_axioms.lean',source_sha256=sha(source),
        log='fiber_axioms.log',log_sha256=sha(DATA/'fiber_axioms.log'),theorems=len(closures),
        allowed_axioms=['propext','Classical.choice','Quot.sound'])
    report['status'] = 'PASS'
    save()
    print('PASS: concrete Fourier/fiber connection and standard-axiom audit.',flush=True)


if __name__ == '__main__':
    main()
