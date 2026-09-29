"""Replay the concrete map algebra and the finite structural checks.

This is a structural checkpoint, not a claim that all spectrum blocks have
passed or that the distance theorem is closed.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
MODULES = ['PackedMapDefs', 'ConcreteMapData', 'PackedMap', 'ConcreteMaps',
           'ConcreteMapPin', 'MapSpectrumDefs', 'MapSpectrum',
           'MapSpectrumSumDefs', 'MapSpectrumSum']


def main():
    DATA.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['LEAN_PATH'] = os.pathsep.join(str(p) for p in
        [ROOT/'.lake/build/lib/lean'] + sorted((ROOT/'.lake/packages').glob('*/.lake/build/lib/lean')))
    report = dict(status='RUNNING',scope='Concrete map algebra and structural checks; full spectra and distance remain open.', modules=[])
    target = DATA/'structural_verification.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    for name in MODULES:
        source = Path(f'SpinCodes/Structured/{name}.lean')
        digest = hashlib.sha256((ROOT/source).read_bytes()).hexdigest()
        output = (ROOT/'.lake/build/lib/lean'/source).with_suffix('.olean')
        logfile = DATA/f'structural_{name}.log'
        start = time.monotonic()
        with logfile.open('w',encoding='utf-8') as stream:
            result = subprocess.run(['lean','-o',str(output),str(source)],cwd=ROOT,env=env,
                stdout=stream,stderr=subprocess.STDOUT,timeout=900)
        assert digest == hashlib.sha256((ROOT/source).read_bytes()).hexdigest()
        report['modules'].append(dict(path=source.as_posix(),source_sha256=digest,
            exit_code=result.returncode,seconds=round(time.monotonic()-start,2),log=logfile.name))
        if result.returncode:
            report['status'] = 'FAIL'
        target.write_text(json.dumps(report,indent=2)+'\n')
        print(f'{name}: {result.returncode}',flush=True)
        if result.returncode:
            print(logfile.read_text(encoding='utf-8'))
            raise SystemExit(result.returncode)
    report['status'] = 'PASS'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: concrete map structural checkpoint.',flush=True)


if __name__ == '__main__':
    main()
