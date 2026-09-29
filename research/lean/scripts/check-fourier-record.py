"""Record checked Fourier modules and audit the dedicated pin's proof dependencies.

Run peach-lean.py on modules in this order before recording. This records those
individual kernel checks; it does not claim a fresh dependency replay or any
completed numerical Fourier boxes.
"""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
STEMS = [
    'ConcreteFourierProduct', 'ConcreteFourierProbability', 'ConcreteProductTilt',
    'ConcreteFourierTilt', 'ConcreteFourierOverlap', 'ConcreteFourierSpectrum',
    'ConcreteFourierRow', 'ConcreteFourierTransfer', 'ConcreteFourierRouted',
    'ConcreteFourierPin',
]
DATA = ROOT/'scripts/map_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records = {}
    for stem in STEMS:
        src = ROOT/f'SpinCodes/Structured/{stem}.lean'
        obj = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
        log = DATA/f'peach_{stem}.log'
        assert obj.exists() and obj.stat().st_mtime >= src.stat().st_mtime, stem
        assert log.exists() and ': error:' not in log.read_text(encoding='utf-8'), stem
        assert not re.search(r'\b(sorry|admit|axiom|native_decide|trustCompiler)\b',
                             src.read_text(encoding='utf-8')), stem
        records[src.relative_to(ROOT).as_posix()] = {
            'source_sha256': sha(src), 'object_sha256': sha(obj), 'log': log.name,
        }
    audits = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",
                        (DATA/'peach_ConcreteFourierPin.log').read_text(encoding='utf-8'))
    assert len(audits) == 15, len(audits)
    for theorem, axioms in audits:
        assert set(a.strip() for a in axioms.split(',')) <= {
            'propext', 'Classical.choice', 'Quot.sound'}, (theorem, axioms)
    result = {
        'status': 'PASS',
        'scope': 'Ten modules individually compiled on Peach; checked dependencies reused. '
                 'Actual Fourier transfer and routed probability bridge. Zero numerical Fourier boxes checked.',
        'modules': records, 'axiom_audits': len(audits),
        'axiom_theorems': [theorem for theorem, _ in audits],
        'numerical_fourier_boxes_checked': 0,
    }
    (DATA/'fourier_transfer_verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: ten Fourier modules, fifteen standard-only axiom audits; zero numerical Fourier boxes.')


if __name__ == '__main__':
    main()
