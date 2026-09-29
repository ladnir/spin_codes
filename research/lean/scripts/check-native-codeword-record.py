"""Record individually compiled native assembly and actual codeword modules."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
STEMS = ['EventualRegimes','ConcreteNativeClosure','ConcreteNativeClosurePin',
         'ConcreteEncoderOutput','ConcreteRoutedOutput','ConcreteNativeCodeword',
         'ConcreteNativeCodewordPin','ConcreteNativeTwoRegimes']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records, audits = {}, []
    for stem in STEMS:
        src = ROOT/f'SpinCodes/Structured/{stem}.lean'
        obj = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
        log = DATA/f'peach_{stem}.log'
        assert obj.exists() and obj.stat().st_mtime >= src.stat().st_mtime, stem
        output = log.read_text(encoding='utf-8')
        assert ': error:' not in output, stem
        assert not re.search(r'\b(sorry|admit|axiom|native_decide|trustCompiler)\b',src.read_text(encoding='utf-8')), stem
        records[src.relative_to(ROOT).as_posix()] = {
            'source_sha256':sha(src),'object_sha256':sha(obj),'log':log.name}
        audits.extend(re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",output))
    assert len(audits)==22, len(audits)
    for theorem, axioms in audits:
        assert set(a.strip() for a in axioms.split(',')) <= {'propext','Classical.choice','Quot.sound'}, theorem
    (DATA/'native_codeword_verification.json').write_text(json.dumps({
        'status':'PASS','scope':'Eight modules individually compiled on Peach; checked dependencies reused. '
        'Actual emitted codeword/injectivity/weight/cardinality and eventual two-regime assembly. '
        'Fixed and positive occupation remain hypotheses. Separate outer worker modules provide full linearity.',
        'modules':records,'axiom_audits':len(audits),'axiom_theorems':[t for t,_ in audits],
    },indent=2)+'\n',encoding='utf-8')
    print('PASS: eight modules, twenty-two standard-only audits; fixed/dense hypotheses remain.')


if __name__=='__main__':
    main()
