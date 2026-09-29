"""Record successful individual Lean checks without claiming a fresh dependency replay."""
from pathlib import Path
import hashlib,json,re
R=Path(__file__).resolve().parents[1]
D=R/'scripts/map_data'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def record(name,stems,pin,n,scope):
    rows={}
    for stem in stems:
        p=R/f'SpinCodes/Structured/{stem}.lean'
        o=R/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
        l=D/f'peach_{stem}.log'
        assert o.exists() and o.stat().st_mtime>=p.stat().st_mtime,stem
        assert l.exists() and ': error:' not in l.read_text(encoding='utf-8'),stem
        rows[p.relative_to(R).as_posix()]={'source_sha256':sha(p),'object_sha256':sha(o),'log':l.name}
    audits=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",(D/f'peach_{pin}.log').read_text(encoding='utf-8'))
    assert len(audits)==n,(pin,len(audits))
    for t,a in audits:
        assert set(x.strip() for x in a.split(','))<={'propext','Classical.choice','Quot.sound'},t
    (D/name).write_text(json.dumps({'status':'PASS','scope':scope+' Modules individually checked on Peach; existing checked dependencies reused, not a fresh closure replay.','modules':rows,'axiom_audits':n,'axiom_theorems':[t for t,a in audits]},indent=2)+'\n')
    print(name,len(rows),'modules,',n,'standard-only audits')
record('outer_tail_closure_verification.json', sorted(p.stem for p in (R/'SpinCodes/Structured').glob('ConcreteOuterTail*.lean')),'ConcreteOuterTailClosurePin',10,'Unconditional actual Golay-BAA outerTail and nativeTail limits, eventual positive good event, selection failure limit, and growing sparse regime.')
record('outer_linearity_verification.json',['ConcreteOuterLinearity','ConcreteOuterGolayLinearity','ConcreteOuterNativeLinearity','ConcreteOuterLinearityPin'],'ConcreteOuterLinearityPin',7,'Actual Golay generator fold, accumulator, shared outer and native-row binary linearity.')
record('native_linear_code_verification.json',['ConcreteNativeLinearCode','ConcreteNativeLinearCodeDistance','ConcreteNativeLinearCodePin'],'ConcreteNativeLinearCodePin',9,'Actual native linear map, range submodule, exact dimension/rate, minimum/pairwise Hamming distance and bad-event equivalence.')
record('marked_uniform_subset_verification.json',['ConcreteMarkedUniformSubset','ConcreteMarkedUniformSubsetPin'],'ConcreteMarkedUniformSubsetPin',2,'Exact fair marked-bit law as pushforward of uniform subsets along any embedding, including zero marks.')
