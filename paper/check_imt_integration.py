"""Check IMT manuscript transcription against retained 11% evidence.

Not an interval replay or an independent analytic proof review.
"""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
THEORY=ROOT/'workstreams/inner_design/imt_asymptotic'
sys.path.insert(0,str(THEORY))
import screen


def main():
    if not __debug__:raise RuntimeError('Do not disable checks with -O.')
    model=screen.Model().engine
    paper=ROOT/'paper'
    text=(paper/'structured_imt_appendix.tex').read_text()
    words=re.findall(r'^(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}',text,re.M)
    assert [int(i) for i,_,_ in words]==list(range(19))
    for j,(_,a,c) in enumerate(words):
        assert len(a)==len(c)==32
        assert int(a,16)==sum(((col>>j)&1)<<i for i,col in enumerate(model.a_columns))
        assert int(c,16)==sum(((col>>j)&1)<<i for i,col in enumerate(model.columns))
    assert model.spectrum=={48:5166,56:110288,64:293455,72:110128,80:5250}
    for count in model.spectrum.values():assert str(count) in text
    records=['d11/ASSEMBLY_FINAL.json','d11/OUTER_REFINED.json','d11/OUTER_REPLAY.json',
             'd11/DENSE_OUTWARD.json','d11/DENSE_REPLAY.json','ASSEMBLY_D1099.json']
    for name in records:
        record=json.loads((THEORY/name).read_text())
        for path,digest in record['source_sha256'].items():
            target=(ROOT/path).resolve()
            assert target.is_relative_to(ROOT)
            assert hashlib.sha256(target.read_bytes()).hexdigest()==digest,path
    final=json.loads((THEORY/records[0]).read_text())
    assert final['delta']=='11/100' and final['block_constant']=='39/4'
    assert final['inner']==model.identity()['inner']
    assert final['dense_boxes']==1023 and final['outer_segments']==39
    assert final['sparse_cutoff']==4096
    for value in (r'1{,}023',r'5{,}363','4096','-0.006','-0.0008679'):
        assert value in text or value in (paper/'structured_appendix.tex').read_text(),value
    body=(paper/'structured_spin.tex').read_text()
    assert r'Q_{i+1}:=M_iQ_i+C(X_i)' in body
    assert r'r_i=A^\top U_i+M_i^\top r_{i+1}' in body
    assert 'RM2Sub-S19' not in body
    assert r'\sum_{a=1}^b\sum_{c=0}^b' in (paper/'structured_appendix.tex').read_text()
    assert r'Q_{i+1}=M_iQ_i+C(X_i)' in (paper/'finite_certificates.tex').read_text()
    assert '10.110' in (paper/'implementation.tex').read_text()
    print('IMT manuscript integration PASS: 38 map words, 11% evidence bindings, and shared finite inner.')


if __name__=='__main__':main()
