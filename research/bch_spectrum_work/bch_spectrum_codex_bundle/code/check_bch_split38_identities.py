"""Exhaustive small-code check of split transforms and conditional moments."""
import json
import math
from pathlib import Path
from prepare_bch_split38_probe import kraw
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch_split38_identity_toy.json'


def main():
    # RM(1,3), with coordinates indexed by the three-bit vectors.
    q={sum((a^((linear&x).bit_count()%2))<<x for x in range(8))
       for a in range(2) for linear in range(8)}
    assert len(q)==16 and min(v.bit_count() for v in q if v)==4
    dual={v for v in range(256) if all((v&u).bit_count()%2==0 for u in q)}
    assert dual==q
    pivot=3
    partition=lambda v:((v&3).bit_count(),(v>>2).bit_count())
    split={(i,j):sum(partition(v)==(i,j) for v in q) for i in range(3) for j in range(7)}
    ka,kb=kraw(2),kraw(6)
    for u in range(3):
        for v in range(7):
            transform=sum(split[i,j]*ka[i][u]*kb[j][v] for i in range(3) for j in range(7))
            actual=sum(partition(word)==(u,v) for word in dual)
            assert transform==len(q)*actual
    expected_coset=[sum((word^pivot).bit_count()==w for word in q) for w in range(9)]
    for v in range(256):
        if v.bit_count()%2==0 and v not in q:
            assert [sum((word^v).bit_count()==w for word in q) for w in range(9)]==expected_coset
    for w in range(9):
        assert sum(split[i,j] for i in range(3) for j in range(7) if 2-i+j==w)==expected_coset[w]
    for pattern in range(4):
        fiber=[v>>2 for v in q if v&3==pattern]
        assert len(fiber)==4
        for j in range(6):
            assert sum((v>>j)&1 for v in fiber)==2
    for i in range(3):
        assert sum(split[i,j] for j in range(7))==math.comb(2,i)*4
        assert sum(split[i,j]*kb[j][1] for j in range(7))==0
    result=dict(classification='Exhaustive exact small-code split identity check',
                code='RM(1,3) [8,4,4]',codewords=16,split_transform_entries_checked=21,
                all_nonzero_even_cosets_have_same_spectrum=True,all_four_projection_fibers_checked=True,
                local_moments_checked=True,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                                                        (Path(__file__),ROOT/'code/prepare_bch_split38_probe.py')})
    if OUTPUT.exists():
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
