"""Exact complementary-support bounds for shortened BCH dimensions.

Replays the dual-distance premise and checks Q <= C <= P against the
production generator. No files or production defaults are changed.
"""
import json
from pathlib import Path
import re
import sys
from joint_support import span, gf2_rank
from shortened_bound import dimension_caps


def complement_caps(primal, dual, k):
    n = len(primal)-1
    assert len(dual) == n+1
    result = [min(primal[v], v-(n-k)+dual[n-v]) for v in range(n+1)]
    assert all(0 <= x <= min(v,k) for v,x in enumerate(result))
    # Maximum shortened dimension is nondecreasing and grows by at most one.
    for v in range(n-1,-1,-1):
        result[v] = min(result[v],result[v+1])
    for v in range(1,n+1):
        result[v] = min(result[v],result[v-1]+1)
    return result


def verify_bch_premise():
    repo = Path(__file__).resolve().parents[3]
    bundle = repo/'research/bch_spectrum_work/bch_spectrum_codex_bundle'
    sys.path.insert(0,str(bundle/'code'))
    from certify_bch_shift_rank_split import build
    from certify_bch_shift_rank import in_span
    from bch_quotient import generator_polynomial, binary_poly_divmod
    proof = build()
    assert proof == json.loads((bundle/'generated/bch256_shift_rank_q30_refined.json').read_text())
    assert proof['Qdual_extended_distance_lower'] >= 30
    header = repo/'spin/src/kernels/generated/BchCircuit.h'
    words = [int(x,16) for x in re.findall(r'0x([0-9a-f]+)ULL',header.read_text())]
    assert len(words) == 512
    rows = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    assert gf2_rank(rows) == 128
    p,q = generator_polynomial(37),generator_polynomial(39)
    extend = lambda w: w | ((w.bit_count()%2) << 255)
    assert all(in_span(extend(q << i),rows) for i in range(123))
    assert all(w.bit_count()%2 == 0 and binary_poly_divmod(w&((1<<255)-1),p)[1] == 0 for w in rows)
    assert in_span((1<<256)-1,rows)  # C dual is even; dual LP bounds apply.
    print('Replayed dual distance >=30; production rank128, Q<=C<=P and all-ones containment verified',flush=True)


def improve_dimensions(primal, last_lp=104):
    verify_bch_premise()
    dual = dimension_caps(256,128,30,last_lp=last_lp)
    improved = complement_caps(primal,dual,128)
    print('Dual-complement shortening improved',sum(a>b for a,b in zip(primal,improved)),
          'lengths;',[(v,primal[v],improved[v]) for v in (128,144,160,176,192,200,224,240,256)],flush=True)
    return improved


def self_test():
    checks = 0
    for basis,n in (([],4),([15],4),([1,2,4],4),([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8)):
        code = span(basis); k = gf2_rank(basis)
        dual = [w for w in range(1<<n) if all((w&r).bit_count()%2 == 0 for r in basis)]
        caps,dcaps = [0]*(n+1),[0]*(n+1)
        for s in range(1<<n):
            v = s.bit_count(); outside = ((1<<n)-1)^s
            a = sum(w&~s == 0 for w in code).bit_length()-1
            b = sum(w&~outside == 0 for w in dual).bit_length()-1
            assert a == v-(n-k)+b
            caps[v] = max(caps[v],a); dcaps[n-v] = max(dcaps[n-v],b)
            checks += 1
        assert complement_caps(caps,dcaps,k) == caps
        assert complement_caps([min(k,v) for v in range(n+1)],dcaps,k) == caps
    print('Exact shortening/complement identities:',checks,'passed',flush=True)


if __name__ == '__main__':
    self_test()
    improve_dimensions(dimension_caps())
