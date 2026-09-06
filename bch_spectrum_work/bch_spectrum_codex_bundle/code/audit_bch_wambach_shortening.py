"""Exact projection rank and exhaustive small dual-support exclusion.

Use one explicit minimum word to shorten Q on its 38-coordinate support.
No spectrum cap is inferred. NumPy performs only unsigned integer XOR,
lexicographic sorting, and comparisons; every subset of size <=6 is covered.
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from affine_wambach import gf_pow
from bch_quotient import generator_polynomial,binary_poly_divmod
from verify_locator_reduction import wambach_support
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_wambach_shortening.json'
DTYPE=np.dtype([('hi','<u8'),('lo','<u8')])


def rank(words):
    pivots={}
    for word in words:
        while word:
            bit=word.bit_length()-1
            if bit not in pivots:
                pivots[bit]=word
                break
            word^=pivots[bit]
    return len(pivots)


def projection():
    field_support=wambach_support()+[0]
    field_to_coordinate={gf_pow(2,j):j for j in range(255)}
    field_to_coordinate[0]=255
    support=sorted(field_to_coordinate[x] for x in field_support)
    mask=sum(1<<j for j in support)
    p,q=generator_polynomial(37),generator_polynomial(39)
    punctured=mask&((1<<255)-1)
    assert mask.bit_count()==38 and mask.bit_count()%2==0
    assert binary_poly_divmod(punctured,p)[1]==0
    assert binary_poly_divmod(punctured,q)[1]!=0
    rows=[]
    for shift in range(123):
        word=q<<shift
        word|=(word.bit_count()%2)<<255
        rows.append(word)
    assert rank(rows)==123
    pivots={}
    kernel=[]
    for original in rows:
        word=original
        projected=sum(((word>>j)&1)<<i for i,j in enumerate(support))
        while projected:
            bit=projected.bit_length()-1
            if bit not in pivots:
                pivots[bit]=(projected,word)
                break
            v,full=pivots[bit]
            projected^=v
            word^=full
        if not projected:
            assert word and word&mask==0
            kernel.append(word)
    assert len(pivots)==38 and len(kernel)==85 and rank(kernel)==85
    assert all(binary_poly_divmod(v&((1<<255)-1),q)[1]==0 and v.bit_count()%2==0 for v in kernel)
    outside=[j for j in range(256) if j not in support]
    columns=[sum(((word>>j)&1)<<i for i,word in enumerate(kernel)) for j in outside]
    return dict(support_coordinates=support,support_field_elements=field_support,
                Q_generator_hex=hex(q),projection_rank=38,shortened_dimension=85,
                kernel_basis_hex=[hex(v) for v in kernel],outside_coordinates=outside,
                shortened_columns_hex=[hex(v) for v in columns]),columns


def packed(values):
    result=np.empty(len(values),dtype=DTYPE)
    result['hi']=[v>>64 for v in values]
    result['lo']=[v&((1<<64)-1) for v in values]
    return result


def membership(left,right):
    # right is sorted; appending no sentinel avoids altering the set.
    indices=np.searchsorted(right,left)
    valid=indices<len(right)
    return bool(np.any(left[valid]==right[indices[valid]]))


def check_columns(columns):
    n=len(columns)
    assert n==218 and len(set(columns))==n and 0 not in columns
    single=packed(columns)
    single.sort(order=('hi','lo'))
    pairs=packed([columns[i]^columns[j] for i in range(n) for j in range(i+1,n)])
    pairs.sort(order=('hi','lo'))
    assert len(pairs)==math.comb(n,2)
    assert not np.any(pairs[1:]==pairs[:-1])
    assert not membership(pairs,single)
    triples=np.empty(math.comb(n,3),dtype=DTYPE)
    base=packed(columns)
    offset=0
    for i in range(n):
        for j in range(i+1,n):
            width=n-j-1
            triples['hi'][offset:offset+width]=base['hi'][i]^base['hi'][j]^base['hi'][j+1:]
            triples['lo'][offset:offset+width]=base['lo'][i]^base['lo'][j]^base['lo'][j+1:]
            offset+=width
    assert offset==math.comb(n,3)
    triples.sort(order=('hi','lo'))
    assert not np.any((triples['hi']==0)&(triples['lo']==0))
    assert not membership(triples,pairs)
    assert not np.any(triples[1:]==triples[:-1])
    # Explicitly verify the sorting order, which is required by searchsorted.
    for values in (single,pairs,triples):
        assert np.all((values['hi'][1:]>values['hi'][:-1]) |
                      ((values['hi'][1:]==values['hi'][:-1]) & (values['lo'][1:]>values['lo'][:-1])))
    return dict(singletons_checked=n,pairs_checked=len(pairs),triples_checked=len(triples),
                zero_singletons_or_triples=False,duplicate_singletons_pairs_or_triples=False,
                pair_singleton_or_triple_pair_matches=False,shortened_dual_minimum_at_least=7)


def build():
    data,columns=projection()
    checks=check_columns(columns)
    paths=[Path(__file__),ROOT/'code/affine_wambach.py',ROOT/'code/bch_quotient.py',ROOT/'code/verify_locator_reduction.py']
    return dict(classification='Exact code-specific projection and exhaustive shortened-dual support exclusion',
                Q_parameters=[256,123,40],shortened_parameters=[218,85,40],
                projection_fiber_size=str(1<<85),projection=data,checks=checks,
                new_split_dual_zero_region='All (u,v) with 0<=u<=38 and 0<=v<=6, except (0,0)',
                A38_upper_cap_proved=False,original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:result[k] for k in ('classification','shortened_parameters','projection_fiber_size',
                                          'checks','new_split_dual_zero_region','A38_upper_cap_proved')},indent=2))
