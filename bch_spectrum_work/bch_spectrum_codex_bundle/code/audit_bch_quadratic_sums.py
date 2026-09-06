"""Exact signed weight sums for P, Q, and a nonzero Q-coset.

For even binary words, wt(x)/2 modulo two is quadratic with polar form x dot y.
Split off hyperbolic pairs, then evaluate the remaining linear radical.
"""
import itertools
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from bch_quotient import generator_polynomial
from audit_bch_wambach_shortening import rank
from verify_locator_reduction import wambach_support
from affine_wambach import gf_pow
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_quadratic_weight_sums.json'


def phase(word,offset):
    assert word.bit_count()%2==offset.bit_count()%2==0
    return ((word.bit_count()//2)+(word&offset).bit_count())%2


def signed_sum(basis,offset=0):
    original=len(basis)
    basis=basis.copy()
    sign=-1 if (offset.bit_count()//2)%2 else 1
    pairs=0
    while True:
        pair=next(((i,j) for i in range(len(basis)) for j in range(i+1,len(basis))
                   if (basis[i]&basis[j]).bit_count()%2),None)
        if pair is None:
            break
        i,j=pair
        u,v=basis[i],basis[j]
        sign*=(-1 if phase(u,offset)*phase(v,offset) else 1)
        next_basis=[]
        for k,x in enumerate(basis):
            if k in (i,j):
                continue
            y=x
            if (x&v).bit_count()%2:
                y^=u
            if (x&u).bit_count()%2:
                y^=v
            assert (y&u).bit_count()%2==(y&v).bit_count()%2==0
            next_basis.append(y)
        basis=next_basis
        pairs+=1
    assert 2*pairs+len(basis)==original and rank(basis)==len(basis)
    nonzero=next((x for x in basis if phase(x,offset)),None)
    value=0 if nonzero is not None else sign*(1<<(pairs+len(basis)))
    return dict(signed_sum=str(value),bilinear_rank=2*pairs,radical_dimension=len(basis),
                radical_phase_witness_hex=hex(nonzero) if nonzero is not None else None,
                zero_sum=nonzero is not None)


def rows(delta,dimension):
    g=generator_polynomial(delta)
    result=[]
    for i in range(dimension):
        x=g<<i
        result.append(x|((x.bit_count()%2)<<255))
    assert rank(result)==dimension
    return result


def toy():
    checks=0
    for n in (4,6,8):
        basis=[(1<<j)|(1<<(n-1)) for j in range(n-1)]
        for dimension in range(n):
            b=basis[:dimension]
            for offset in (0,3,((1<<n)-1)):
                exact=0
                for mask in range(1<<dimension):
                    word=offset
                    for j,x in enumerate(b):
                        if (mask>>j)&1:
                            word^=x
                    exact+=(-1)**(word.bit_count()//2)
                assert int(signed_sum(b,offset)['signed_sum'])==exact
                checks+=1
    return checks


def build():
    coordinates={gf_pow(2,j):j for j in range(255)}
    coordinates[0]=255
    offset=sum(1<<coordinates[x] for x in wambach_support()+[0])
    p,q=rows(37,131),rows(39,123)
    results={'P':signed_sum(p),'Q':signed_sum(q),'H':signed_sum(q,offset)}
    assert int(results['P']['signed_sum'])==int(results['Q']['signed_sum'])+255*int(results['H']['signed_sum'])
    return dict(classification='Exact quadratic-form signed weight enumeration',
                meaning='signed_sum = number of words of weight 0 mod 4 minus number of words of weight 2 mod 4',
                results=results,exhaustive_toy_cases=toy(),new_A38_cap_proved=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/bch_quotient.py',
                  ROOT/'code/audit_bch_wambach_shortening.py',ROOT/'code/verify_locator_reduction.py',ROOT/'code/affine_wambach.py')})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
