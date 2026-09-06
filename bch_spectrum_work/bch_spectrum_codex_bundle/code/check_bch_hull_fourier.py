"""Exhaustive small-code checks of the signed-transform hull inequalities."""
import json
import random
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from audit_bch_hulls import echelon, hull, nullspace
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_hull_fourier_toy.json'


def span(basis):
    values=[0]
    for x in basis:
        values += [y^x for y in values]
    return set(values)


def build():
    rng=random.Random(20260905)
    checked=0
    for n in (4,6,8):
        even=[x for x in range(1<<n) if x.bit_count()%2==0]
        for trial in range(120):
            basis=list(echelon([rng.choice(even) for _ in range(1+trial%(n-1))]).values())
            e=span(basis)
            signs={x:(-1)**(x.bit_count()//2) for x in e}
            signed_sum=sum(signs.values())
            if signed_sum>=0:
                continue
            r=span(hull(basis))
            dual=span(nullspace(basis,n))
            hull_dual=span(nullspace(list(r),n))
            magnitude=-signed_sum
            assert all(x.bit_count()%4==0 for x in r)
            values={u:sum(signs[x]*((-1)**((x&u).bit_count()%2)) for x in e)
                    for u in range(1<<n)}
            assert {u for u,v in values.items() if v}==hull_dual
            assert {abs(v) for v in values.values() if v}=={magnitude}
            for j in range(n+1):
                pos={u for u,v in values.items() if u.bit_count()==j and v>0}
                neg={u for u,v in values.items() if u.bit_count()==j and v<0}
                ej={x for x in e if x.bit_count()==j}
                dj={x for x in dual if x.bit_count()==j}
                rj={x for x in r if x.bit_count()==j}
                assert (ej if j%4==2 else set()) <= pos
                assert dj | (ej if j%4==0 else set()) <= neg
                assert len(dj | (ej if j%4==0 else set())) == len(dj)+(len(ej) if j%4==0 else 0)-len(rj)
            checked+=1
    assert checked>20
    return dict(classification='Exhaustive toy verification, not a BCH spectrum certificate',
                negative_Gauss_sum_codes_checked=checked,lengths=[4,6,8],
                all_Fourier_support_magnitude_and_signed_inclusion_checks_passed=True,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/audit_bch_hulls.py')})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps(result,indent=2))
