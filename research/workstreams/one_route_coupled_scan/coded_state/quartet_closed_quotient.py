"""Exact seven-type quotient of the twelve pending-input types.

No age maximum, population truncation, or changed construction. Exhaustive
integer comparison of every local row, plus exact small regional replay.
"""
from collections import defaultdict
from functools import lru_cache
import itertools
import math
from quartet_closed_transfer import TYPES,local,terminal


# Z:U0; V0,V1:unselected allowed runs; B:must see a non-full quartet;
# L0,L1,L2:selected run stages, with L2 mature enough to stop.
NAMES=('Z','V0','V1','B','L0','L1','L2')
REPRESENTATIVES=(0,1,4,7,2,11,14)
PROJECT={0:0,3:0,6:0,15:0,1:1,4:2,7:3,16:3,2:4,11:5,14:6,17:6}
assert set(PROJECT)==set(TYPES)


def project_row(hist):
    out=defaultdict(int)
    for (dest,u,v,a,d),ways in local(hist):
        out[tuple(sorted(PROJECT[t] for t in dest)),u,v,a,d]+=ways
    return dict(out)


@lru_cache(maxsize=210)
def kernel(hist):
    return tuple(project_row(tuple(sorted(REPRESENTATIVES[t] for t in hist))).items())


def accepted(hist):
    return all(t in (0,4,5,6) for t in hist)


def coefficient_table(logb=2.5):
    """Exact quartet coefficients, evaluated in FLOAT for proposals.

    This has the same interface as joint.caps, but keeps all seven types.
    There is no coefficientwise maximization over hidden run ages.
    """
    import numpy as np
    out={}
    for hist in itertools.combinations_with_replacement(range(7),4):
        old=tuple(hist.count(t) for t in range(7))
        for (dest,u,v,a,d),ways in kernel(hist):
            new=tuple(dest.count(t) for t in range(7))
            for paid in range(v+1):
                key=old,new,u+paid
                term=math.log2(ways*math.comb(v,paid))+logb*((707*a-1905*d)/2560-6*paid/5)
                out[key]=np.logaddexp2(out.get(key,-math.inf),term)
    return out


def regional():
    from quartet_closed_paid_sets import tests as closed_tests
    expected=closed_tests();states={((0,0,0,0),0,0,0):1};peak=1
    for step in range(4):
        after=defaultdict(int)
        for (hist,running,paid,score),count in states.items():
            for (dest,u,v,a,d),ways in kernel(hist):
                for c in range(min(v,2-running-u)+1):
                    quota=running+u+c
                    if step%2 and quota!=2:continue
                    after[dest,0 if step%2 else quota,paid+c,score+2*a-5*d]+=count*ways*math.comb(v,c)
        states=dict(after);peak=max(peak,len(states))
    result=defaultdict(int)
    for (hist,quota,paid,score),count in states.items():
        if accepted(hist):result[paid,score]+=count
    assert dict(result)==expected
    print(dict(seven_type_regional_coefficients=len(result),peak_states=peak),flush=True)


def tests():
    checked=0
    for hist in itertools.combinations_with_replacement(TYPES,4):
        reduced=tuple(sorted(PROJECT[t] for t in hist))
        assert project_row(hist)==dict(kernel(reduced));checked+=1
    assert all(terminal((t,))==accepted((PROJECT[t],)) for t in TYPES)
    identities=0
    for hist in itertools.combinations_with_replacement(range(7),4):
        for (dest,u,v,a,d),ways in kernel(hist):
            assert dest.count(5)==hist.count(4)
            assert dest.count(6)==hist.count(5)+hist.count(6)-d
            assert dest.count(3)>=d
            assert a-sum(t>=4 for t in hist)==dest.count(4)-d
            identities+=1
    print(dict(exact_local_quotient_rows=checked,quotient_rows=kernel.cache_info().currsize,
               original_types=len(TYPES),exact_types=7,terminal_types_checked=len(TYPES),
               exact_age_flow_rows=identities),flush=True)
    regional()


def diagnostic():
    """Selected stationary-population bounds; not coverage or a certificate."""
    from quartet_closed_joint import bound
    table=coefficient_table()
    profiles=((104,156,60,64,40,40,48),(104,156,60,64,32,32,64),
              (128,140,52,64,40,40,48),(104,64,24,64,48,48,160),
              (104,176,72,32,16,16,96))
    for pop in profiles:
        z,v0,v1,b,l0,l1,l2=pop
        assert sum(pop)==512 and l0==l1 and l0<=min(b,l2)
        values=[bound(table,pop,pop,i) for i in (131,132)]
        score=12*values[0][0]+4*values[1][0]+math.log2(math.comb(16,4))-math.log2(math.comb(8192,2100))
        print(dict(population=pop,conditional_region_log2=float(score),
                   maximum_fit_residual=max(v[2] for v in values),full_certificate=False),flush=True)


if __name__=='__main__':
    import sys
    tests()
    if '--diagnostic' in sys.argv:diagnostic()
