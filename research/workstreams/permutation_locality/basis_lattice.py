"""Exact support-CDF improvement using attainable basis weights and orderings."""
from fractions import Fraction
from itertools import accumulate,product
from math import factorial,prod,log2
from flint import fmpz_poly
from joint_support import gf2_rank,span,griesmer,joint_bound


def power_data(spectrum,h):
    coefficients=[int(v) for v in (fmpz_poly([0]+spectrum[1:])**h).coeffs()]
    return list(accumulate(coefficients)),[w for w,c in enumerate(coefficients) if c]


def bound(spectrum,g,h,u,data=None):
    d=next(w for w,c in enumerate(spectrum) if w and c)
    if h==1:
        return ((1<<g)-1)*sum(spectrum[1:u+1])
    if u<griesmer(d,h):
        return 0
    cumulative,possible=power_data(spectrum,h) if data is None else data
    mean=Fraction(h*(1<<(h-1))*u,(1<<h)-1)
    minimum=h*d
    if mean<minimum:
        return 0
    bases=prod((1<<h)-(1<<j) for j in range(h))
    embeddings=prod((1<<g)-(1<<j) for j in range(h))
    orderings=factorial(h)
    best=None
    for index,threshold in enumerate(possible):
        if best is not None and embeddings*cumulative[threshold]//bases>=best:
            break
        if index+1==len(possible):
            good=bases
        else:
            next_weight=possible[index+1]
            if next_weight<=mean:
                continue
            # Every cheap unordered basis contributes h! ordered bases.
            lower=Fraction(bases)*(next_weight-mean)/(next_weight-minimum)/orderings
            good=orderings*(-(-lower.numerator//lower.denominator))
        assert 0<good<=bases
        candidate=embeddings*cumulative[threshold]//good
        best=candidate if best is None else min(best,candidate)
    assert best is not None
    return best


def improve_caps(caps,spectrum):
    result=[row[:] for row in caps]
    for h,row in enumerate(result,1):
        data=power_data(spectrum,h)
        for u in range(len(row)):
            row[u]=min(row[u],bound(spectrum,len(caps),h,u,data))
        for u in range(len(row)-2,-1,-1):
            row[u]=min(row[u],row[u+1])
    return result


def self_test():
    checks=0
    for rows,n,g in (([1,2,4],3,4),([15,51,85],7,4),([15,3],4,4),([0x97,0x4b,0x2d,0x1e],8,3)):
        words=span(rows)
        spectrum=[0]*(n+1)
        for w in words:
            spectrum[w.bit_count()]+=1
        shells=[[0]*(n+1) for _ in range(g+1)]
        for values in product(words,repeat=g):
            union=0
            for value in values:
                union|=value
            shells[gf2_rank(values)][union.bit_count()]+=1
        for h in range(1,min(len(rows),g)+1):
            data=power_data(spectrum,h)
            cumulative=list(accumulate(shells[h]))
            for u,actual in enumerate(cumulative):
                upper=bound(spectrum,g,h,u,data)
                assert actual<=upper
                checks+=1
    print('Attainable-basis-weight and h! rounding tests:',checks,'exact small-code inequalities',flush=True)


if __name__=='__main__':
    from bch_joint_support import authenticated_caps,support_caps
    from shortened_bound import dimension_caps
    self_test()
    spectrum=authenticated_caps()
    old=support_caps(spectrum,g=4,dimensions=dimension_caps())
    new=improve_caps(old,spectrum)
    for h in (2,3,4):
        print('rank',h,'support, old log2, new log2')
        for u in (67,68,70,71,72,74,80,96,128,192,256):
            if new[h-1][u]:
                print(u,log2(old[h-1][u]),log2(new[h-1][u]))
