"""Joint total-weight/all-one weighting using exact marginal count caps.

For a>=1 and 0<rho<=1, bound sum a^-W rho^-J over four-word groups.
Pack the W marginal downward and the J marginal upward, then pair their
largest weights. This is a relaxed coupling bound, not a joint enumerator.
"""
from fractions import Fraction as Q
from collections import Counter
from functools import lru_cache
from itertools import product
from flint import arb,fmpz_poly
from occupancy_allones import full_column_caps,prefix_flag_caps
from occupancy_weighted import allocate


def coupling(weight_cdf,j_caps,total):
    assert all(a<=b for a,b in zip(weight_cdf,weight_cdf[1:]))
    assert total>=0 and all(v>=0 for v in j_caps)
    total=min(total,weight_cdf[-1],sum(j_caps))
    previous=0;ws=[]
    for w,cap in enumerate(weight_cdf):
        current=min(total,cap)
        if current>previous:ws.append([w,current-previous])
        previous=current
    js=[list(x) for x in allocate(j_caps,total,Q(1,2))]
    result=[];i=k=0
    while i<len(ws) and k<len(js):
        count=min(ws[i][1],js[k][1])
        result.append((ws[i][0],js[k][0],count))
        ws[i][1]-=count;js[k][1]-=count
        if not ws[i][1]:i+=1
        if not js[k][1]:k+=1
    assert sum(c for _,_,c in result)==total
    return result


@lru_cache(maxsize=4)
def j_marginals(spectrum,caps,dimensions,prefix_flags):
    n=len(spectrum)-1
    totals=[sum(row[u] for row in caps) for u in range(n+1)]
    prefix=prefix_flag_caps(spectrum,caps,dimensions) if prefix_flags else None
    shells=[0]*(n+1);result=[(0,)]
    for u in range(1,n+1):
        current=full_column_caps(spectrum,caps,dimensions,u,totals[u])
        for j,value in enumerate(current):shells[j]+=value
        result.append(tuple([totals[u]]+[min(shells[j],prefix[u][j]) if prefix else shells[j]
                                         for j in range(1,u+1)]))
    return tuple(result)


def caps(spectrum,support_caps,dimensions,a,rho,prefix_flags=True):
    a=arb(a);rho=arb(rho)
    assert a>=1 and 0<rho<=1
    n=len(spectrum)-1
    counts=[sum(row[u] for row in support_caps) for u in range(n+1)]
    coefficients=[int(x) for x in fmpz_poly(spectrum)**4]
    coefficients += [0]*(4*n+1-len(coefficients))
    coefficients[0]-=1
    cumulative=[];total=0
    for value in coefficients:total+=value;cumulative.append(total)
    j_caps=j_marginals(tuple(spectrum),tuple(tuple(row) for row in support_caps),tuple(dimensions),prefix_flags)
    wp=[a**(-w) for w in range(4*n+1)]
    jp=[rho**(-j) for j in range(n+1)]
    result=[]
    for u,limit in enumerate(counts):
        pairs=coupling(cumulative[:4*u+1],j_caps[u],limit)
        bound=sum((count*wp[w]*jp[j] for w,j,count in pairs),arb(0))
        result.append(int(bound.upper().ceil().unique_fmpz()))
    assert all(x<=y for x,y in zip(result,result[1:]))
    return result


def self_test():
    checked=0
    # Exhaust all 2-by-3 count tables, including inflated marginal caps.
    for values in product(range(3),repeat=6):
        table=[values[:3],values[3:]]
        wm=[sum(row) for row in table]
        jm=[sum(table[w][j] for w in range(2)) for j in range(3)]
        for inflation in (0,1):
            cdf=[wm[0]+inflation,sum(wm)+2*inflation]
            pairs=coupling(cdf,[v+inflation for v in jm],sum(wm)+inflation)
            for a,rho in ((Q(1),Q(1)),(Q(11,10),Q(3,4)),(Q(2),Q(1,2))):
                actual=sum(table[w][j]*a**(-w)*rho**(-j) for w in range(2) for j in range(3))
                upper=sum(c*a**(-w)*rho**(-j) for w,j,c in pairs)
                assert actual<=upper
                checked+=1
    print('Joint weight/all-one rearrangement:',checked,'exact table checks passed',flush=True)
    from joint_support import span,gf2_rank
    code_checks=0
    for basis,n in (([1,2,4],3),([15,51,85],7),([3,5],4)):
        words=span(basis)
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        shells=[[0]*(n+1) for _ in range(4)];observed=Counter()
        for xs in product(words,repeat=4):
            if not any(xs):continue
            u=(xs[0]|xs[1]|xs[2]|xs[3]).bit_count()
            w=sum(x.bit_count() for x in xs)
            j=(xs[0]&xs[1]&xs[2]&xs[3]).bit_count()
            shells[gf2_rank(xs)-1][u]+=1;observed[u,w,j]+=1
        support_caps=[[sum(row[:u+1]) for u in range(n+1)] for row in shells]
        dimensions=[0]*(n+1)
        for mask in range(1<<n):
            k=sum(w & ~mask==0 for w in words).bit_length()-1
            dimensions[mask.bit_count()]=max(dimensions[mask.bit_count()],k)
        for a,rho in (('1','1'),('1.03','.75'),('1.2','.5')):
            upper=caps(spectrum,support_caps,dimensions,a,rho)
            for u in range(n+1):
                actual=sum(count*Q(a)**(-w)*Q(rho)**(-j)
                           for (v,w,j),count in observed.items() if v<=u)
                assert actual<=upper[u]
                code_checks+=1
    print('Joint weight/all-one complete outer pipeline:',code_checks,'exact small-code inequalities',flush=True)


if __name__=='__main__':
    self_test()
