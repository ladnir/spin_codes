"""Exact count bounds retaining the weight of the XOR of four outer rows.

The XOR support is exactly the set of odd-weight columns. For a nonzero
XOR f, count subcodes with a basis extending f. This is a count component,
not an inner probability bound or a complete SPIN certificate.
"""
from fractions import Fraction as Q
from itertools import product
from math import factorial,prod,log2,comb
from functools import lru_cache
from basis_lattice import power_data
from bch_joint_support import rank_total
from joint_support import span,gf2_rank,griesmer


def assignments(h,g=4):
    """Spanning tuples with sum zero, and with sum a specified nonzero f."""
    total=rank_total(h,g,h)
    zero=rank_total(h,g-1,h) if h<g else 0
    nonzero,remainder=divmod(total-zero,(1<<h)-1)
    assert not remainder
    return total,zero,nonzero


@lru_cache(maxsize=16)
def powers(spectrum,h):
    return power_data(list(spectrum),h)


def pointed_bound(spectrum,h,u,v,g=4):
    """Count rank-h spanning g-tuples, support <=u, nonzero XOR weight v.

    Ordinary spectrum coefficients may be upper bounds. The argument
    averages ordered bases (f,b_1,...,b_{h-1}) within each actual subcode.
    """
    assert 1<=h<=g and 0<v<len(spectrum)
    if v>u or not spectrum[v]:return 0
    d=next(w for w,c in enumerate(spectrum) if w and c)
    if u<griesmer(d,h):return 0
    multiplicity=assignments(h,g)[2]
    if h==1:return multiplicity*spectrum[v]
    r=h-1
    mean=Q(r*((1<<(h-1))*u-v),(1<<h)-2)
    minimum=r*d
    if mean<minimum:return 0
    bases=prod((1<<h)-(1<<i) for i in range(1,h))
    orderings=factorial(r)
    cumulative,possible=powers(tuple(spectrum),r)
    best=None
    for i,threshold in enumerate(possible):
        numerator=multiplicity*spectrum[v]*cumulative[threshold]
        if best is not None and numerator//bases>=best:break
        if i+1==len(possible):good=bases
        else:
            next_weight=possible[i+1]
            if next_weight<=mean:continue
            lower=Q(bases)*(next_weight-mean)/(next_weight-minimum)/orderings
            good=orderings*(-(-lower.numerator//lower.denominator))
        assert 0<good<=bases
        candidate=numerator//good
        best=candidate if best is None else min(best,candidate)
    assert best is not None
    return best


def marginal(spectrum,caps,u,g=4,dimensions=None):
    """Upper counts at each XOR weight, for union support <=u, by rank."""
    result=[]
    for h,row in enumerate(caps,1):
        total,zero,_=assignments(h,g)
        out=[row[u]*zero//total]
        out.extend(min(row[u],pointed_bound(spectrum,h,u,v,g)) for v in range(1,u+1))
        if dimensions is not None:
            from shortening_moments import gaussian
            n=len(spectrum)-1
            for v in range(1,u+1):
                for t in range(u,n+1):
                    spaces=gaussian(dimensions[t]-1,h-1) if dimensions[t]>=h else 0
                    numerator=spectrum[v]*comb(n-v,t-v)*spaces*assignments(h,g)[2]
                    out[v]=min(out[v],numerator//comb(n-u,t-u))
        result.append(out)
    return result


def weighted_support(spectrum,caps,dimensions,u,eta,rho='1'):
    """Integer upper sum eta^-V rho^-J, relaxing to the two marginals."""
    from flint import arb
    from joint_weight_caps import coupling,j_marginals
    eta=arb(eta);rho=arb(rho)
    assert eta>0 and 0<rho<=1
    rows=marginal(spectrum,caps,u,dimensions=dimensions)
    reflected=eta<1
    if reflected:rows=[list(reversed(row)) for row in rows]
    cumulative=[]
    for v in range(u+1):
        cumulative.append(sum(min(cap[u],sum(row[:v+1])) for cap,row in zip(caps,rows)))
    total=sum(cap[u] for cap in caps)
    js=j_marginals(tuple(spectrum),tuple(tuple(row) for row in caps),tuple(dimensions),True)[u]
    pairs=coupling(cumulative,js,total)
    value=sum((count*eta**(-(u-v) if reflected else -v)*rho**(-j) for v,j,count in pairs),arb(0))
    return int(value.upper().ceil().unique_fmpz())


def self_test():
    from collections import Counter
    checks=weighted_checks=mean_checks=0
    for basis,n in (([1,2,4],4),([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8)):
        words=span(basis)
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        # Independently check the exact mean used by the basis argument.
        h=len(basis);whole=0
        for w in words:whole|=w
        expected_bases=prod((1<<h)-(1<<i) for i in range(1,h))
        for f in words[1:]:
            count=weight_sum=0
            for extension in product(words[1:],repeat=h-1):
                if gf2_rank((f,)+extension)!=h:continue
                count+=1;weight_sum+=sum(w.bit_count() for w in extension)
            mean=Q((h-1)*((1<<(h-1))*whole.bit_count()-f.bit_count()),(1<<h)-2)
            assert count==expected_bases and weight_sum==count*mean
            mean_checks+=1
        observed=Counter();shells=[[0]*(n+1) for _ in range(4)]
        for xs in product(words,repeat=4):
            h=gf2_rank(xs)
            if not h:continue
            u=(xs[0]|xs[1]|xs[2]|xs[3]).bit_count()
            v=(xs[0]^xs[1]^xs[2]^xs[3]).bit_count()
            observed[h,u,v]+=1;shells[h-1][u]+=1
        caps=[[sum(row[:u+1]) for u in range(n+1)] for row in shells]
        for inflated in (False,True):
            spec=[c*(2 if inflated and w else 1) for w,c in enumerate(spectrum)]
            for u in range(n+1):
                upper=marginal(spec,caps,u)
                for h in range(1,5):
                    for v in range(u+1):
                        actual=sum(observed[h,x,v] for x in range(u+1))
                        assert actual<=upper[h-1][v],(basis,h,u,v,actual,upper[h-1][v])
                        checks+=1
                    # Sum-zero is an exact fraction within every subcode.
                    total,zero,_=assignments(h)
                    assert sum(observed[h,x,0] for x in range(u+1))*total==caps[h-1][u]*zero
        dimensions=[]
        for u in range(n+1):
            dimensions.append(max((sum(w&~mask==0 for w in words).bit_length()-1
                                   for mask in range(1<<n) if mask.bit_count()==u),default=0))
        for inflated in (False,True):
            spec=[c*(2 if inflated and w else 1) for w,c in enumerate(spectrum)]
            for u in range(n+1):
                upper=marginal(spec,caps,u,dimensions=dimensions)
                for h in range(1,5):
                    for v in range(u+1):
                        assert sum(observed[h,x,v] for x in range(u+1))<=upper[h-1][v]
                        checks+=1
        for eta,rho in (('1','1'),('1.1','.75'),('1.5','.5'),('.9','.75'),('.75','1')):
            for u in range(n+1):
                actual=Q(0)
                for xs in product(words,repeat=4):
                    if not any(xs):continue
                    if (xs[0]|xs[1]|xs[2]|xs[3]).bit_count()>u:continue
                    v=(xs[0]^xs[1]^xs[2]^xs[3]).bit_count()
                    j=(xs[0]&xs[1]&xs[2]&xs[3]).bit_count()
                    actual+=Q(eta)**(-v)*Q(rho)**(-j)
                assert actual<=weighted_support(spectrum,caps,dimensions,u,eta,rho)
                weighted_checks+=1
    print('Odd-column bounds:',checks,'exact/inflated marginal checks;',mean_checks,
          'basis-mean identities;',weighted_checks,'joint weighted inequalities passed',flush=True)


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test-only',action='store_true')
    args=parser.parse_args();self_test()
    if args.test_only:return
    from bch_joint_support import authenticated_caps,support_caps
    from shortened_bound import dimension_caps
    from basis_lattice import improve_caps
    from shortening_moments import improve
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    caps,_=improve(improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum),dimensions)
    for u in (80,96,128,160):
        rows=marginal(spectrum,caps,u,dimensions=dimensions)
        print('Support',u,'old total log2',log2(sum(row[u] for row in caps)),flush=True)
        for v in (0,38,48,64,80,96,128):
            if v>u:continue
            values=[row[v] for row in rows];total=sum(values)
            print('  XOR weight',v,'log2 upper',log2(total) if total else None,
                  'rank4',log2(values[3]) if values[3] else None,flush=True)
        total=sum(min(cap[u],sum(row)) for cap,row in zip(caps,rows))
        print('  summed marginal log2',log2(total),flush=True)


if __name__=='__main__':main()
