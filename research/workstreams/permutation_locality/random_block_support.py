"""Exact delayed-activation bound for uniform shared coordinates, g=c=4.

This bounds one active group's completely zero state prefix, not distance.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,log2
from flint import fmpz_poly
from bch_joint_support import authenticated_caps,support_caps,weighted_cdf_bound
from canonical_census import census
from canonical_support import rank_one_zero_test
from shortened_bound import dimension_caps


def prefix_weights(n,width,length,q):
    assert n%width==0 and 0<=length<=n//width and 0<=q<=1
    a,b=q.numerator,q.denominator
    active=fmpz_poly([b]+[a*comb(width,j) for j in range(1,width+1)])
    free=fmpz_poly([comb(n-width*length,j) for j in range(n-width*length+1)])
    polynomial=active**length*free
    result=[F(int(polynomial[u]),b**length*comb(n,u)) for u in range(n+1)]
    assert result[0]==1 and all(x>=y>=0 for x,y in zip(result,result[1:]))
    return result


def self_test():
    for n,width in ((4,2),(8,2),(8,4),(12,4)):
        for length in range(n//width+1):
            for q in (F(0),F(1,3),F(1)):
                weights=prefix_weights(n,width,length,q)
                for u in range(n+1):
                    exact=sum((q**len({i//width for i in support if i<width*length})
                               for support in combinations(range(n),u)),F(0))/comb(n,u)
                    assert exact==weights[u]
    print('Uniform-support occupancy polynomial matches exhaustive subsets',flush=True)


def main():
    self_test()
    totals,zero=census()
    rank_one_zero_test(totals,zero)
    q=max(F(zero[s],v) for s,v in totals.items() if any(s))
    assert q==F(1,9216)
    caps=support_caps(authenticated_caps(),g=4,dimensions=dimension_caps())
    logarithm=lambda v: log2(v.numerator)-log2(v.denominator) if v else float('-inf')
    for length in (40,44,46,48,49,50,51,52,53,54,55,56,57,60,63):
        bounds=[2048*weighted_cdf_bound(row,prefix_weights(256,4,length,F(0) if h==1 else q))
                for h,row in enumerate(caps,1)]
        total=sum(bounds,F(0))
        print('prefix',length,'all ranks log2',round(logarithm(total),6),
              'rank log2',[round(logarithm(v),6) for v in bounds],
              'below 2^-40',total<F(1,1 << 40),flush=True)
        if length==51:
            assert 0<total<F(1,1 << 45)
            print('EXACT CHECK: one-active-group zero-prefix union < 2^-45 at 51 macroregions',flush=True)
    print('No full-distance claim: later cancellation, output weight, and multiple active groups remain open.')


if __name__=='__main__':
    main()
