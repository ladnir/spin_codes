"""Independent small/exact checks for the two-group proof components."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb, exp

import numpy as np
from flint import arb, ctx

from group_moment import maps
from two_group_moment import collision_census, self_test
from two_group_screen import log_power_moment
from two_group_verify import collision_transfers
from two_column_moment import census as single_census
from two_group_moment import worst_regions
from two_group_weight import weighted_regions


def test_census_direct(data):
    images,columns,_ = maps()
    bundle_columns=len(next(iter(data[1]))[0])
    width=4*bundle_columns
    selections=(((0,1),(0,1)),((4,4),(4,4)),((3,3),(3,3))) if bundle_columns==2 else (((1,),(1,)),((4,),(4,)),((3,),(3,)))
    for a,b in selections:
        masks = {s: [mask for mask in range(1,1<<width)
                     if tuple(sorted(((mask>>(4*j))&15).bit_count() for j in range(bundle_columns))) == s]
                 for s in (a,b)}
        atoms = Counter()
        cancel = {v:Counter() for v in data[0]}
        for j in range(128//width):
            for k in range(128//width):
                if j == k:
                    continue
                for ma,mb in product(masks[a],masks[b]):
                    word = (ma << (width*j)) | (mb << (width*k))
                    q = 0
                    for bit in range(128):
                        if (word >> bit)&1:
                            q ^= columns[bit]
                    atoms[q] += 1
                    if q:
                        cancel[images[q].bit_count()][(images[q]^word).bit_count()] += 1
        choices,zero,maximum,expected = data[1][a,b]
        assert sum(atoms.values()) == choices and atoms[0] == zero
        assert max(count for q,count in atoms.items() if q) == maximum
        assert cancel == expected
    print('Direct 128-bit scalar census agrees for three shape pairs; columns',bundle_columns,flush=True)


def test_coefficients():
    # Two macroregions, each with two support coordinates per group.
    # Enumerate support sets directly, without a polynomial recurrence.
    region = {(a,b): np.array([[1+a,1+b],[a+b,2+a*b]],dtype=np.int64)
              for a in range(3) for b in range(3)}
    coefficients = [[0]*5 for _ in range(5)]
    for a,b in product(range(16),repeat=2):
        row = np.array([1,0],dtype=np.int64)
        for j in range(2):
            row = row@region[((a>>(2*j))&3).bit_count(),((b>>(2*j))&3).bit_count()]
        coefficients[a.bit_count()][b.bit_count()] += int(row.sum())
    for p,q in ((F(1,3),F(2,5)),(F(3,4),F(1,2))):
        bernoulli = sum(coefficients[u][v]*p**u*(1-p)**(4-u)*q**v*(1-q)**(4-v)
                        for u in range(5) for v in range(5))
        for lo in range(5):
            for hi in range(lo,5):
                min_p = min(F(comb(4,u))*p**u*(1-p)**(4-u) for u in (lo,hi))
                assert min_p == min(F(comb(4,u))*p**u*(1-p)**(4-u) for u in range(lo,hi+1))
                for low in range(5):
                    for high in range(low,5):
                        min_q = min(F(comb(4,v))*q**v*(1-q)**(4-v) for v in (low,high))
                        bound = bernoulli/(min_p*min_q)
                        for u in range(lo,hi+1):
                            for v in range(low,high+1):
                                actual = F(coefficients[u][v],comb(4,u)*comb(4,v))
                                assert actual <= bound
    for n in (1,2,7,16):
        matrix = np.array([[.1,.3],[.2,.4]])
        expected = np.linalg.matrix_power(matrix,n)[0].sum()
        assert abs(exp(log_power_moment(matrix,n))/expected-1) < 1e-13
    print('Exact positive coefficient/rectangle bounds and scaled matrix powers checked',flush=True)


def main():
    self_test()
    test_coefficients()
    data = collision_census()
    test_census_direct(data)
    ctx.prec = 192
    for tilt in ('0','0.0004','0.0025'):
        transfers = collision_transfers(data,tilt)
        for (a,b),(choices,zero,maximum,cancel) in data[1].items():
            matrix = transfers[a,b]
            for i,v in enumerate(sorted(data[0]),2):
                # The exact lazy cancellation component must be bounded
                # even before adding the positive uniform-refresh term.
                exact = sum((count*(-arb(tilt)*w).exp() for w,count in cancel[v].items()),arb(0))
                exact /= 2*choices*data[0][v]
                assert matrix[i,0] >= exact
            if tilt == '0':
                assert all(sum((matrix[i,j] for j in range(7)),arb(0)).upper() >= 1 for i in range(7))
    print('Outward cancellation components and zero-tilt mass checked',flush=True)
    one=collision_census(1)
    assert len(one[1])==16 and all(row[1]==0 for row in one[1].values())
    test_census_direct(one)
    for width,pair in ((1,one),(2,data)):
        single=single_census(width)
        raw=worst_regions(single,pair,.0004,return_shapes=True)
        reference=worst_regions(single,pair,.0004)
        for statistic in ('weight','all-ones'):
            actual=weighted_regions(raw,1.,1.,statistic)
            assert actual.keys()==reference.keys()
            assert all(np.allclose(actual[key],reference[key],rtol=1e-14,atol=1e-280) for key in actual)
    print('Unweighted limits of both auxiliary statistics match both geometries',flush=True)


if __name__ == '__main__':
    main()
