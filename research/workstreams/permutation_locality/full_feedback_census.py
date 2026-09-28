"""Exact full feedback distributions for up to ten distinct windows.

Walsh inversion of compressed character polynomials. All counts and
intermediate integer bounds are exact; no independence approximation.
"""
import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, prod
import numpy as np

from feedback_character_census import character_polynomials
from pair_tail import hadamard
from group_moment import maps


def census(maximum=6):
    shapes,coefficients,denoms,multiplicities,inverse=character_polynomials(maximum)
    n=1<<19
    images,columns,spectrum=maps()
    assert all(c.bit_count()%2==1 for c in columns)
    parity=np.bitwise_count(np.arange(n,dtype=np.uint32))&1
    expansion=np.array([w.bit_count() for w in images],dtype=np.uint8)
    results={}
    for shape,row in zip(shapes,coefficients):
        weights=tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        # Every Walsh butterfly is bounded by the sum of absolute inputs.
        # Evaluate that bound with Python integers before the int64 WHT.
        denominator=denoms[shape]
        assert int(np.abs(row).max())<=denominator
        absolute=sum(abs(int(x))*int(m) for x,m in zip(row,multiplicities))
        assert absolute<1<<63,('Walsh integer width insufficient',weights,absolute)
        counts=hadamard(row[inverse])
        assert not np.any(counts%n)
        counts//=n
        assert int(counts.min())>=0 and int(counts.sum())==denominator
        assert not np.any(counts[parity!=(sum(weights)%2)])
        assert int(counts[0])*n==sum(int(x)*int(m) for x,m in zip(row,multiplicities))
        # Verify selected characters again from the recovered distribution.
        for character in (0,1,17,524287):
            signs=1-2*(np.bitwise_count(np.arange(n,dtype=np.uint32)&character)&1).astype(np.int64)
            assert int(counts@signs)==int(row[inverse[character]])
        classes={v:int(counts[expansion==v].sum()) for v in spectrum}
        assert sum(classes.values())+int(counts[0])==denominator
        results[weights]=(int(counts[0]),int(counts[1:].max()),denominator,classes)
        if len(weights)<=2:
            # Independent enumeration: unordered windows, all assignments of
            # the weight multiset, and every corresponding lane mask.
            direct=np.zeros(n,dtype=np.int64)
            for windows in combinations(range(32),len(weights)):
                for masks in product(range(1,16),repeat=len(weights)):
                    if tuple(sorted(m.bit_count() for m in masks))!=weights:
                        continue
                    syndrome=0
                    for window,mask in zip(windows,masks):
                        for bit in range(4):
                            if mask>>bit&1:syndrome^=columns[4*window+bit]
                    direct[syndrome]+=1
            assert np.array_equal(counts,direct),weights
    for j in range(1,maximum+1):
        rows=[(Fraction(peak,den),w,Fraction(z,den)) for w,(z,peak,den,_) in results.items() if len(w)==j]
        print('Exact full feedback: windows',j,'worst nonzero atom/shape/zero',max(rows),flush=True)
    print('All full distributions pass integer bounds, positivity, total, zero, and character checks;',
          'every one-/two-window distribution matches independent enumeration',flush=True)
    return results


def check_prior_zeros(results):
    from occupancy_model import local_data
    data=local_data(4);checked=0
    for j,rows in data[2].items():
        for weights in product(range(1,5),repeat=j):
            key=tuple(sorted(weights))
            multiplicity=prod(factorial(weights.count(b)) for b in range(1,5))
            assert results[key][0]*multiplicity==rows.get(weights,0)
            checked+=1
    print('Full-feedback vs direct-rank zero checks:',checked,'ordered shapes',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,choices=range(2,11),default=6)
    args=parser.parse_args()
    results=census(args.maximum)
    if args.maximum>=4:check_prior_zeros(results)
