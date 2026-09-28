"""Exact zero-feedback counts for one-column groups sharing an IMT epoch.

Enumerates window subsets and their kernel vectors, not all input masks.
This is a local activation check, not an output-distance certificate.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations,permutations
from math import comb,factorial,prod
import argparse

from group_moment import maps
from joint_support import gf2_rank


def kernel_basis(columns):
    basis={}
    kernel=[]
    for j,column in enumerate(columns):
        value,word=column,1<<j
        while value:
            pivot=value.bit_length()-1
            if pivot not in basis:
                basis[pivot]=value,word
                break
            old,mask=basis[pivot]
            value^=old
            word^=mask
        if not value:
            kernel.append(word)
    assert len(basis)+len(kernel)==len(columns)
    return len(basis),kernel


def vectors(basis):
    result=[0]
    for word in basis:
        result += [old^word for old in result]
    return result


def self_test():
    for columns in ([1,2,4,8],[1,2,3,1,2,3],[0,1,1,2,3,4,5,6]):
        rank,basis=kernel_basis(columns)
        actual=[]
        for mask in range(1<<len(columns)):
            image=0
            for j,column in enumerate(columns):
                if (mask>>j)&1:
                    image^=column
            if not image:
                actual.append(mask)
        assert set(vectors(basis))==set(actual)
        assert rank==gf2_rank(columns)
    print('Kernel basis enumeration agrees with direct exhaustive maps',flush=True)


def census(columns,q):
    ranks=Counter()
    zero=Counter()
    count=0
    for windows in combinations(range(32),q):
        local=[column for j in windows for column in columns[4*j:4*j+4]]
        rank,basis=kernel_basis(local)
        ranks[rank]+=1
        for word in vectors(basis)[1:]:
            weights=tuple(((word>>(4*j))&15).bit_count() for j in range(q))
            if not all(weights):
                continue
            count+=1
            # Every assignment of these distinct physical windows to labeled
            # groups counts once, even when some input weights coincide.
            zero.update(permutations(weights))
    assert sum(ranks.values())==comb(32,q)
    assert sum(zero.values())==factorial(q)*count
    denominator=prod(range(32-q+1,33))
    worst=max(((Fraction(n,denominator*prod(comb(4,w) for w in shape)),shape)
               for shape,n in zero.items()),default=(Fraction(0),None))
    print('q',q,'rank histogram',dict(sorted(ranks.items())),
          'fully active kernel masks',count,'worst conditional zero probability',worst,flush=True)
    return ranks,zero


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-groups',type=int,default=4,choices=(1,2,3,4))
    args=parser.parse_args()
    self_test()
    _,columns,_=maps()
    for q in range(1,args.max_groups+1):
        census(columns,q)
    print('Conditional on sharing one epoch; not a full-code failure bound.')


if __name__=='__main__':
    main()
