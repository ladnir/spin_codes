"""Exact upper bounds for canonical-block supports and zero prefixes.

This is a necessary-event diagnostic, not a full-distance certificate.
All accepted LP witnesses, counts, and probability bounds use exact rationals.
"""
from fractions import Fraction as F
from itertools import accumulate, combinations, product
from math import comb, log2

from bch_joint_support import (authenticated_caps, prefix_probabilities,
    rank_total, support_caps, weighted_cdf_bound)
from canonical_census import census
from joint_support import gf2_rank, span
from shortened_bound import dimension_caps


def block_caps(bit_caps, dimensions, n=256, width=4, g=4, k=128):
    """CDF of rank-h g-tuples whose union touches at most u fixed blocks."""
    assert n % width == 0
    blocks = n // width
    result = []
    for h, row in enumerate(bit_caps, 1):
        caps = []
        for u in range(blocks+1):
            # A tuple in u blocks is supported on width*u bits. Alternatively
            # union-bound over all u-block shortened codes, counting ranks.
            dim = dimensions[width*u]
            shortened = comb(blocks,u)*rank_total(dim,g,h) if dim >= h else 0
            caps.append(min(row[width*u], shortened, rank_total(k,g,h)))
        for u in range(blocks-1,-1,-1):
            caps[u] = min(caps[u],caps[u+1])
        result.append(caps)
    return result


def rank_one_zero_test(totals,zero):
    # Rank-one tuples have rows in {0,c}. Every nonzero column has the
    # same weight a (the number of copies of c), independently of c's bits.
    # Neither nonempty zero-feedback orbit has this form.
    for a in range(1,5):
        for active_columns in range(1,5):
            shape = (0,)*(4-active_columns)+(a,)*active_columns
            assert shape in totals and zero[shape] == 0


def self_test():
    # Exhaust every triple of words of a small even code and every fixed
    # pair-block support. Use actual shortened dimensions, then test the
    # block CDF and its decreasing-weight functional against exact counts.
    rows = [0b11110000,0b11001100,0b10101010]
    n,g,width = 8,3,2
    words = span(rows)
    dimensions = []
    for size in range(n+1):
        maximum = 1
        for support in combinations(range(n),size):
            mask = sum(1 << j for j in support)
            maximum = max(maximum,sum(w & ~mask == 0 for w in words))
        assert maximum & (maximum-1) == 0
        dimensions.append(maximum.bit_length()-1)
    bits = [[0]*(n+1) for _ in range(g)]
    blocks = [[0]*(n//width+1) for _ in range(g)]
    for tup in product(words,repeat=g):
        h = gf2_rank(tup)
        if not h:
            continue
        union = tup[0] | tup[1] | tup[2]
        bits[h-1][union.bit_count()] += 1
        blocks[h-1][sum(bool((union >> j)&3) for j in range(0,n,2))] += 1
    bounds = block_caps([list(accumulate(row)) for row in bits],dimensions,
                        n,width,g,len(rows))
    for row,cap in zip(blocks,bounds):
        assert all(a <= b for a,b in zip(accumulate(row),cap))
        for length in range(5):
            for q in (F(0),F(1,9216),F(1)):
                p = prefix_probabilities(4,length,q)
                assert sum(a*b for a,b in zip(row,p)) <= weighted_cdf_bound(cap,p)
    print('Exhaustive small-code block-support/CDF checks passed',flush=True)


def main():
    self_test()
    totals,zero = census()
    rank_one_zero_test(totals,zero)
    q = max(F(zero[s],v) for s,v in totals.items() if any(s))
    assert q == F(1,9216)
    print('Rank-one zero feedback: impossible; other ranks: at most',q,flush=True)
    spectrum = authenticated_caps()
    dimensions = dimension_caps()
    bits = support_caps(spectrum,g=4,dimensions=dimensions)
    caps = block_caps(bits,dimensions)
    for h,row in enumerate(caps,1):
        print('rank',h,'block support lower bound',next(u for u,v in enumerate(row) if v),
              'log2 CDF',[(u,log2(row[u]) if row[u] else None) for u in (10,15,20,24,32,48,64)],flush=True)
    groups = (1 << 20)//128//4
    logarithm = lambda v: log2(v.numerator)-log2(v.denominator) if v else float('-inf')
    for length in (40,44,48,50,52,53,54,55,56,57,58,60,63):
        bounds = [groups*weighted_cdf_bound(row,prefix_probabilities(64,length,F(0) if h==1 else q))
                  for h,row in enumerate(caps,1)]
        total = sum(bounds,F(0))
        print('prefix',length,'all ranks log2',round(logarithm(total),6),
              'rank log2',[round(logarithm(v),6) for v in bounds],
              'below 2^-40',total < F(1,1 << 40),flush=True)
        if length == 57:
            assert 0 < total < F(1,1 << 42)
            print('EXACT CHECK: one-active-group zero-prefix union < 2^-42 at 57 macroregions',flush=True)
    print('No full-distance claim: these events only bound a completely zero state prefix.')


if __name__ == '__main__':
    main()
