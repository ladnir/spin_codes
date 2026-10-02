"""Exact support comparison for two independently shuffled BCH row pairs.

The CDF bounds expected message counts over the two pair shuffles. It is
not the deterministic spectrum of a realized four-row outer. See
PAIRWISE_SHUFFLE.md for the coupling and conditional GF16 argument.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb, lcm, log2
from pathlib import Path

from shared_support import shared_counts
from bch_joint_support import rank_total


def pair_cdf(ranks, k=128):
    """Rescale valid rank-h four-tuple CDF caps into two-tuple caps.

    Each h-dimensional subspace has rank_total(h,g,h) ordered spanning
    g-tuples. Integer division first bounds the number of subspaces.
    The returned CDF includes the unique zero pair.
    """
    if (type(k) is not int or k < 2 or len(ranks) != 4
            or not ranks[0] or any(len(row) != len(ranks[0]) for row in ranks)):
        raise ValueError('four rank CDFs and dimension at least two required')
    result = [1] * len(ranks[0])
    for h in (1, 2):
        row = ranks[h-1]
        if (row[0] != 0 or row[-1] != rank_total(k, 4, h)
                or any(type(c) is not int or c < 0 for c in row)
                or any(a > b for a, b in zip(row, row[1:]))):
            raise ValueError('valid rank CDF with exact total required')
        source = rank_total(h, 4, h)
        target = rank_total(h, 2, h)
        result = [a + (b // source)*target for a, b in zip(result, row)]
    assert result[0] == 1 and result[-1] == 1 << (2*k)
    return result


def union_shells(left, right):
    """Exact union distribution of two independent exchangeable subsets.

    Inputs are nonnegative integer *comparison masses* by subset size.
    Outputs are rational masses, without rounding or a floating proposal.
    A binomial zeta transform computes subset-containment masses, which
    multiply under independent union. Binomial inversion recovers shells.
    This takes O(n^2) integer operations rather than enumerating triples.
    """
    if (not left or len(left) != len(right)
            or any(type(v) is not int or v < 0 for v in (*left, *right))):
        raise ValueError('equal nonempty nonnegative integer mass arrays required')
    n = len(left)-1
    choose = [comb(n, a) for a in range(n+1)]
    denominator = lcm(*choose)
    x = [v*(denominator//c) for v, c in zip(left, choose)]
    y = [v*(denominator//c) for v, c in zip(right, choose)]
    contained = []
    for u in range(n+1):
        contained.append(sum(comb(u, a)*x[a] for a in range(u+1))
                         * sum(comb(u, b)*y[b] for b in range(u+1)))
    result = []
    for u in range(n+1):
        exact = sum((-1)**(u-j)*comb(u, j)*contained[j] for j in range(u+1))
        assert exact >= 0
        result.append(Q(choose[u]*exact, denominator**2))
    assert sum(result) == sum(left)*sum(right)
    return result


def combine_pair_cdfs(left, right):
    """Upper expected nonzero four-row CDF, from two pair upper CDFs.

    Differences are masses of stochastically smaller comparison laws,
    NOT bounds on the actual pair shell counts. Monotone union coupling
    justifies using these masses. The unique all-zero tuple is removed.
    """
    for cdf in (left, right):
        if (not cdf or cdf[0] != 1 or any(type(v) is not int or v < 0 for v in cdf)
                or any(a > b for a, b in zip(cdf, cdf[1:]))):
            raise ValueError('integer pair CDF including the unique zero required')
    masses = lambda cdf: [cdf[0]] + [b-a for a, b in zip(cdf, cdf[1:])]
    shells = union_shells(masses(left), masses(right))
    cumulative = Q(-1)
    result = []
    for value in shells:
        cumulative += value
        result.append(cumulative)
    assert result[0] == 0 and result[-1] == left[-1]*right[-1]-1
    return result


def pairwise_counts(refined=False, coupled=False):
    shared, ranks = shared_counts(refined=refined, coupled=coupled)
    pairs = pair_cdf(ranks)
    exact = combine_pair_cdfs(pairs, pairs)
    # Existing CDF fold engines take integer caps. Ceiling only enlarges
    # each expected CDF; minimum support and the exact total are unchanged.
    counts = [-(-v.numerator//v.denominator) for v in exact]
    assert counts[0] == 0 and counts[-1] == (1 << 512)-1
    assert all(v == 0 for v in counts[:38])
    return counts, dict(pair_cdf=pairs, expected_cdf=exact, shared_cdf=shared)


def tighten_cdf(cdf,shell_caps):
    """Intersect a CDF cap with prefix sums of genuine shell caps."""
    if (len(cdf)!=len(shell_caps) or not cdf or cdf[0]!=1 or shell_caps[0]!=1
            or any(type(v) is not int or v<0 for v in (*cdf,*shell_caps))
            or any(a>b for a,b in zip(cdf,cdf[1:]))):
        raise ValueError('matching nonnegative pair bounds including zero required')
    result=[];total=0
    for cap,shell in zip(cdf,shell_caps):
        total+=shell;result.append(min(cap,total))
    if result[-1]!=cdf[-1]:
        raise ArithmeticError('shell bounds contradict the exact total')
    return result


def shell_refined_counts():
    """Regenerate moment/transform caps before independent pair union."""
    from pairwise_macwilliams import actual
    model=actual()
    shells=[1]+[model.transform_cap(u) for u in range(1,257)]
    pairs=tighten_cdf(model.cdfs[0],shells)
    exact=combine_pair_cdfs(pairs,pairs)
    counts=[-(-v.numerator//v.denominator) for v in exact]
    for u in (64,80,96,112,128,144,160,192):
        print('PAIR SHELL-CDF support',u,'gain bits',log2(model.cdfs[0][u])-log2(pairs[u]),flush=True)
    return counts,dict(pair_cdf=pairs,original_pair_cdf=model.cdfs[0],
                      pair_shell_caps=shells,expected_cdf=exact)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refined-counts', action='store_true')
    parser.add_argument('--coupled-counts', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts, data = pairwise_counts(args.refined_counts, args.coupled_counts)
    for u in (38, 57, 80, 96, 104, 108, 112, 116, 120, 128, 144, 192, 256):
        old, new = data['shared_cdf'][u], counts[u]
        print('support', u, 'shared log2', log2(old), 'pairwise log2', log2(new),
              'reduction bits', log2(old)-log2(new), flush=True)
    if args.output:
        record = dict(schema='pairwise-gf16-cdf-1', refined_counts=args.refined_counts,
                      coupled_counts=args.coupled_counts, counts=counts,
                      pair_cdf=data['pair_cdf'],
                      expected_cdf=[[v.numerator, v.denominator] for v in data['expected_cdf']],
                      note='Exact comparison CDF arithmetic. No distance certificate.')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    main()
