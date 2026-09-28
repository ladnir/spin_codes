"""Joint-support CDF bounds from moments of shortened-code dimensions.

Exact integer/rational proof component, not a full SPIN certificate.
Count r-dimensional subcodes inside each t-coordinate set, then use a
uniform shortening dimension bound to bound its h-dimensional subcodes.
"""
from fractions import Fraction
import argparse
from itertools import combinations, product
from math import comb, log2, prod

from bch_joint_support import rank_total
from joint_support import span, gf2_rank


def gaussian(k, h):
    if h < 0 or h > k:
        return 0
    numerator = prod((1 << k) - (1 << j) for j in range(h))
    denominator = prod((1 << h) - (1 << j) for j in range(h))
    value, remainder = divmod(numerator, denominator)
    assert not remainder
    return value


def containment_moment(cdf, t):
    """Upper sum over t-sets of contained subcodes from their support CDF."""
    n = len(cdf)-1
    weights = [comb(n-v, t-v) for v in range(t+1)]
    assert all(x >= y for x, y in zip(weights, weights[1:]))
    return sum((cdf[v]-(cdf[v-1] if v else 0))*weights[v] for v in range(t+1))


def improve(caps, dimensions, g=4):
    """Input/output count ordered g-tuples, separated by rank."""
    n = len(dimensions)-1
    spaces = [[x//rank_total(h, g, h) for x in row]
              for h, row in enumerate(caps, 1)]
    moments = [[containment_moment(row, t) for t in range(n+1)] for row in spaces]
    result = [row[:] for row in caps]
    witnesses = {}
    for h in range(2, len(caps)+1):
        multiplicity = rank_total(h, g, h)
        for u in range(n+1):
            for t in range(u, n+1):
                d = dimensions[t]
                if d < h:
                    candidate = 0
                    r = 0
                else:
                    candidates = [(Fraction(gaussian(d,h)*moments[r-1][t], gaussian(d,r)),r)
                                  for r in range(1,h)]
                    bound, r = min(candidates)
                    # Each subcode supported on <=u coordinates lies in at
                    # least binomial(n-u,t-u) different t-coordinate sets.
                    bound *= Fraction(multiplicity, comb(n-u,t-u))
                    candidate = bound.numerator//bound.denominator
                if candidate < result[h-1][u]:
                    result[h-1][u] = candidate
                    witnesses[h,u] = t,r
        for u in range(n-1,-1,-1):
            result[h-1][u] = min(result[h-1][u], result[h-1][u+1])
        # Later ranks may reuse the already improved lower-rank CDF.
        spaces[h-1] = [x//multiplicity for x in result[h-1]]
        moments[h-1] = [containment_moment(spaces[h-1],t) for t in range(n+1)]
    return result, witnesses


def self_test():
    checks = 0
    for k in range(1,12):
        for h in range(2,min(k,4)+1):
            for r in range(1,h):
                ratios = [Fraction(gaussian(d,h),gaussian(d,r)) for d in range(h,k+1)]
                assert all(a<=b for a,b in zip(ratios,ratios[1:]))
    for basis,n in (([1,2,4],4),([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8)):
        words = span(basis); g = 4
        shells = [[0]*(n+1) for _ in range(4)]
        for xs in product(words,repeat=g):
            rank = gf2_rank(xs)
            if rank:
                shells[rank-1][(xs[0]|xs[1]|xs[2]|xs[3]).bit_count()] += 1
        actual = [[sum(row[:u+1]) for u in range(n+1)] for row in shells]
        dimensions = []
        for t in range(n+1):
            maximum = 0
            total_by_rank = [0]*4
            for coordinates in combinations(range(n),t):
                mask = sum(1<<j for j in coordinates)
                size = sum(w & ~mask == 0 for w in words)
                k = size.bit_length()-1
                maximum = max(maximum,k)
                for h in range(1,5):
                    total_by_rank[h-1] += gaussian(k,h)
            dimensions.append(maximum)
            for h in range(1,5):
                cdf = [v//rank_total(h,g,h) for v in actual[h-1]]
                assert containment_moment(cdf,t) == total_by_rank[h-1]
                checks += 1
        # Exact, inflated, and highly loose but valid input CDFs.
        for caps in (actual, [[2*x for x in row] for row in actual],
                     [actual[0]]+[[row[-1]]*(n+1) for row in actual[1:]]):
            improved,_ = improve(caps,dimensions,g)
            assert all(a<=b<=c for ar,br,cr in zip(actual,improved,caps)
                       for a,b,c in zip(ar,br,cr))
    print('Shortening moments: exact containment identities and CDF checks passed;',checks,'identities',flush=True)


def main():
    from bch_joint_support import authenticated_caps, support_caps
    from shortened_bound import dimension_caps
    from basis_lattice import improve_caps
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--positive-polynomials',action='store_true')
    parser.add_argument('--dual-shortening',action='store_true')
    args=parser.parse_args()
    self_test()
    spectrum = authenticated_caps(); dimensions = dimension_caps()
    old = improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    starting=old
    if args.positive_polynomials:
        from shortening_polynomial import improve_dimensions
        dimensions=improve_dimensions(dimensions)
        starting=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    if args.dual_shortening:
        from dual_shortening import improve_dimensions, self_test as dual_test
        dual_test();dimensions=improve_dimensions(dimensions)
        starting=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    new, witnesses = improve(starting,dimensions)
    for u in (72,80,96,112,128,144,160,176,192,224,256):
        before = sum(row[u] for row in old); after = sum(row[u] for row in new)
        print('u',u,'old/new log2',log2(before),log2(after),
              'rank witnesses',[(h,witnesses[h,u]) for h in range(2,5) if (h,u) in witnesses],flush=True)
    changed = [(log2(a)-log2(b),h,u,witnesses.get((h,u)))
               for h,(ar,br) in enumerate(zip(old,new),1) for u,(a,b) in enumerate(zip(ar,br)) if 0<b<a]
    print('Best rank-specific savings:',sorted(changed,reverse=True)[:12],flush=True)
    print('No inner bound or full occupancy certificate was changed.',flush=True)


if __name__ == '__main__':
    main()
