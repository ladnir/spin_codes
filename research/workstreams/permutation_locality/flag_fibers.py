"""Exact bounds for all-one-column flags using averaged extension fibers.

Counts only: no change to SPIN or its probability model.
"""
from itertools import combinations, product
from math import comb, log2, prod
from fractions import Fraction as Q
from bch_joint_support import rank_total
from shortening_moments import containment_moment
from joint_support import span, gf2_rank


def fiber_cdfs(caps,dimensions,g=4):
    """Bound sum_H 2^(dim C_support(H)-r), rank(H)=r, support(H)<=u."""
    n = len(dimensions)-1
    spaces = [[x//rank_total(r,g,r) for x in row] for r,row in enumerate(caps,1)]
    moments = [[containment_moment(row,t) for t in range(n+1)] for row in spaces]
    result = []
    for r in range(1,g):
        row = []
        for u in range(n+1):
            best = spaces[r-1][u]*(1 << (dimensions[u]-r)) if dimensions[u]>=r else 0
            for t in range(u,n+1):
                # 2^(d-r) G(d,r) = G(d,r)+(2^(r+1)-1)G(d,r+1).
                value = (moments[r-1][t]+((1<<(r+1))-1)*moments[r][t])//comb(n-u,t-u)
                best = min(best,value)
            row.append(best)
        for u in range(n-1,-1,-1):row[u] = min(row[u],row[u+1])
        result.append(row)
    return result


def prefix_caps(spectrum,caps,dimensions,g=4):
    """Count tuples with support<=u and exactly j all-one columns (j>0)."""
    n = len(dimensions)-1
    fibers = fiber_cdfs(caps,dimensions,g)
    result = [[0]*(n+1) for _ in range(n+1)]
    for j in range(1,n+1):
        for u in range(j,n+1):
            total = spectrum[j]  # Rank one: (w,w,w,w).
            for r,row in enumerate(fibers,1):
                # binom(n-v,j) decreases in v; CDF summation is safe.
                value = sum((row[v]-(row[v-1] if v else 0))*comb(n-v,j) for v in range(u-j+1))
                tuples = (1<<r)*prod((1<<(g-1))-(1<<i) for i in range(r))
                total += min(caps[r][u],tuples*value)
            result[u][j] = total
    return result


def self_test():
    checks = 0
    for basis,n in (([1,2,4],4),([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8)):
        code = span(basis);g=4
        spectrum = [sum(w.bit_count()==v for w in code) for v in range(n+1)]
        actual = [[0]*(n+1) for _ in range(n+1)]
        shells = [[0]*(n+1) for _ in range(g)]
        for xs in product(code,repeat=g):
            u = (xs[0]|xs[1]|xs[2]|xs[3]).bit_count()
            j = (xs[0]&xs[1]&xs[2]&xs[3]).bit_count()
            r = gf2_rank(xs)
            if r:shells[r-1][u]+=1
            actual[u][j]+=1
        caps = [[sum(row[:u+1]) for u in range(n+1)] for row in shells]
        dimensions = [0]*(n+1)
        exact_dims = {}
        for mask in range(1<<n):
            d = sum(w&~mask == 0 for w in code).bit_length()-1
            exact_dims[mask] = d
            dimensions[mask.bit_count()] = max(dimensions[mask.bit_count()],d)
        fibers = fiber_cdfs(caps,dimensions,g)
        for r,row in enumerate(fibers,1):
            subspaces = {tuple(sorted(span(xs))) for xs in combinations(code[1:],r) if gf2_rank(xs)==r}
            exact = [0]*(n+1)
            for subspace in subspaces:
                support=0
                for x in subspace:support|=x
                exact[support.bit_count()] += 1 << (exact_dims[support]-r)
            for u in range(n+1):
                assert sum(exact[:u+1])<=row[u];checks+=1
        for upper in (prefix_caps(spectrum,caps,dimensions),
                      prefix_caps([2*x for x in spectrum],[[2*x for x in row] for row in caps],dimensions)):
            for u in range(n+1):
                for j in range(1,n+1):
                    assert sum(actual[v][j] for v in range(u+1))<=upper[u][j];checks+=1
    print('Averaged flag fibers:',checks,'exact small-code count inequalities passed',flush=True)


def main():
    from bch_joint_support import authenticated_caps,support_caps
    from shortened_bound import dimension_caps
    from shortening_polynomial import improve_dimensions
    from dual_shortening import improve_dimensions as dual_dimensions
    from basis_lattice import improve_caps
    from shortening_moments import improve
    from dual_moments import refine_bch
    from occupancy_allones import weighted_cdf
    self_test()
    spectrum=authenticated_caps()
    dimensions=dual_dimensions(improve_dimensions(dimension_caps()))
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    caps,_=improve(caps,dimensions);caps=refine_bch(caps,dimensions)
    extra=prefix_caps(spectrum,caps,dimensions)
    for penalty in ('1/2','3/4','9/10'):
        old=weighted_cdf(spectrum,caps,dimensions,penalty,prefix_flags=True)
        new=weighted_cdf(spectrum,caps,dimensions,penalty,prefix_flags=True,extra_prefix=extra)
        assert all(a<=b for a,b in zip(new,old))
        print('Penalty',penalty,'log2 weighted counts before/after',
              [(u,log2(old[u]),log2(new[u])) for u in (80,96,112,128,144,160,176,192)],flush=True)
    print('Outer-count refinement only; no full occupancy certificate.',flush=True)


if __name__=='__main__':main()
