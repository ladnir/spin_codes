"""Positive dual-complement identities for average shortening moments.

An exact outer-count experiment, not a full SPIN certificate. No writes.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, log2
from bch_joint_support import rank_total, support_caps, authenticated_caps
from basis_lattice import improve_caps
from joint_support import span, gf2_rank
from shortened_bound import dimension_caps, krawtchouk
from shortening_moments import gaussian, containment_moment, improve
from dual_shortening import complement_caps, verify_bch_premise


def split_coefficient(a,h,r):
    value = gaussian(a,h-r)
    return value*(1 << (r*(a-h+r))) if value else 0


def kernel_cap(n,k,degree,w):
    ks = krawtchouk(n,w)
    kernel = sum((Q(ks[j]**2,comb(n,j)) for j in range(degree+1)),Q(0))
    bound = Q(1<<k)/kernel
    return bound.numerator//bound.denominator


def dual_shell_caps():
    # C has distance >=38, hence D=C^perp is an OA of strength >=37.
    # The squared degree-18 reproducing kernel uses only moments through36.
    result = [0]*257
    result[0] = result[256] = 1
    for w in range(30,227,2):
        result[w] = min(1<<128,comb(256,w),kernel_cap(256,128,18,w))
    assert result == result[::-1]
    return result


def improve_from_dual(primal,dual,k,g=4):
    """Rank-specific g-tuple CDFs; both codes have length n, primal dim k."""
    n = len(primal[0])-1
    moments = [[comb(n,t) for t in range(n+1)]]
    for r,row in enumerate(dual,1):
        spaces = [x//rank_total(r,g,r) for x in row]
        moments.append([containment_moment(spaces,t) for t in range(n+1)])
    result = [row[:] for row in primal]; witnesses = {}
    for h,row in enumerate(result,1):
        multiplicity = rank_total(h,g,h)
        # a=t-(n-k)>=0 ensures the Gaussian-binomial expansion is positive.
        totals = {t: sum(split_coefficient(t-(n-k),h,r)*moments[r][n-t]
                         for r in range(h+1)) for t in range(n-k,n+1)}
        for u in range(n+1):
            for t in range(max(u,n-k),n+1):
                upper = multiplicity*totals[t]//comb(n-u,t-u)
                if upper < row[u]:
                    row[u] = upper; witnesses[h,u] = t
        for u in range(n-1,-1,-1):
            row[u] = min(row[u],row[u+1])
    return result,witnesses


def exact_tuple_cdfs(words,n,g):
    shells = [[0]*(n+1) for _ in range(g)]
    for xs in product(words,repeat=g):
        rank = gf2_rank(xs)
        if rank:
            union = 0
            for x in xs: union |= x
            shells[rank-1][union.bit_count()] += 1
    return [[sum(row[:u+1]) for u in range(n+1)] for row in shells]


def self_test():
    checks = 0
    for a in range(9):
        for b in range(9):
            for h in range(1,5):
                assert gaussian(a+b,h) == sum(split_coefficient(a,h,r)*gaussian(b,r) for r in range(h+1))
                checks += 1
    for basis,n in (([3,5],4),([1,2,4],4),([15,51,85],7)):
        code = span(basis)
        dual = [w for w in range(1<<n) if all((w&r).bit_count()%2 == 0 for r in basis)]
        minimum = min(w.bit_count() for w in code if w)
        degree = (minimum-1)//2
        for w in range(n+1):
            assert sum(x.bit_count()==w for x in dual) <= kernel_cap(n,n-len(basis),degree,w)
            checks += 1
        g = 3
        exact = exact_tuple_cdfs(code,n,g)
        dcdf = exact_tuple_cdfs(dual,n,g)
        loose = [[row[-1]]*(n+1) for row in exact]
        for dinput in (dcdf,[[2*x for x in row] for row in dcdf]):
            bounded,_ = improve_from_dual(loose,dinput,len(basis),g)
            assert all(x<=y<=z for xs,ys,zs in zip(exact,bounded,loose) for x,y,z in zip(xs,ys,zs))
        # Verify the averaged Gaussian identity before any CDF upper bounds.
        for t in range(n-len(basis),n+1):
            for h in range(1,g+1):
                lhs = containment_moment([x//rank_total(h,g,h) for x in exact[h-1]],t)
                rhs = 0
                for r in range(h+1):
                    moment = (comb(n,n-t) if r==0 else
                              containment_moment([x//rank_total(r,g,r) for x in dcdf[r-1]],n-t))
                    rhs += split_coefficient(t-n+len(basis),h,r)*moment
                assert lhs == rhs; checks += 1
    print('Positive dual-moment identities:',checks,'exact checks passed',flush=True)


def refine_bch(primal,dimensions,iterations=3,dual_spectrum=None):
    verify_bch_premise()
    dspectrum = dual_shell_caps() if dual_spectrum is None else dual_spectrum
    assert len(dspectrum)==257 and dspectrum==dspectrum[::-1]
    ddims = complement_caps(dimension_caps(256,128,30,last_lp=104),dimensions,128)
    dual = improve_caps(support_caps(dspectrum,g=4,dimensions=ddims),dspectrum)
    dual,_ = improve(dual,ddims)
    result = primal
    for _ in range(iterations):
        result,_ = improve_from_dual(result,dual,128)
        dual,_ = improve_from_dual(dual,result,128)
    print('Dual moments log2 CDF before/after:',
          [(u,log2(sum(row[u] for row in primal)),log2(sum(row[u] for row in result)))
           for u in (96,128,144,160,176)],flush=True)
    return result


def main():
    self_test();verify_bch_premise()
    spectrum = authenticated_caps(); dspectrum = dual_shell_caps()
    from shortening_polynomial import improve_dimensions
    dims = improve_dimensions(dimension_caps())
    ddims = dimension_caps(256,128,30,last_lp=104)
    dims = complement_caps(dims,ddims,128)
    ddims = complement_caps(ddims,dims,128)
    primal = improve_caps(support_caps(spectrum,g=4,dimensions=dims),spectrum)
    primal,_ = improve(primal,dims)
    dual = improve_caps(support_caps(dspectrum,g=4,dimensions=ddims),dspectrum)
    dual,_ = improve(dual,ddims)
    for iteration in range(3):
        revised,witness = improve_from_dual(primal,dual,128)
        for u in (80,96,112,128,144,160,176,192,224):
            print('Iteration',iteration,'support',u,'log2 before/after',
                  log2(sum(row[u] for row in primal)),log2(sum(row[u] for row in revised)),
                  'witnesses',[(h,witness[h,u]) for h in range(1,5) if (h,u) in witness],flush=True)
        primal = revised
        dual,_ = improve_from_dual(dual,primal,128)
    print('Exact outer bounds only; no inner probability or full-code certificate.',flush=True)


if __name__ == '__main__':main()
