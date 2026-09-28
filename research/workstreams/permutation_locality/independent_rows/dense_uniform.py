"""Outward dense-component bound for independently shuffled BCH rows.

This is a partial first-moment decomposition, not an occupancy certificate.
The number j counts uniform components among 8192 binary outer rows, not
active four-row groups. The components with smaller j remain unbounded.
"""
import argparse
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys

from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bch_joint_support import authenticated_caps
from group_rank_one_verify import up


def decomposition(caps, exponent=126):
    """Return c and residual shell masses for mu <= c U_even + D.

    A row's averaged counting density at a word of weight w is
    A[w]/binom(n,w). The supplied caps dominate A coefficientwise.
    U_even is uniform on the even-weight subspace, of dimension n-1.
    The density of c U_even is 2**(-exponent) at each even word.
    """
    n = len(caps)-1
    if n < 1 or exponent < 0 or any(a < 0 for a in caps):
        raise ValueError('positive length and nonnegative caps/exponent required')
    if any(caps[w] for w in range(1, n+1, 2)):
        raise ValueError('the outer code must be even')
    density = Q(1, 1 << exponent)
    residual = [max(Q(0), Q(a)-density*comb(n, w)) if w % 2 == 0 else Q(0)
                for w, a in enumerate(caps)]
    coefficient = Q(2)**(n-1-exponent)
    # Check the pointwise majorant on every weight shell exactly.
    assert all(Q(a) <= (coefficient*Q(comb(n, w), 1 << (n-1))
                        if w % 2 == 0 else Q(0)) + residual[w]
               for w, a in enumerate(caps))
    return coefficient, residual


def rational(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def component_bound(j, coefficient, residual_mass, rows=8192,
                    row_dimension=255, threshold=209715):
    """Whole labeled-j contribution; directed Arb upper.

    The arbitrary residual coordinates give an affine shift of an
    injective image of j independent even-row subspaces. Its dimension
    is row_dimension*j. An information set gives the tilted moment below.
    """
    if not 0 <= j <= rows or threshold < 0:
        raise ValueError('invalid component count or weight threshold')
    dimension = row_dimension*j
    if dimension <= 2*threshold:
        z = Q(1)
    else:
        z = Q(threshold, dimension-threshold)
    if threshold == 0:
        # z=0 is the limiting Chernoff witness; avoid 0**0 in Arb.
        moment = arb(2)**(-dimension)
    else:
        za = rational(z)
        moment = za**(-threshold)*((1+za)/2)**dimension
    return up(arb(comb(rows, j))*rational(coefficient)**j
              *rational(residual_mass)**(rows-j)*moment)


def certified_tail(coefficient, residual_mass, budget=52, rows=8192,
                   row_dimension=255, threshold=209715):
    """Largest downward-contiguous tail satisfying the aggregate budget."""
    limit = arb(2)**(-budget)
    total = arb(0)
    first = None
    for j in range(rows, -1, -1):
        term = component_bound(j, coefficient, residual_mass, rows,
                               row_dimension, threshold)
        candidate = up(total+term)
        if not candidate < limit:
            return first, total, j, candidate
        first, total = j, candidate
    return first, total, None, None


def span(basis):
    result = [0]
    for vector in basis:
        result += [x ^ vector for x in result]
    assert len(set(result)) == len(result)
    return result


def self_test():
    """Exact small-code checks, including nonzero affine shifts."""
    from itertools import product
    checks = 0
    for n, basis in ((3, [3, 5]), (4, [3, 5, 9]),
                     (5, [7, 19]), (6, [3, 12, 48])):
        code = span(basis)
        caps = [sum(x.bit_count() == w for x in code) for w in range(n+1)]
        if all(x.bit_count() % 2 == 0 for x in code):
            for exponent in (0, 1, n-1, n+1):
                coefficient, residual = decomposition(caps, exponent)
                density = Q(coefficient, 1 << (n-1))
                for x in range(1 << n):
                    w = x.bit_count()
                    actual = Q(caps[w], comb(n, w))
                    dominant = (density if w % 2 == 0 else 0) + residual[w]/comb(n, w)
                    assert actual <= dominant
                    checks += 1
        for z in (Q(0), Q(1, 8), Q(1, 2), Q(7, 8), Q(1)):
            unshifted = sum(z**x.bit_count() for x in code)/len(code)
            information = ((1+z)/2)**len(basis)
            for offset in range(1 << n):
                shifted = sum(z**(x ^ offset).bit_count() for x in code)/len(code)
                assert shifted <= unshifted <= information
                checks += 1
    # Exact positive tensor expansion, without assuming row-label events.
    for rows in (1, 2, 3, 4):
        for coefficient, residual_mass in product((Q(1, 2), Q(3), Q(7, 3)), repeat=2):
            assert sum(Q(comb(rows, j))*coefficient**j*residual_mass**(rows-j)
                       for j in range(rows+1)) == (coefficient+residual_mass)**rows
            checks += 1
    print('Dense-component exact domination, affine, and tensor checks:', checks, 'passed', flush=True)


def feedback_rank_check():
    """Exact precursor for a future uniform-input transfer, not a tilt bound."""
    from feedback_character_census import census
    histograms = census()
    zero = (32, 0, 0, 0, 0)
    assert histograms[zero] == 1
    nonzero = {h: count for h, count in histograms.items() if h != zero}
    maximum = max(h[0] for h in nonzero)
    assert maximum == 13
    print('VERIFIED: every subset of 14 feedback windows has full rank 19.', flush=True)
    for windows in (5, 6, 7, 8, 10, 12, 14):
        annihilators = Q(sum(count*comb(h[0], windows) for h, count in nonzero.items()
                            if h[0] >= windows), comb(32, windows))
        zero_probability = (1+annihilators)/(1 << 19)
        print('Uniform-window count', windows, 'rank-deficiency upper', min(Q(1), annihilators),
              'exact average zero-feedback probability', zero_probability, flush=True)
    print('These are untilted feedback facts; they do not assert uniformity after output tilting.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, default=126)
    parser.add_argument('--budget', type=int, default=52)
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--self-test-only', action='store_true')
    parser.add_argument('--feedback-rank', action='store_true')
    args = parser.parse_args()
    if args.precision < 128 or args.budget < 0 or args.exponent < 0:
        parser.error('precision >=128 and nonnegative budget/exponent required')
    self_test()
    if args.feedback_rank:
        feedback_rank_check()
    if args.self_test_only:
        return
    caps = authenticated_caps()
    # Imported legacy verifiers can change the process-wide Arb precision.
    ctx.prec = args.precision
    coefficient, residual = decomposition(caps, args.exponent)
    residual_mass = sum(residual, Q(0))
    print('Uniform even-row coefficient:', coefficient, flush=True)
    print('Residual mass log2:', rational(residual_mass).log()/arb(2).log(), flush=True)
    print('Residual support weights:', [w for w, a in enumerate(residual) if a], flush=True)
    first, total, excluded, excluded_total = certified_tail(coefficient, residual_mass, args.budget)
    if first is None:
        print('No nonempty component tail meets the requested budget.', flush=True)
        return
    assert total < arb(2)**(-args.budget)
    print('VERIFIED uniform-component count j in', first, '.. 8192, precision', args.precision,
          'upper', total, 'margin', -total.log()/arb(2).log(), flush=True)
    if excluded is not None:
        print('Next excluded count:', excluded, 'aggregate log2 upper',
              excluded_total.log()/arb(2).log(), flush=True)
    print('This is a partial positive first-moment decomposition, NOT active-group occupancy coverage.', flush=True)
    print('All lower-j components remain; NOT a whole-code distance certificate.', flush=True)


if __name__ == '__main__':
    main()
