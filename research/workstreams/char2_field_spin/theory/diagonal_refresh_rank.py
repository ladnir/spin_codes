"""Rank-moment bound for a cheap diagonal state refresh.

Requires a superregular emission matrix. This is NOT the original monomial
emission map with a diagonal inserted. No native encoder or certificate is
provided. The exact tests enumerate every GF8 subspace in dimension three
and every nonzero diagonal, and compare the analytic potential bound.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json


def rank_potential_bound(rank, roots, tau, q=256, state=4):
    """Uniform E[tau**(-rank(A_Z D S)) ; remaining state rank>0].

    S has rank `rank`; all roots-row/column minors of A_Z are nonzero.
    Uses only a first moment of the surviving subspace size, not a claim
    that diagonal refresh samples a uniform Grassmannian subspace.
    """
    if not (1 <= rank <= state and 0 <= roots <= state and 1 <= tau <= q-1):
        raise ValueError('invalid rank, roots, or tilt')
    if isinstance(tau, int):
        tau = Fraction(tau)
    if roots == 0:
        return Fraction(1)
    maximum = min(rank, state-roots)
    if maximum == 0:
        return Fraction(0)
    moment = min(Fraction(q**rank-1, (q-1)**roots), Fraction(q**maximum-1))
    left_x, left_y = Fraction(0), Fraction(0)
    for h in range(1, maximum+1):
        right_x = Fraction(q**h-1)
        right_y = tau**(h-rank)
        if moment <= right_x:
            return left_y+(right_y-left_y)*(moment-left_x)/(right_x-left_x)
        left_x, left_y = right_x, right_y
    raise AssertionError('capped moment must lie on the polygon')


def root_factors(tau, q=256):
    import numpy as np
    return np.asarray([float(max(rank_potential_bound(r, z, tau, q) for r in range(1, 5)))
                       for z in range(4)])


def gf8_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return result


MUL = tuple(tuple(gf8_mul(a, b) for b in range(8)) for a in range(8))
INV = (0,)+tuple(next(b for b in range(1, 8) if MUL[a][b] == 1) for a in range(1, 8))


def matrix_rank(matrix):
    work = [list(row) for row in matrix]
    if not work:
        return 0
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = INV[work[rank][column]]
        work[rank] = [MUL[inverse][x] for x in work[rank]]
        for row in range(rank+1, len(work)):
            factor = work[row][column]
            work[row] = [x ^ MUL[factor][y] for x, y in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def subspaces(dimension, rank):
    """Canonical RREF bases, each GF8 subspace exactly once."""
    for pivots in combinations(range(dimension), rank):
        slots = [(i, j) for i, pivot in enumerate(pivots)
                 for j in range(pivot+1, dimension) if j not in pivots]
        for entries in product(range(8), repeat=len(slots)):
            rows = [[0]*dimension for _ in range(rank)]
            for i, pivot in enumerate(pivots):
                rows[i][pivot] = 1
            for (i, j), value in zip(slots, entries):
                rows[i][j] = value
            yield rows


def exact_checks():
    from regional_refresh_screen import gaussian
    # Three input columns, two rows: disjoint Cauchy nodes {0,1,2}/{3,4}.
    emission = [[INV[x ^ y] for y in range(3)] for x in (3, 4)]
    for z in (1, 2):
        for columns in combinations(range(3), z):
            assert matrix_rank([[row[j] for j in columns] for row in emission[:z]]) == z
    cases, diagonals = 0, 0
    for rank in range(1, 4):
        spaces = list(subspaces(3, rank))
        assert len(spaces) == gaussian(3, rank, 8)
        for basis in spaces:
            histogram = [[0]*(rank+1) for _ in range(2)]
            for diagonal in product(range(1, 8), repeat=3):
                coefficients = [[MUL[row[j]][diagonal[j]] for j in range(3)] for row in emission]
                equations = [[MUL[c[0]][u[0]] ^ MUL[c[1]][u[1]] ^ MUL[c[2]][u[2]]
                              for u in basis] for c in coefficients]
                for z in (1, 2):
                    h = rank-matrix_rank(equations[:z])
                    histogram[z-1][h] += 1
                diagonals += 1
            for z in (1, 2):
                frequencies = [Fraction(count, 7**3) for count in histogram[z-1]]
                moment = sum(p*(8**h-1) for h, p in enumerate(frequencies))
                assert moment <= Fraction(8**rank-1, 7**z)
                for tau in (Fraction(1), Fraction(2), Fraction(7)):
                    actual = sum(frequencies[h]*tau**(h-rank) for h in range(1, rank+1))
                    assert actual <= rank_potential_bound(rank, z, tau, 8, 3)
                    cases += 1
    # Regression: the formula needs a superregular matrix, not just row rank.
    # A=[1,0], U=span(e2) has H=1 always despite roots=1.
    assert rank_potential_bound(1, 1, Fraction(255), 256, 2) < 1
    return {'exact_potential_cases': cases, 'exhausted_diagonal_subspace_pairs': diagonals}


def compare():
    from regional_refresh_screen import root_factors as gl_factors
    return [dict(tau=tau, diagonal=root_factors(Fraction(tau)).tolist(),
                 uniform_gl=gl_factors(tau).tolist()) for tau in (2, 14, 255)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--compare', action='store_true')
    args = parser.parse_args()
    if args.check:
        print(json.dumps(exact_checks()), flush=True)
    if args.compare:
        for result in compare():
            print(json.dumps(result), flush=True)
