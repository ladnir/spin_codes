"""Independent finite subspace enumeration for the fresh-GL rank law."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import unittest

from regional_grassmann_gate import (
    alive_root_contraction, check, cost, gaussian_binomial,
    root_rank_distribution, zero_phase_action,
)


def multiply(a, b, q):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & q:
            a ^= {2: 3, 4: 7}[q]
    return result


def rank(matrix, n, q):
    work = [list(row) for row in matrix]
    current = 0
    for column in range(n):
        pivot = next((i for i in range(current, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[current], work[pivot] = work[pivot], work[current]
        inverse = next(b for b in range(1, q) if multiply(work[current][column], b, q) == 1)
        work[current] = [multiply(inverse, value, q) for value in work[current]]
        for i in range(current+1, len(work)):
            coefficient = work[i][column]
            work[i] = [a ^ multiply(coefficient, b, q) for a, b in zip(work[i], work[current])]
        current += 1
    return current


def subspaces(n, r, q):
    """Enumerate each row space once via its unique reduced echelon basis."""
    for pivots in combinations(range(n), r):
        slots = [(i, j) for i, pivot in enumerate(pivots)
                 for j in range(pivot+1, n) if j not in pivots]
        for values in product(range(q), repeat=len(slots)):
            matrix = [[int(j == pivot) for j in range(n)] for pivot in pivots]
            for (i, j), value in zip(slots, values):
                matrix[i][j] = value
            yield tuple(map(tuple, matrix))


def dot(a, b, q):
    value = 0
    for x, y in zip(a, b):
        value ^= multiply(x, y, q)
    return value


class GrassmannTests(unittest.TestCase):
    def test_all_small_subspaces(self):
        for q in (2, 4):
            for r in range(5):
                spaces = tuple(subspaces(4, r, q))
                self.assertEqual(len(spaces), gaussian_binomial(4, r, q))
                for z in range(5):
                    # K has its first z coordinates zero.  The restricted
                    # projection has rank r-h, computed without the formula.
                    counts = Counter(r-rank([row[:z] for row in basis], z, q)
                                     for basis in spaces)
                    expected = root_rank_distribution(q, 4, r, z)
                    self.assertEqual(tuple(Fraction(counts[h], len(spaces)) for h in range(5)), expected)

    def test_exact_q256_endpoints_and_slack(self):
        largest, _, _ = check()
        self.assertLess(largest, Fraction(257, 256))
        self.assertGreater(largest, 1)
        for r in range(1, 5):
            self.assertEqual(alive_root_contraction(256, 4, r, 0, 255), 1)
            self.assertEqual(alive_root_contraction(256, 4, r, 7, 255), 0)

    def test_arbitrary_continuations(self):
        # In particular, these continuation vectors are not all monotone.
        for values in ([1, 17, 3, 2, 41], [31, 2, 101, 4, 1], [0]*5):
            for r, z, f0, d0 in product(range(1, 5), range(5), range(3), range(3)):
                exact = zero_phase_action(256, 4, r, z, f0, d0, 255, values)
                plain = zero_phase_action(256, 4, r, z, f0, d0, 255, values, fresh=False)
                self.assertLessEqual(exact, plain)
                self.assertEqual(zero_phase_action(256, 4, r, z, f0, d0, 255,
                                                  [7*x for x in values]), 7*exact)

    def test_one_output_excludes_injection_column(self):
        # Every GF4 subspace satisfying a nonempty set of cubic root rows.
        # Here F_j=I+eta*c_j*a_j; a_j*c_j=0 and F_j*c_j=c_j.
        q = 4
        curves = tuple((1, u, multiply(u, u, q), multiply(multiply(u, u, q), u, q))
                       for u in range(q))
        checked = 0
        for j in range(q):
            injection = curves[j]
            emission = tuple(reversed(injection))
            self.assertEqual(dot(emission, injection, q), 0)
            for z in range(q):
                if z == j:
                    continue
                root = tuple(reversed(curves[z]))
                self.assertNotEqual(dot(root, injection, q), 0)
                for r in range(1, 4):
                    for basis in subspaces(4, r, q):
                        if any(dot(root, row, q) for row in basis):
                            continue
                        for eta in range(1, q):
                            transformed = tuple(tuple(x ^ multiply(multiply(eta, dot(emission, row, q), q), c, q)
                                                      for x, c in zip(row, injection)) for row in basis)
                            self.assertEqual(rank(transformed + (injection,), 4, q), r+1)
                            # Any dependent scalar form changes basis columns
                            # by multiples of c_j without losing their rank.
                            for coefficients in product(range(q), repeat=r):
                                updated = tuple(tuple(x ^ multiply(coefficient, c, q)
                                                      for x, c in zip(row, injection))
                                                for row, coefficient in zip(transformed, coefficients))
                                self.assertEqual(rank(updated, 4, q), r)
                                checked += 1
        self.assertGreater(checked, 1000)

    def test_costs_and_guards(self):
        self.assertEqual(cost(64)[2], Fraction(189, 88))
        self.assertEqual(cost(16)[2], Fraction(213, 88))
        for arguments in ((6, 4, 1, 1), (256, 4, 5, 0), (256, 4, 1, -1)):
            with self.assertRaises(ValueError):
                root_rank_distribution(*arguments)
        with self.assertRaises(ValueError):
            zero_phase_action(256, 4, 0, 0, 0, 0, 255, [1]*5)
        with self.assertRaises(ValueError):
            zero_phase_action(256, 4, 1, 0, 0, 0, 256, [1]*5)


if __name__ == "__main__":
    unittest.main()
