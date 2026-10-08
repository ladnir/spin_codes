from fractions import Fraction as F
from itertools import combinations, product
import unittest

from live_kernel_rank import live_kernel_bound


def mul(a, b):
    v = 0
    while b:
        if b & 1:
            v ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return v


M = tuple(tuple(mul(a, b) for b in range(4)) for a in range(4))
INV = (0, 1, 3, 2)


def rank(rows):
    a = [list(row) for row in rows]
    if not a:
        return 0
    k = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(k, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        a[k] = [M[INV[a[k][column]]][x] for x in a[k]]
        for i in range(k+1, len(a)):
            scale = a[i][column]
            a[i] = [x ^ M[scale][y] for x, y in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def subspaces(m, d):
    for pivots in combinations(range(m), d):
        slots = [(i, j) for i, p in enumerate(pivots)
                 for j in range(p+1, m) if j not in pivots]
        for values in product(range(4), repeat=len(slots)):
            rows = [[0]*m for _ in range(d)]
            for i, p in enumerate(pivots):
                rows[i][p] = 1
            for (i, j), v in zip(slots, values):
                rows[i][j] = v
            yield tuple(tuple(row) for row in rows)


def dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= M[x][y]
    return result


def vectors(basis):
    return [tuple(dot(coeff, column) for column in zip(*basis))
            for coeff in product(range(4), repeat=len(basis))]


class LiveKernelRankTests(unittest.TestCase):
    def test_exhaustive_gf4_dimension_three(self):
        # All source subspaces and all equation row spaces in F4^3.
        # Different bases of the same row space impose the same kernel.
        checked = 0
        for d in (1, 2, 3):
            for basis in subspaces(3, d):
                live = [v for v in vectors(basis) if all(v)]
                for r in (1, 2, 3):
                    for b in subspaces(3, r):
                        histogram = [0]*(d+1)
                        for diag in product(range(1, 4), repeat=3):
                            bd = [tuple(M[x][y] for x, y in zip(row, diag)) for row in b]
                            if not any(all(dot(row, v) == 0 for row in bd) for v in live):
                                continue
                            eq = [[dot(row, v) for v in basis] for row in bd]
                            histogram[rank(eq)] += 1
                        for tau in (F(1), F(2), F(3)):
                            actual = sum(F(count, 27)*tau**(-gain)
                                         for gain, count in enumerate(histogram))
                            sharp = live_kernel_bound(d, 3, r, tau, 4)
                            self.assertLessEqual(actual, sharp)
                            self.assertLessEqual(sharp, live_kernel_bound(d, 3, r, tau, 4, False))
                            checked += 1
        self.assertEqual(checked, 5547)

    def test_monotone_dimension_and_coordinate_bounds(self):
        for tau in (F(1), F(6), F(255)):
            for r in (1, 4, 8):
                for m in range(max(r, 1), 17):
                    previous = F(0)
                    for d in range(1, m+1):
                        value = live_kernel_bound(d, m, r, tau)
                        self.assertGreaterEqual(value, previous)
                        previous = value
                        if m < 16:
                            self.assertLessEqual(value, live_kernel_bound(d, m+1, r, tau))

    def test_sparse_return_gain_and_full_rank(self):
        tau = F(6)
        bound = live_kernel_bound(2, 9, 8, tau)
        self.assertEqual(bound, F(256**2-1, 255**9)/tau)
        self.assertLess(bound, tau**-8/F(10**10))
        self.assertEqual(live_kernel_bound(3, 3, 3, tau), 0)


if __name__ == '__main__':
    unittest.main()
