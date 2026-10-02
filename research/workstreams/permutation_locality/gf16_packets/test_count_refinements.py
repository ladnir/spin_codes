import unittest
from fractions import Fraction as Q
from itertools import product

from count_refinements import basis_laplace_cap
from joint_support import span, gf2_rank


class BasisLaplaceTests(unittest.TestCase):
    def test_exhaustive_small_codes(self):
        for basis, n in (([1, 2, 4], 3), ([3, 5], 3), ([15, 51, 85], 7)):
            code = span(basis)
            spectrum = [sum(w.bit_count() == u for w in code) for u in range(n+1)]
            shells = [[0]*(n+1) for _ in range(4)]
            for words in product(code, repeat=3):
                shells[gf2_rank(words)][(words[0] | words[1] | words[2]).bit_count()] += 1
            for h in range(1, 4):
                for u in range(n+1):
                    for tilt in (Q(1, 8), Q(1, 3), Q(7, 8), Q(1)):
                        self.assertLessEqual(sum(shells[h][:u+1]), basis_laplace_cap(spectrum, 3, h, u, tilt))

    def test_spectrum_upper_bounds(self):
        exact = [1, 0, 3, 0]
        inflated = [1, 0, 6, 0]
        for h in range(1, 4):
            for u in range(4):
                self.assertGreaterEqual(basis_laplace_cap(inflated, 3, h, u, Q(1, 2)),
                                        basis_laplace_cap(exact, 3, h, u, Q(1, 2)))


if __name__ == '__main__':
    unittest.main()
