"""Independent bounded checks of the uniform permitted-pair upper bound."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from flint import arb, ctx
import hd_component_pair as pairs
from hd_component_rank import gap_moment_exact, exact_kernels, arb_fraction


class PairTests(unittest.TestCase):
    def setUp(self):
        self.previous = ctx.prec
        ctx.prec = 128

    def tearDown(self):
        ctx.prec = self.previous

    def assert_upper(self, bound, exact):
        self.assertTrue(bound+arb(2)**-95 > arb_fraction(exact))

    def test_vandermonde_first_return_sum(self):
        for rows, base in product(range(2, 11), (F(1, 2), F(9, 10))):
            actual = sum((h-1)*base**(rows-2)*gap_moment_exact(rows, h, 1, base)
                         for h in range(2, rows+1))
            self.assertEqual(actual, pairs.first_return_sum(rows, base))
            for gap in range(rows-1):
                lhs = sum(h*(h-1)*comb(rows-h, gap) for h in range(2, rows-gap+1))
                self.assertEqual(lhs, 2*comb(rows+1, gap+3))

    def test_two_point_range_formula(self):
        for rows, base in product(range(2, 11), (F(1, 2), F(9, 10))):
            actual = sum(base**(right-left-1) for left, right in combinations(range(rows), 2))/comb(rows, 2)
            formula = 2*((rows-1)-rows*base+base**rows)/(rows*(rows-1)*(1-base)**2)
            self.assertEqual(actual, formula)

    def test_uniform_entries_for_all_small_occupancies(self):
        for rows, q, mark in product((2, 5, 8), (4, 16), (F(1, 2), F(9, 10))):
            kernels = exact_kernels(rows, q, mark, collapsed=True)
            for minimum in (1, 2):
                bound = pairs.uniform_positive(rows, q, mark, minimum)
                for h in range(minimum, rows+1):
                    for i, j in product(range(2), repeat=2):
                        self.assert_upper(bound[i][j], kernels[h][i][j])

    def test_pair_matrix_all_small_cut_positions(self):
        for rows in range(2, 9):
            q, mark = 16, F(9, 10)
            kernels = exact_kernels(rows, q, mark, collapsed=True)
            bound = pairs.pair_matrix(rows, q, mark)
            for left, right in combinations(range(rows), 2):
                for incoming, outgoing in product(range(2), repeat=2):
                    for i, j in product(range(2), repeat=2):
                        actual = sum((q-1)**(middle+outgoing)*kernels[
                            left*incoming+(right-left)*middle+(rows-right)*outgoing][i][j]
                                     for middle in range(2))
                        self.assert_upper(bound[2*incoming+i][2*outgoing+j], actual)

    def test_directed_k20_pair_gate(self):
        result = pairs.receipt()
        self.assertTrue(result['all_pairs_below_2_to_minus_40'])
        self.assertTrue(result['any_fixed_regionwise_row_order'])
        self.assertFalse(result['whole_distance_certificate'])


if __name__ == '__main__':
    unittest.main()
