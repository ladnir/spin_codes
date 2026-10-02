"""Exact checks for uniform outer domination and regional averaging."""
from fractions import Fraction as Q
from itertools import product
import unittest

from flint import fmpq, fmpq_mat
import packet_uniform_tail as screen


class UniformTailTests(unittest.TestCase):
    def test_both_outer_shell_envelopes(self):
        for n, k, packets in ((8, 4, 8), (16, 8, 4)):
            beta, counts = screen.uniform_envelope(n, k, packets)
            self.assertEqual(sum(counts), 2**128 - 1)
            self.assertGreater(beta, 2**128)

    def test_binomial_matrix_average_matches_labeled_patterns(self):
        regional = [fmpq_mat([[1, j], [j + 1, 2]]) for j in range(5)]
        for q in range(5):
            actual = screen.regional_uniform(regional, q, matrix=fmpq_mat,
                rational=lambda x: fmpq(x.numerator, x.denominator), rounding=lambda x: x)
            expected = fmpq_mat(2, 2)
            for pattern in product((0, 1), repeat=q):
                weight = fmpq(15**sum(pattern), 16**q)
                expected += regional[sum(pattern)] * weight
            self.assertEqual(actual, expected)

    def test_missing_occupancy_rejected(self):
        with self.assertRaises(ValueError):
            screen.regional_uniform([fmpq_mat([[1]])], 1)


if __name__ == '__main__':
    unittest.main()
