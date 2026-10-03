"""Small exact-count comparisons for the bounded-memory tuple merge."""
from itertools import product
from math import comb
import unittest

import numpy as np

import long_grouped_gate as long


class LongGroupedTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, 0.], [0., .3]],
            [[.15, .25], [.03, .08]], [[.08, .02], [.09, .04]]])

    def test_matches_existing_fine_g4(self):
        for alpha in (.35, .7, 1.):
            actual, metadata = long.long_operators(self.local, 4, alpha)
            expected = long.fine.fine_operators(long.fine.tuple_products(self.local, 4), alpha)
            np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=1e-15)
            self.assertLess(metadata['alpha1_max_relative_error'], 3e-13)

    def test_eight_step_merge_matches_explicit_binary_counts(self):
        # Each physical step now has one potential slot, hence only256 tuples.
        local = self.local[:2]
        actual, metadata = long.long_operators(local, 8, .35)
        expected = np.zeros_like(actual)
        for counts in product(range(2), repeat=8):
            matrix = np.eye(2)
            for j in counts:
                matrix = matrix@local[j]
            expected[sum(counts)] += matrix**.35/comb(8, sum(counts))
        np.testing.assert_allclose(actual, expected, rtol=5e-13, atol=1e-15)
        self.assertEqual(metadata['tuple_count'], 256)

    def test_vandermonde_mass_normalization(self):
        local = np.ones((3, 1, 1))
        for group in (2, 4, 8):
            actual, _ = long.long_operators(local, group, .35)
            np.testing.assert_allclose(actual, 1., rtol=2e-14, atol=2e-14)

    def test_invalid_parameters(self):
        for group in (0, 3, 16):
            with self.assertRaises(ValueError):
                long.long_operators(self.local, group, .35)


if __name__ == '__main__':
    unittest.main()
