"""Small independent normalization and chronology tests for mixed proposals."""
from fractions import Fraction
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_xor_joint_columns import categories, joint_log_count
from hd_xor_mixed_gate import (dense_density_log2, integrate_dense_regions,
    mixed_column_log_coefficients, ordered_region_tables, profile_counts,
    source_coloring_count, sparse_log_count)
from hd_xor_paired_joint_gate import ordered_tables
from hd_xor_column_screen import group_density_bits


class MixedGateTests(unittest.TestCase):
    def test_mixed_column_multiplicities_match_all_binary_supports(self):
        for n1, n2 in ((1, 1), (2, 2), (0, 3), (3, 0), (0, 0)):
            caps = {(s, d): -(s+3*d)/4 for s, d in categories(n1+n2)}
            expected = np.zeros(n1+2*n2+1)
            for bits in product(range(2), repeat=n1+2*n2):
                sizes = [bits[n1+2*i]+bits[n1+2*i+1] for i in range(n2)]
                single = sum(bits[:n1])+sizes.count(1)
                expected[sum(bits)] += 2**caps[single, sizes.count(2)]
            np.testing.assert_allclose(np.exp2(mixed_column_log_coefficients(caps, n1, n2)), expected)

    def test_sparse_tail_matches_exhaustive_toy_arrays(self):
        n1, n2, n, k, q = 1, 1, 3, 2, 4
        caps = {(s, d): Fraction(1, 3**s*7**d) for s, d in categories(2)}
        exact = Fraction(0)
        for columns in product(range(8), repeat=n):
            degree = sum(column.bit_count() for column in columns)
            if degree < 6:
                continue
            value = Fraction(3**(degree-3))
            for column in columns:
                one = column & 1
                pair_size = ((column >> 1) & 3).bit_count()
                value *= caps[one+(pair_size == 1), int(pair_size == 2)]
            exact += value
        result = sparse_log_count({key: log2(value) for key, value in caps.items()}, n1, n2, n, k, q)
        self.assertAlmostEqual(result, log2(exact), places=12)

    def test_two_row_only_reduces_to_existing_joint_tail(self):
        caps = {(s, d): -3*s-4*d for s, d in categories(3)}
        self.assertAlmostEqual(sparse_log_count(caps, 0, 3), joint_log_count(caps, 3), places=10)
        self.assertAlmostEqual(sparse_log_count({(0, 0): -3.}, 0, 0), -96.)

    def test_dense_mixture_matches_labelled_bernoulli_subsets(self):
        eta = .3*255/256
        matrices = np.arange(1, 8*5*4+1, dtype=float).reshape(8, 5, 2, 2)/1024
        mixed = integrate_dense_regions(np.log2(matrices), 2, 2, .3)
        expected = np.zeros((8, 3, 2, 2))
        for bits in product(range(2), repeat=2):
            occupied = sum(bits)
            probability = eta**occupied*(1-eta)**(2-occupied)
            expected += probability*matrices[:, occupied:occupied+3]
        np.testing.assert_allclose(np.exp2(mixed), expected, rtol=1e-14)
        np.testing.assert_allclose(integrate_dense_regions(np.log2(matrices), 4, 0, .3), np.log2(matrices))

    def test_source_geometry_preserves_root_labelled_order(self):
        raw = {(m, f): np.array([[-m-f, -2*m-f], [-3*m+f, -4*m+f]], dtype=float)
               for m in range(18) for f in range(m//2+1)}
        for a, b in ((0, 8), (8, 0), (3, 5)):
            np.testing.assert_array_equal(ordered_region_tables(raw, 17, a, b), ordered_tables(raw, a, b))
        self.assertFalse(np.array_equal(ordered_region_tables(raw, 17, 0, 8), ordered_region_tables(raw, 17, 8, 0)))

    def test_cohort_and_dense_type_choices_are_counted_once(self):
        self.assertEqual(source_coloring_count((1, 16, 0)), 17)
        self.assertEqual(source_coloring_count((0, 16, 1)), 17)
        self.assertEqual(source_coloring_count((2, 3, 4)), comb(9, 2)*comb(7, 3))
        p = .98
        expected = np.logaddexp2.reduce([group_density_bits(active, p) for active in range(3, 9)])
        self.assertEqual(dense_density_log2(p), expected)
        with self.assertRaises(ValueError):
            profile_counts((1, 17, 0))


if __name__ == '__main__':
    unittest.main()
