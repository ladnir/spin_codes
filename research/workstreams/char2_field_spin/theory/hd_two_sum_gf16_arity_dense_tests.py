"""Exact tiny and labelled-g8 controls; no production regional bank in tests."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_gf16_arity_dense import weight_vectors, arity_weights, run
from hd_two_sum_affine_dense import dense_mask_vectors, affine_class_moments
from hd_two_sum_gf16_full_dense import column_vectors
from hd_xor_dense_holder_classes import mask_branch_caps, class_log_weights
from hd_xor_dense_two_feature_tests import exact_tables, logs, action
from hd_xor_joint_columns import categories


class Gf16ArityDenseTests(unittest.TestCase):
    def test_independent_labelled_holes_and_weight_mixture(self):
        tables, q, potential = exact_tables(length=3, maximum=3), 7, F(2, 3)
        actual = weight_vectors(logs(tables), np.log2(float(potential)), 2, q=q)
        for index, (single, double) in enumerate(categories(2)):
            expected = [[F(0), F(0)] for _ in range(4)]
            for mask in range(8):
                for holes in product(range(-1, 3), repeat=double):
                    mass = F((q-4)**holes.count(-1), (q-1)**double)
                    occupancy = [single+double-holes.count(i)+((mask >> i) & 1) for i in range(3)]
                    vector = action(tables, occupancy, potential)
                    for state in range(2):
                        expected[mask.bit_count()][state] += mass*vector[state]
            np.testing.assert_allclose(np.exp2(actual[index]), np.asarray(expected, dtype=float), rtol=3e-14)

    def test_scalar_identity_counts_each_mask_once(self):
        table = np.tile([[0., -np.inf], [-np.inf, 0.]], (8, 4, 1, 1))
        actual = np.exp2(weight_vectors(table, -1., 2))
        expected = np.array([[comb(8, w), comb(8, w)/2] for w in range(9)])
        for row in actual:
            np.testing.assert_allclose(row, expected, rtol=2e-14)

    def test_g8_labelled_reproduction_all_arity_caps(self):
        tables, logp = logs(exact_tables(length=8, maximum=3)), -2.
        labelled = dense_mask_vectors(tables, logp, sparse_groups=2)
        reduced = weight_vectors(tables, logp, sparse_groups=2)
        for w in range(9):
            selected = [mask for mask in range(256) if mask.bit_count() == w]
            np.testing.assert_allclose(reduced[:, w], np.logaddexp2.reduce(labelled[:, selected], axis=1), atol=3e-14)
        for arity in range(3, 9):
            masks = tuple(mask for mask in range(256) if mask.bit_count() == arity)
            nums, denominator = affine_class_moments(masks, columns=2, dimension=1)
            old = mask_branch_caps(labelled, logp, class_log_weights(nums, denominator, 2))
            weights, _ = arity_weights(arity, 8, 2, 1)
            new = mask_branch_caps(reduced, logp, weights)
            for a, b in zip(old, new):
                np.testing.assert_allclose(a, b, atol=4e-14)

    def test_full_arity_agrees_with_independent_seen_count_dp(self):
        tables, logp = logs(exact_tables(length=16, maximum=3)), -1.
        vectors = weight_vectors(tables, logp, 2)
        weights, _ = arity_weights(16)
        actual = mask_branch_caps(vectors, logp, weights)
        seen = column_vectors(tables, logp, 2, 1)
        expected = np.max(seen-[0., logp], axis=2)-8*np.log2(255)
        np.testing.assert_allclose(np.column_stack(actual), expected, atol=3e-13)

    def test_exact_zero_weights_remain_negative_infinity(self):
        tables = np.full((3, 4, 2, 2), -np.inf)
        actual = weight_vectors(tables, -1., 2, q=7)
        self.assertTrue(np.all(np.isneginf(actual)))
        weights, _ = arity_weights(3)
        self.assertTrue(np.any(np.isneginf(weights)))
        self.assertFalse(np.any(np.isnan(weights)))

    def test_arity_partition_scope_and_timeout_before_bank(self):
        with patch('hd_two_sum_gf16_arity_dense.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertEqual(result['arities'], list(range(3, 17)))
        self.assertEqual(result['fixed_dense_row_sets'], 65399)
        self.assertEqual(result['cohort_counts'], [0, 16, 1])
        self.assertEqual(result['sparse_normalizer_exponent'], 256)
        self.assertFalse(result['extra_dense_rowset_factor'])
        self.assertFalse(result['extra_target_mask_binomial'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['complete_selected_cohort_and_source'])
        self.assertNotIn('failure_log2_screen', result)

    def test_bad_dimensions_and_numerics_rejected(self):
        table = logs(exact_tables(length=3, maximum=3))
        for invalid in (True, -1, 17):
            with self.assertRaises(ValueError):
                weight_vectors(table, -1., invalid)
        for invalid in (True, 0, 121, float('nan')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)
        table[0, 0, 0, 0] = np.inf
        with self.assertRaises(ValueError):
            weight_vectors(table, -1., 2)


if __name__ == '__main__':
    unittest.main()
