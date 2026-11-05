"""Tiny exact labelled controls for the two-weight h2 reduction."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_gf16_dense_two import (arity_pairs, aggregate_moments, two_weight_vectors,
                                      pair_branch_caps, run)
from hd_two_sum_gf16_arity_dense import arity_weights
from hd_two_sum_gf16_full_dense import column_vectors
from hd_two_sum_affine_outer import uniform_mask_counts
from hd_xor_dense_holder_classes import class_log_weights
from hd_xor_dense_two_feature_tests import exact_tables, logs, action
from hd_xor_joint_columns import categories


class Gf16DenseTwoTests(unittest.TestCase):
    def test_direct_labelled_noncommuting_masks_and_holes(self):
        tables, potential, q = exact_tables(length=3, maximum=4), F(2, 3), 7
        actual = two_weight_vectors(logs(tables), np.log2(float(potential)), 2, q=q)
        for index, (single, double) in enumerate(categories(2)):
            expected = [[[F(0), F(0)] for _ in range(4)] for _ in range(4)]
            for left, right in product(range(8), repeat=2):
                for holes in product(range(-1, 3), repeat=double):
                    mass = F((q-4)**holes.count(-1), (q-1)**double)
                    occupancy = [single+double-holes.count(i)+((left >> i) & 1)+((right >> i) & 1)
                                 for i in range(3)]
                    vector = action(tables, occupancy, potential)
                    for state in range(2):
                        expected[left.bit_count()][right.bit_count()][state] += mass*vector[state]
            np.testing.assert_allclose(np.exp2(actual[index]), np.asarray(expected, dtype=float), rtol=5e-14)

    def test_scalar_identity_and_group_swap_symmetry(self):
        identity = np.tile([[0., -np.inf], [-np.inf, 0.]], (4, 5, 1, 1))
        actual = two_weight_vectors(identity, -2., 2, q=7)
        expected = np.array([[[comb(4, a)*comb(4, b), comb(4, a)*comb(4, b)/4]
                              for b in range(5)] for a in range(5)])
        for row in actual:
            np.testing.assert_allclose(np.exp2(row), expected, rtol=3e-14)
        self.assertTrue(np.allclose(actual, actual.swapaxes(1, 2), atol=1e-14))

    def test_mask_weights_mixed_before_potential_and_marker_split(self):
        vectors = np.log2(np.array([[[[2., 1.], [1., 10.]], [[10., 1.], [3., 2.]]]]))
        caps = np.exp2(pair_branch_caps(vectors, 0., np.log2([1., 2.]), np.log2([3., 5.])))[0]
        np.testing.assert_allclose(caps, [6., 65., 30.])
        # The two nonempty alternatives are mixed before max:65, not60+50.
        self.assertLess(caps[1], 110.)

    def test_full_arity_against_independent_seen_count(self):
        tables, logp = logs(exact_tables(length=16, maximum=4)), -1.
        actual = two_weight_vectors(tables, logp, 2)
        weights, _ = arity_weights(16)
        caps = pair_branch_caps(actual, logp, weights, weights)
        old = np.max(column_vectors(tables, logp, 2, 2)-[0., logp], axis=2)-16*np.log2(255)
        np.testing.assert_allclose(caps, old, atol=4e-13)

    def test_105_classes_partition_ordered_rowsets_once(self):
        pairs = arity_pairs()
        self.assertEqual(len(pairs), 105)
        self.assertEqual(sum(mult*comb(16, a)*comb(16, b) for a, b, mult in pairs), 65399**2)
        self.assertEqual(pairs[-1], (16, 16, 1))
        self.assertEqual(next(mult for a, b, mult in pairs if (a, b) == (15, 16)), 2)

    def test_aggregate_moments_against_exact_g8_labelled_sums(self):
        nums, denominator = aggregate_moments(8, 2, 1)
        counts = [(active, uniform_mask_counts(active)) for active in range(256) if active.bit_count() >= 3]
        for target in range(256):
            expected = sum((F(values[target]**2, 255**active.bit_count()) for active, values in counts), F(0))
            self.assertEqual(F(nums[target.bit_count()], denominator), expected)
        weights = class_log_weights(nums, denominator, 2)
        self.assertFalse(np.any(np.isnan(weights)))

    def test_zero_matrices_preserve_zero(self):
        vectors = two_weight_vectors(np.full((3, 4, 2, 2), -np.inf), -1., 1, q=7)
        caps = pair_branch_caps(vectors, -1., np.zeros(4), np.zeros(4))
        self.assertTrue(np.all(np.isneginf(caps)))

    def test_scope_and_timeout_before_regional_bank(self):
        with patch('hd_two_sum_gf16_dense_two.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertEqual(result['cohort_counts'], [0, 15, 2])
        self.assertEqual(result['source_colorings'], 136)
        self.assertEqual(result['sparse_pair_assignments'], 120**15)
        self.assertEqual(result['sparse_normalizer_exponent'], 240)
        self.assertEqual(result['dense_nonempty_column_floor'], 18)
        self.assertFalse(result['complete_selected_cohort_and_source'])
        self.assertFalse(result['directed'])
        self.assertNotIn('failure_log2_screen', result)

    def test_invalid_inputs_rejected(self):
        tables = logs(exact_tables(length=3, maximum=4))
        for invalid in (True, -1, 16):
            with self.assertRaises(ValueError):
                two_weight_vectors(tables, -1., invalid)
        for invalid in (True, 0, 121, float('inf')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)
        with self.assertRaises(ValueError):
            pair_branch_caps(np.zeros((1, 2, 2, 2)), 0., [0., np.nan], [0., 0.])


if __name__ == '__main__':
    unittest.main()
