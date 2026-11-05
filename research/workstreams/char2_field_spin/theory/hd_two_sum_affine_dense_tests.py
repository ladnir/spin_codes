"""Exact/light affine dense-gate controls, without a production screen."""
from fractions import Fraction as F
from itertools import product
from math import comb, prod
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_affine_outer import uniform_mask_counts, column_action_exact
from hd_two_sum_affine_dense import (ALL_DENSE, ARITY_CLASSES, affine_class_moments,
    offset_hole_suffix, dense_mask_vectors, run)
from hd_xor_dense_two_feature_tests import exact_tables, action, logs
from hd_xor_joint_columns import categories


class AffineDenseTests(unittest.TestCase):
    def test_new_affine_kernel_not_old_xor_singleton_or_pair(self):
        for coordinate in range(8):
            expected = [0]*256
            expected[0], expected[255] = 1, 255
            self.assertEqual(uniform_mask_counts(1 << coordinate), tuple(expected))
        for active in (3, 5, 129, 192):
            expected = [0]*256
            expected[0], expected[255] = 1, 255*249
            for region in range(8):
                expected[255 ^ (1 << region)] = 255
            self.assertEqual(uniform_mask_counts(active), tuple(expected))
        self.assertEqual(sum(uniform_mask_counts(7)), 256**3)
        self.assertEqual(uniform_mask_counts(7)[0], 1)

    def test_exact_affine_moments_and_fixed_class_holder(self):
        members = (1, 3, 7)
        nums, denominator = affine_class_moments(members, columns=2, dimension=1)
        for target in range(256):
            expected = sum(F(uniform_mask_counts(a)[target]**2, 255**a.bit_count()) for a in members)
            self.assertEqual(F(nums[target], denominator), expected)
        for left, right in product((0, 7, 127, 254, 255), repeat=2):
            fixed_sum = sum(F(uniform_mask_counts(a)[left]*uniform_mask_counts(a)[right],
                              255**a.bit_count()) for a in members)
            self.assertLessEqual(fixed_sum**2, F(nums[left]*nums[right], denominator**2))
        # The affine singleton has full support, unlike XOR2's missing diagonal.
        nums, denominator = affine_class_moments((1,), columns=2, dimension=1)
        self.assertEqual(F(nums[255], denominator), 255)
        self.assertEqual(nums[254], 0)

    def test_labelled_dense_offsets_equal_independent_exact_hole_recurrence(self):
        tables, potential = exact_tables(maximum=3), F(2, 3)
        result = dense_mask_vectors(logs(tables), np.log2(float(potential)), sparse_groups=2, q=8)
        for index, (single, double) in enumerate(categories(2)):
            for mask in range(16):
                shifted = [tuple(tuple(value for row in table[m+((mask >> i) & 1)] for value in row)
                                 for m in range(3)) for i, table in enumerate(tables)]
                exact = column_action_exact(shifted, (F(1), potential), single, double, q=8)
                np.testing.assert_allclose(np.exp2(result[index, mask]), np.asarray(exact, dtype=float), rtol=2e-14)

    def test_direct_labelled_hole_assignments_noncommuting(self):
        tables, potential = exact_tables(maximum=3), F(2, 3)
        result = dense_mask_vectors(logs(tables), np.log2(float(potential)), sparse_groups=2, q=8)
        index = categories(2).index((0, 2))
        for mask in range(16):
            total = [F(0), F(0)]
            for holes in product(range(-1, 4), repeat=2):
                probability = prod(F(3 if hole == -1 else 1, 7) for hole in holes)
                occupancy = tuple(2-holes.count(i)+((mask >> i) & 1) for i in range(4))
                vector = action(tables, occupancy, potential)
                for state in range(2):
                    total[state] += probability*vector[state]
            np.testing.assert_allclose(np.exp2(result[index, mask]), np.asarray(total, dtype=float), rtol=2e-14)
        identity = np.tile(np.array([[0., -np.inf], [-np.inf, 0.]]), (4, 4, 1, 1))
        normalized = dense_mask_vectors(identity, -1., sparse_groups=2, q=8)
        np.testing.assert_allclose(normalized, np.broadcast_to([0., -1.], normalized.shape), atol=2e-15)

    def test_arity_partition_and_once_only_scope(self):
        self.assertEqual(len(ALL_DENSE), 219)
        self.assertEqual(sorted(mask for _, masks in ARITY_CLASSES for mask in masks), list(ALL_DENSE))
        self.assertEqual([len(masks) for _, masks in ARITY_CLASSES], [comb(8, a) for a in range(3, 9)])
        namespace = 'hd_two_sum_affine_dense.'
        with patch(namespace+'affine_class_moments', return_value=((1,)*256, 1)), \
                patch(namespace+'pair_regions_screen', return_value=None), \
                patch(namespace+'ordered_region_tables', return_value=None), \
                patch(namespace+'dense_mask_vectors', return_value=np.zeros((153, 256, 2))), \
                patch(namespace+'bivariate_log_tail', return_value=(0., [dict(log2_mass=0.,
                    nonempty_dense_columns=17, largest_sparse_degree=544)])) as tail:
            result = run()
        self.assertEqual(tail.call_count, 6)
        self.assertEqual(result['cohort_counts'], [0, 16, 1])
        self.assertEqual(result['source_matching_counts'], [8, 0])
        self.assertEqual(result['sparse_pair_assignments'], 28**16)
        self.assertEqual(result['source_colorings'], 17)
        self.assertEqual(result['fixed_row_sets'], 219)
        self.assertEqual(result['sparse_normalizer_exponent'], 512)
        self.assertEqual(result['labels'], list(range(8, 16)))
        self.assertFalse(result['old_xor_mask_counts_used'])
        self.assertFalse(result['old_certificates_transferred'])
        self.assertFalse(result['directed'])
        for call in tail.call_args_list:
            self.assertEqual(call.kwargs, dict(columns=32, x_floor=544, y_floor=17))


if __name__ == '__main__':
    unittest.main()
