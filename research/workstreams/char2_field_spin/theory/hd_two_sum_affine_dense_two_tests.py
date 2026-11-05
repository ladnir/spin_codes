"""Exact/light ternary affine h2 controls; no production operator run."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_affine_dense import offset_hole_suffix as binary_suffix
from hd_two_sum_affine_outer import column_action_exact
from hd_two_sum_affine_dense_two import (two_mask_masses, offset_hole_suffix,
    occupancy_vectors, arity_pairs, branch_caps, run_one)
from hd_xor_dense_two_feature_tests import exact_tables, logs
from hd_xor_joint_columns import categories


class AffineDenseTwoTests(unittest.TestCase):
    def test_exact_asymmetric_mask_pair_convolution(self):
        left = [F((3*i+2) % 7, 5) for i in range(16)]
        right = [F((5*i+3) % 11, 13) for i in range(16)]
        actual = two_mask_masses(left, right)
        expected = [[F(0)]*81 for _ in range(3)]
        for a, b in product(range(16), repeat=2):
            code = sum((((a >> i) & 1)+((b >> i) & 1))*3**i for i in range(4))
            expected[int(a != 0)+int(b != 0)][code] += left[a]*right[b]
        self.assertEqual(actual, expected)
        self.assertEqual(sum(map(sum, actual)), sum(left)*sum(right))
        self.assertEqual(actual, two_mask_masses(right, left))

    def test_h1_specialization_reproduces_existing_affine_suffix(self):
        tables = logs(exact_tables(maximum=3))
        for occupied in range(3):
            np.testing.assert_array_equal(offset_hole_suffix(tables, -1., occupied, dense_groups=1),
                                          binary_suffix(tables, -1., occupied))

    def test_all_ternary_offsets_match_exact_noncommuting_hole_law(self):
        tables, potential = exact_tables(maximum=4), F(2, 3)
        actual = occupancy_vectors(logs(tables), np.log2(float(potential)), sparse_groups=2, q=8)
        for index, (single, double) in enumerate(categories(2)):
            for code in range(81):
                offset = [(code // 3**i) % 3 for i in range(4)]
                shifted = [tuple(tuple(value for row in table[m+offset[i]] for value in row)
                                 for m in range(3)) for i, table in enumerate(tables)]
                exact = column_action_exact(shifted, (F(1), potential), single, double, q=8)
                np.testing.assert_allclose(np.exp2(actual[index, code]), np.asarray(exact, dtype=float), rtol=3e-14)

    def test_masks_are_mixed_before_one_potential_quotient(self):
        vectors = np.log2(np.array([[[10., 1.], [1., 10.], [2., 3.]]]))
        with np.errstate(divide='ignore'):
            masses = np.log2(np.array([[1., 0., 0.], [0., 0., 1.], [2., 5., 0.]]))
        actual = np.exp2(branch_caps(vectors, masses, 0.))[0]
        np.testing.assert_allclose(actual, [10., 3., 52.])
        self.assertLess(actual[2], 70.)

    def test_unordered_arity_types_count_ordered_rowsets_once(self):
        pairs = arity_pairs()
        self.assertEqual(len(pairs), 21)
        self.assertEqual(sum(mult*comb(8, a)*comb(8, b) for a, b, mult in pairs), 219**2)
        self.assertEqual(next(mult for a, b, mult in pairs if (a, b) == (8, 8)), 1)
        self.assertEqual(next(mult for a, b, mult in pairs if (a, b) == (7, 8)), 2)

    def test_bounded_scope_and_timeout_without_claim(self):
        namespace = 'hd_two_sum_affine_dense_two.'
        with patch(namespace+'pair_regions_screen', return_value=None), \
                patch(namespace+'ordered_region_tables', return_value=None), \
                patch(namespace+'occupancy_vectors') as bank:
            result = run_one(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertNotIn('failure_log2_screen', result)
        self.assertEqual(result['source_colorings'], 136)
        self.assertEqual(result['sparse_pair_assignments'], 28**15)
        self.assertEqual(result['sparse_normalizer_exponent'], 480)
        self.assertEqual(result['dense_nonempty_column_floor'], 34)
        self.assertFalse(result['complete_dense_rowset_coverage'])
        self.assertFalse(result['directed'])
        with self.assertRaises(ValueError):
            run_one(7, 8)


if __name__ == '__main__':
    unittest.main()
