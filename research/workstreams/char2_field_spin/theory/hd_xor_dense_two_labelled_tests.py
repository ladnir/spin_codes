"""Exact small labelled-mask controls; no production operator bank."""
from fractions import Fraction as F
from itertools import product
import unittest
from unittest.mock import patch

import numpy as np

from hd_xor_dense_two_feature_tests import exact_tables, action, logs
from hd_xor_dense_two_labelled import (occupancy_convolution, reachable_sums,
                                      offset_vectors, labelled_caps, evaluate)
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_floor import offset_resource_vectors


def direct_full_dp(tables, potential, lower, upper, offset):
    """Independent positive rational recurrence for one complete offset."""
    width = upper-lower
    values = [(F(1), potential)]
    for region in reversed(range(len(tables))):
        following = [(F(0), F(0)) for _ in range(len(values)+width)]
        for j, vector in enumerate(values):
            for delta in range(width+1):
                matrix = tables[region][lower+delta+offset[region]]
                candidate = tuple(sum(matrix[i][k]*vector[k] for k in range(2)) for i in range(2))
                following[j+delta] = tuple(max(following[j+delta][i], candidate[i]) for i in range(2))
        values = following
    return values


def digits(code, base, width):
    return tuple((code // base**i) % base for i in range(width))


class DenseTwoLabelledTests(unittest.TestCase):
    def test_exact_whole_mask_pair_convolution_and_nonempty_flags(self):
        weights = [F((5*i+3) % 13, 7) for i in range(16)]
        result = occupancy_convolution(weights, 2)
        expected = [[F(0) for _ in range(81)] for _ in range(3)]
        for left, right in product(range(16), repeat=2):
            occupancy = tuple(((left >> i) & 1)+((right >> i) & 1) for i in range(4))
            code = sum(value*3**i for i, value in enumerate(occupancy))
            expected[int(left != 0)+int(right != 0)][code] += weights[left]*weights[right]
        self.assertEqual(result, expected)
        self.assertEqual(sum(map(sum, result)), sum(weights)**2)
        self.assertEqual(result[0][0], weights[0]**2)
        # Code1 permits e1 and e2 only when distinct nonempty masks can
        # partition its support; one occupied coordinate forces e1.
        self.assertEqual(result[1][1], 2*weights[0]*weights[1])
        self.assertEqual(result[2][1], 0)
        self.assertEqual(result[2][2], weights[1]**2)

    def test_reachability_matches_all_suffix_completions(self):
        for width in range(4):
            targets = tuple(sorted({3*width-k for k in range(width+1)}))
            keep = reachable_sums(4, width, targets)
            for stage in range(5):
                expected = {sum(prefix) for prefix in product(range(width+1), repeat=stage)
                            if any(sum(prefix)+sum(suffix) in targets
                                   for suffix in product(range(width+1), repeat=4-stage))}
                self.assertEqual(set(keep[stage]), expected)

    def test_shared_ternary_suffixes_equal_full_exact_noncommuting_dp(self):
        tables, potential = exact_tables(maximum=4), F(2, 3)
        targets = [0, 2, 3, 6, 8]
        result, stats = offset_vectors(logs(tables), np.log2(float(potential)), 0, 2, 2, targets)
        self.assertEqual(stats['target_sums'], targets)
        for code in range(81):
            offset = digits(code, 3, 4)
            exact = direct_full_dp(tables, potential, 0, 2, offset)
            np.testing.assert_allclose(np.exp2(result[code]),
                                       np.array([exact[j] for j in targets], dtype=float), rtol=2e-14)
            for occupancy in product(range(3), repeat=4):
                if sum(occupancy) in targets:
                    vector = action(tables, tuple(occupancy[i]+offset[i] for i in range(4)), potential)
                    bound = np.exp2(result[code, targets.index(sum(occupancy))])
                    self.assertTrue(np.all(bound*(1+1e-12) >= np.asarray(vector, dtype=float)))

    def test_h1_same_recurrence_equals_existing_unpruned_mask_bank(self):
        tables = logs(exact_tables(maximum=3))
        for lower, upper in ((0, 0), (0, 1), (0, 2), (1, 2)):
            targets = sorted({0, upper-lower, 3*(upper-lower), 4*(upper-lower)})
            result, _ = offset_vectors(tables, -1., lower, upper, 1, targets)
            reference = offset_resource_vectors(tables, -1., lower, upper)
            np.testing.assert_array_equal(result, reference[:, targets, :])

    def test_complete_labelled_mixture_matches_direct_mask_pair_sum(self):
        tables, potential = exact_tables(maximum=3), F(2, 3)
        weights = [F((3*i+2) % 7+1, 11) for i in range(16)]
        mass = np.array(occupancy_convolution(weights, 2), dtype=float)
        with np.errstate(divide='ignore'):
            caps, _ = labelled_caps(logs(tables), np.log2(float(potential)), np.log2(mass), 1, 2, q=4)
        for index, (single, double) in enumerate(categories(1)):
            total = [[F(0), F(0)] for _ in range(3)]
            for left, right in product(range(16), repeat=2):
                e = int(left != 0)+int(right != 0)
                offset = tuple(((left >> i) & 1)+((right >> i) & 1) for i in range(4))
                for k in range(double+1):
                    probability = F(2 if k == 0 else 1, 3) if double else F(1)
                    vectors = direct_full_dp(tables, potential, double-k, single+double, offset)
                    target = 3*single+2*k
                    for state in range(2):
                        total[e][state] += weights[left]*weights[right]*probability*vectors[target][state]
            expected = [max(vector[0], vector[1]/potential) for vector in total]
            np.testing.assert_allclose(np.exp2(caps[index]), np.asarray(expected, dtype=float), rtol=3e-14)

    def test_budget_failure_is_recorded_without_a_bound_or_restart(self):
        prefix = 'hd_xor_dense_two_labelled.'
        with patch(prefix+'pair_regions_screen', return_value=None), \
                patch(prefix+'ordered_region_tables', return_value=None), \
                patch(prefix+'labelled_caps') as bank:
            result = evaluate(2, seconds=1e-12)
        bank.assert_not_called()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertNotIn('failure_log2_screen', result)
        self.assertEqual(result['cohort_counts'], [0, 15, 2])
        self.assertEqual(result['source_colorings'], 136)
        self.assertEqual(result['sparse_pair_assignments'], 28**15)
        self.assertEqual(result['sparse_normalizer_exponent'], 480)
        self.assertFalse(result['directed'])
        with self.assertRaises(ValueError):
            evaluate(3)


if __name__ == '__main__':
    unittest.main()
