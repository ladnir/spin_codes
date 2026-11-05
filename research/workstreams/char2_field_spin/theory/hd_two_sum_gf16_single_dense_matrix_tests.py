"""Tiny independent ordered-matrix controls; no production regional bank."""
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_gf16_single_dense_matrix import (log_matmul, weight_matrices,
                                               matrix_branches, matrix_tail, run)
from hd_xor_dense_two_feature_tests import exact_tables, logs


class SingleDenseMatrixTests(unittest.TestCase):
    def test_weight_coefficients_match_all_labelled_products(self):
        tables = np.asarray(exact_tables(length=4, maximum=1), dtype=object)
        actual = np.exp2(weight_matrices(logs(tables)))
        expected = np.zeros((5, 2, 2), dtype=object)
        for bits in product((0, 1), repeat=4):
            matrix = np.eye(2, dtype=object)
            for i, bit in enumerate(bits):
                matrix = matrix @ tables[i, bit]
            expected[sum(bits)] += matrix
        np.testing.assert_allclose(actual, np.asarray(expected, dtype=float), rtol=5e-15)
        backward = np.exp2(weight_matrices(logs(tables[::-1])))
        self.assertGreater(float(np.max(np.abs(backward-actual))), 1.)

    def test_matrix_polynomial_matches_every_noncommuting_word(self):
        branches = np.array([[[2, 1], [1, 3]], [[1, 4], [2, 1]]], dtype=object)
        expected = 0
        for bits in product((0, 1), repeat=5):
            if sum(bits) < 3:
                continue
            matrix = np.eye(2, dtype=object)
            for bit in bits:
                matrix = matrix @ branches[bit]
            expected += matrix[0, 0]+matrix[0, 1]
        actual, _ = matrix_tail(np.log2(np.asarray(branches, dtype=float)), 5, 3)
        self.assertAlmostEqual(actual, log2(expected), places=12)

    def test_scalar_identity_polynomial_and_initial_state(self):
        with np.errstate(divide='ignore'):
            branches = np.log2(np.array([2*np.eye(2), 3*np.eye(2)]))
        actual, _ = matrix_tail(branches, 5, 3)
        expected = sum(comb(5, j)*2**(5-j)*3**j for j in range(3, 6))
        self.assertAlmostEqual(actual, log2(expected), places=12)
        # Huge inaccessible live-state entries cannot be charged to initial0.
        branches[:, 1, 1] = 100.
        again, _ = matrix_tail(branches, 5, 3)
        self.assertEqual(actual, again)

    def test_direct_tail_below_shared_potential_majorant(self):
        branches = np.array([[[2., 1.], [.1, 3.]], [[.3, 4.], [2., .2]]])
        potential = np.array([1., .25])
        caps = np.max((branches @ potential)/potential, axis=1)
        majorant = sum(comb(5, j)*caps[0]**(5-j)*caps[1]**j for j in range(3, 6))/min(potential)
        actual, _ = matrix_tail(np.log2(branches), 5, 3)
        self.assertLessEqual(actual, log2(majorant)+1e-12)

    def test_mixing_weights_before_matrix_power(self):
        matrices = np.log2(np.array([[[1., 2.], [3., 4.]],
                                     [[2., 1.], [4., 3.]], [[5., 6.], [7., 8.]]]))
        actual = np.exp2(matrix_branches(matrices, np.log2([2., 3., 5.])))
        expected = np.array([2*np.exp2(matrices[0]), 3*np.exp2(matrices[1])+5*np.exp2(matrices[2])])
        np.testing.assert_allclose(actual, expected, rtol=2e-15)

    def test_zero_operators_and_bad_inputs(self):
        value, _ = matrix_tail(np.full((2, 2, 2), -np.inf), 3, 1)
        self.assertEqual(value, -np.inf)
        with self.assertRaises(ValueError):
            log_matmul(np.full((2, 2), np.nan), np.eye(2))
        with self.assertRaises(ValueError):
            matrix_tail(np.zeros((2, 2, 2)), 3, 4)

    def test_scope_budget_and_no_endpoint_potential(self):
        with patch('hd_two_sum_gf16_single_dense_matrix.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(result['complete_selected_cohort'])
        self.assertEqual(result['source_supports'], 8192)
        self.assertEqual(result['endpoint_potential_factors'], 0)
        self.assertEqual(result['terminal_vector'], [1, 1])
        self.assertEqual(result['initial_state'], 0)
        self.assertFalse(result['common_potential_used'])
        self.assertFalse(result['directed'])
        for invalid in (True, 0, 21, float('inf')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)


if __name__ == '__main__':
    unittest.main()
