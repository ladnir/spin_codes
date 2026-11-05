"""Bounded exact/overlap tests for degree1152; no controller is run."""
import unittest

import numpy as np

import hd_pair_state_tail as old
import hd_pair_state_tail_general as general


class GeneralTailTests(unittest.TestCase):
    def test_a1_a2_bitwise_overlap_with_frozen_engine(self):
        rng = np.random.default_rng(27)
        for maximum in (2, 4):
            weights = np.array([list(range(maximum+1))]*2)
            counts = rng.uniform(0, 8, weights.shape)
            blocks = [rng.uniform(-40, 2, weights.shape+(2, 2)) for _ in range(2)]
            for exponent in (1, 3):
                a = old.word_polynomial(blocks, counts, weights, exponent)
                b = general.word_polynomial(blocks, counts, weights, exponent)
                np.testing.assert_array_equal(a, b)
                for floor in (0, 3, 2*maximum*exponent, 2*maximum*exponent+1):
                    np.testing.assert_array_equal(old.tail_from_polynomial(a, floor),
                                                  general.tail_from_polynomial(b, floor))
        self.assertEqual(old.MAX_DEGREE, 576)
        self.assertEqual(old.MAX_BRANCH_WEIGHT, 4)

    def test_weight8_and_degree1152_survive_exact_zero_tail(self):
        # One monomial and one state path have a known answer. This checks
        # the whole high-degree interface without evaluating a real bank.
        block = np.full((1, 2, 2), -np.inf)
        block[0, 0, 0] = -18_000.
        polynomial = general.word_polynomial([block]*8, [0.], np.array([8]), 18)
        self.assertEqual(polynomial.shape, (1153, 2, 2))
        self.assertTrue(np.all(np.isneginf(polynomial[:1152])))
        self.assertAlmostEqual(float(polynomial[1152, 0, 0]), -144_000., places=8)
        self.assertAlmostEqual(float(general.tail_from_polynomial(polynomial, 577)), -144_000., places=8)
        self.assertTrue(np.isneginf(general.tail_from_polynomial(polynomial, 1153)))

    def test_high_degree_ordered_noncommuting_product(self):
        left = np.full((578, 2, 2), -np.inf)
        right = np.full((2, 2, 2), -np.inf)
        a, b = np.array([[1., 2.], [3., 1.]]), np.array([[2., 0.], [1., 4.]])
        with np.errstate(divide='ignore'):
            left[577] = np.log2(a)
            right[1] = np.log2(b)
        forward = general.matrix_polynomial_product(left, right)
        backward = general.matrix_polynomial_product(right, left)
        np.testing.assert_allclose(np.exp2(forward[578]), a @ b, rtol=2e-15)
        np.testing.assert_allclose(np.exp2(backward[578]), b @ a, rtol=2e-15)
        self.assertFalse(np.array_equal(a @ b, b @ a))
        self.assertTrue(np.all(np.isneginf(forward[:578])))

    def test_entrywise_and_branch_limits(self):
        coefficient = np.zeros((9, 2, 2))
        powered = general.entrywise_power(coefficient, 2)
        expected = np.array(list(range(1, 10))+list(range(8, 0, -1)))
        for i in range(2):
            for j in range(2):
                np.testing.assert_allclose(np.exp2(powered[:, i, j]), expected)
        with self.assertRaises(ValueError):
            general.branch_polynomial(np.zeros((1, 2, 2)), [0.], np.array([9]), 1)
        with self.assertRaises(ValueError):
            general.matrix_polynomial_product(np.zeros((1153, 2, 2)), np.zeros((2, 2, 2)))
        self.assertFalse(general.metadata()['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
