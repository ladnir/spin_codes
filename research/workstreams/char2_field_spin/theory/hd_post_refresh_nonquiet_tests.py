from fractions import Fraction as Rat
from itertools import product
from math import log2
from unittest import TestCase, main

import numpy as np

import hd_post_refresh_nonquiet as flagged
from hd_coupled_state_moment_gate import column_matrix
from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_two_sum_matrix_moment import matrix_poly_power


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def reference_column(weights, groups, arity, p, b0, b01, b10, a01, a10, a11):
    size = len(weights)-1
    result = np.empty((2, 2, 2, groups*arity+1), dtype=object)
    result.fill(Rat(0))
    for masks in product(range(1 << size), repeat=groups):
        mass, degree, flag = Rat(1), 0, 0
        for mask in masks:
            weight = mask.bit_count()
            mass *= weights[weight]
            degree += min(weight, arity)
            flag |= weight > arity
        value = [[Rat(1), Rat(0)], [Rat(0), Rat(1)]]
        for region in range(size):
            m = sum(mask >> region & 1 for mask in masks)
            step = [[(1-p if m == 0 else 0)+p*b0**m, a01*b01**m],
                    [a10*b10**m, a11]]
            value = multiply(value, step)
        for i, j in product(range(2), repeat=2):
            result[flag, i, j, degree] += mass*value[i][j]
    return result


class PostRefreshNonquietTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(Rat, ('2/5', '1/3', '2/3', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))
        self.weights = tuple(map(Rat, (1, 2, 3, 5)))

    def test_scalar_OR_power_matches_complete_integer_word_enumeration(self):
        terms = ((0, 0, 2), (0, 1, 3), (1, 0, 5), (1, 2, 7))
        base = np.zeros((2, 3))
        for flag, degree, count in terms:
            base[flag, degree] += count
        expected = np.zeros((2, 7), dtype=int)
        for word in product(terms, repeat=3):
            flag, degree, count = 0, 0, 1
            for bit, weight, mass in word:
                flag |= bit
                degree += weight
                count *= mass
            expected[flag, degree] += count
        with np.errstate(divide='ignore'):
            actual = flagged.flag_power(np.log2(base), 3)
        np.testing.assert_allclose(np.exp2(actual), expected, rtol=5e-15, atol=1e-12)

    def test_exact_labelled_masks_and_ordered_state_products(self):
        expected = reference_column(self.weights, 2, 2, *self.constants)
        actual, costs = flagged.column_matrices(np.log2(list(map(float, self.weights))),
            2, self.measure, arity=2)
        np.testing.assert_allclose(np.exp2(actual), np.asarray(expected, dtype=float), rtol=5e-14)
        self.assertTrue(costs['separate_cross_counts_in_cache'])
        p, b0, b01, b10, a01, a10, a11 = self.constants
        first = [[Rat(1), a01], [a10, a11]]
        second = [[p*b0, a01*b01], [a10*b10, a11]]
        self.assertNotEqual(multiply(first, second), multiply(second, first))

    def test_flags_recombine_to_unflagged_column_and_whole_word(self):
        weights = np.log2(list(map(float, self.weights)))
        column, _ = flagged.column_matrices(weights, 2, self.measure, arity=2)
        expected, _ = column_matrix(weights, 2, self.measure, arity=2)
        np.testing.assert_allclose(np.logaddexp2(column[0], column[1]), expected, atol=1e-13)
        word = flagged.matrix_power(column, 3)
        ordinary = matrix_poly_power(expected, 3)
        np.testing.assert_allclose(np.logaddexp2(word[0], word[1]), ordinary, atol=1e-12)
        tails = flagged.tail_logs(word, 7)
        self.assertAlmostEqual(float(np.logaddexp2.reduce(tails)),
            float(np.logaddexp2.reduce(ordinary[0, :, 7:].ravel())), places=11)

    def test_quiet_flag_exactly_matches_quiet_only_engine(self):
        weights = np.log2(list(map(float, self.weights)))
        quiet = weights.copy()
        quiet[3:] = -np.inf
        column, _ = flagged.column_matrices(weights, 2, self.measure, arity=2)
        expected, _ = column_matrix(quiet, 2, self.measure, arity=2)
        np.testing.assert_allclose(column[0], expected, atol=1e-13)
        only, _ = flagged.column_matrices(quiet, 2, self.measure, arity=2)
        self.assertTrue(np.isneginf(only[1]).all())
        self.assertTrue(np.isneginf(flagged.matrix_power(only, 2)[1]).all())

    def test_matrix_OR_product_keeps_noncommuting_factor_order(self):
        left = np.array([[[[1., 2.], [3., 0.]], [[0., 1.], [1., 1.]]],
                         [[[0., 1.], [2., 0.]], [[1., 0.], [0., 3.]]]])
        right = np.array([[[[2., 0.], [0., 1.]], [[1., 3.], [2., 1.]]],
                          [[[1., 1.], [0., 2.]], [[3., 0.], [1., 0.]]]])
        expected = np.zeros((2, 2, 2, 3))
        for fa, fb, i, k, j, da, db in product(range(2), repeat=7):
            expected[fa | fb, i, j, da+db] += left[fa, i, k, da]*right[fb, k, j, db]
        with np.errstate(divide='ignore'):
            actual = flagged.matrix_product(np.log2(left), np.log2(right))
            reversed_product = flagged.matrix_product(np.log2(right), np.log2(left))
        np.testing.assert_allclose(np.exp2(actual), expected, rtol=2e-14, atol=1e-13)
        self.assertFalse(np.allclose(actual, reversed_product))

    def test_actual_half_arity_weight_split_and_zero_exponents(self):
        full, quiet = raw_holder_weights(), raw_holder_weights(quiet=True)
        np.testing.assert_array_equal(full[:9], quiet[:9])
        self.assertTrue(np.isneginf(quiet[9:]).all())
        self.assertTrue(np.isfinite(full[9:]).all())
        scalar = flagged.flag_power(np.zeros((2, 2)), 0)
        np.testing.assert_array_equal(scalar, [[0.], [-np.inf]])
        matrix = flagged.matrix_power(np.zeros((2, 2, 2, 2)), 0)
        self.assertTrue(np.isneginf(matrix[1]).all())
        np.testing.assert_array_equal(np.exp2(matrix[0, :, :, 0]), np.eye(2))

    def test_probability_one_and_guards(self):
        measure = dict(self.measure, log2_probability=0.)
        expected = reference_column(self.weights, 1, 2, Rat(1), *self.constants[1:])
        actual, _ = flagged.column_matrices(np.log2(list(map(float, self.weights))),
            1, measure, arity=2)
        np.testing.assert_allclose(np.exp2(actual), np.asarray(expected, dtype=float), rtol=4e-14)
        with self.assertRaises(ValueError):
            flagged.flag_power(np.zeros((2, 1)), -1)
        with self.assertRaises(ValueError):
            flagged.matrix_product(np.zeros((2, 2, 1)), np.zeros((2, 2, 1)))
        with self.assertRaises(ValueError):
            flagged.tail_logs(np.zeros((2, 2, 2, 1)), 1)
        with self.assertRaises(ValueError):
            flagged.column_matrices([0., 0.], 1, dict(self.measure, log2_probability=1.), arity=1)
        with self.assertRaises(ValueError):
            flagged.post_refresh_remainder(self.measure, np.nan)


if __name__ == '__main__':
    main()
