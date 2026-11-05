"""Tiny exact coefficient and chronology controls; no provider is run."""
from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

import hd_pair_state_tail as tail


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, dtype=float))


def scalar_convolution(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i+j] += x*y
    return result


def scalar_power(values, exponent):
    answer = [1]
    for _ in range(exponent):
        answer = scalar_convolution(answer, values)
    return answer


def matrix_product(left, right):
    left, right = np.asarray(left, object), np.asarray(right, object)
    result = np.zeros((len(left)+len(right)-1, 2, 2), dtype=object)
    for i, j in product(range(len(left)), range(len(right))):
        result[i+j] += left[i] @ right[j]
    return result


class PairStateTailTests(unittest.TestCase):
    def test_branch_weight_grouping_counts_and_roots(self):
        roots = np.array([[[1, 2], [3, 4]], [[2, 1], [1, 3]], [[3, 2], [2, 1]]])
        counts, weights = np.array([2, 3, 5]), np.array([0, 2, 2])
        got = tail.branch_polynomial(logs(roots**3), logs(counts), weights, 3)
        expected = np.zeros((3, 2, 2))
        expected[0] = 2*roots[0]
        expected[2] = 3*roots[1]+5*roots[2]
        np.testing.assert_allclose(np.exp2(got), expected, rtol=2e-15)
        self.assertTrue(np.all(np.isneginf(got[1])))

    def test_entrywise_power_is_not_matrix_power(self):
        coefficients = np.array([[[1, 2], [3, 4]], [[2, 1], [0, 2]], [[1, 0], [2, 1]]])
        for exponent in (0, 1, 3):
            got = tail.entrywise_power(logs(coefficients), exponent)
            expected = np.zeros((2*exponent+1, 2, 2), dtype=object)
            for i, j in product(range(2), repeat=2):
                expected[:, i, j] = scalar_power(list(coefficients[:, i, j]), exponent)
            np.testing.assert_allclose(np.exp2(got), np.asarray(expected, float), rtol=5e-15)

    def test_ordered_noncommuting_matrix_polynomials(self):
        a = np.array([[[1, 2], [0, 1]], [[0, 1], [3, 0]]])
        b = np.array([[[2, 0], [1, 1]], [[1, 3], [0, 1]], [[2, 1], [1, 0]]])
        for left, right in ((a, b), (b, a)):
            got = tail.matrix_polynomial_product(logs(left), logs(right))
            expected = matrix_product(left, right)
            np.testing.assert_allclose(np.exp2(got), np.asarray(expected, float), rtol=4e-15)
        self.assertFalse(np.array_equal(matrix_product(a, b), matrix_product(b, a)))

    def test_coefficientwise_holder_against_exact_positive_measure(self):
        # Two atoms, three branches, and different positive boundary measures.
        f = [[F(1), F(1, 3), F(2, 3)], [F(1), F(3, 4), F(1, 4)]]
        masses = np.array([[[F(1, 2), F(1, 4)], [F(1, 3), F(1, 5)]],
                           [[F(1, 4), F(1, 3)], [F(1, 7), F(1, 2)]]], dtype=object)
        counts, weights = [1, 2, 3], np.array([0, 1, 2])
        for groups in (1, 3):
            moments = np.zeros((3, 2, 2), dtype=object)
            exact = np.zeros((2*groups+1, 2, 2), dtype=object)
            for atom in range(2):
                for branch in range(3):
                    moments[branch] += f[atom][branch]**groups*masses[atom]
                coefficients = scalar_power([counts[j]*f[atom][j] for j in range(3)], groups)
                for degree, coefficient in enumerate(coefficients):
                    exact[degree] += coefficient*masses[atom]
            grouped = tail.branch_polynomial(logs(moments), logs(counts), weights, groups)
            upper = np.exp2(tail.entrywise_power(grouped, groups))
            exact = np.asarray(exact, float)
            if groups == 1:
                np.testing.assert_allclose(upper, exact, rtol=4e-15)
            self.assertTrue(np.all(upper >= exact*(1-2e-14)))

    def test_word_and_tail_against_exact_two_block_polynomial(self):
        counts, weights = [2, 3], np.array([0, 1])
        branch_roots = [np.array([[[1, 2], [1, 1]], [[2, 1], [1, 3]]]),
                        np.array([[[2, 1], [3, 1]], [[1, 2], [2, 1]]])]
        exact_blocks = []
        for roots in branch_roots:
            current = np.zeros((3, 2, 2), dtype=object)
            for i, j in product(range(2), repeat=2):
                current[:, i, j] = scalar_power([counts[k]*roots[k, i, j] for k in range(2)], 2)
            exact_blocks.append(current)
        expected = matrix_product(*exact_blocks)
        blocks = [logs(roots**2) for roots in branch_roots]
        got = tail.word_polynomial(blocks, logs(counts), weights, 2)
        np.testing.assert_allclose(np.exp2(got), np.asarray(expected, float), rtol=8e-15)
        for floor in (0, 1, 3, 4, 5):
            value = tail.word_tail(blocks, logs(counts), weights, 2, floor)
            exact = sum(expected[floor:, 0, :].flat)
            self.assertAlmostEqual(float(np.exp2(value)), float(exact), delta=float(exact)*1e-14)

    def test_batched_fixed_labels_match_separate_calls(self):
        blocks = [logs(np.array([[[[1, 2], [3, 1]], [[2, 1], [1, 2]]],
                                [[[3, 1], [2, 2]], [[1, 3], [2, 1]]]])),
                  logs(np.array([[[[2, 1], [1, 3]], [[1, 2], [2, 1]]],
                                [[[1, 3], [1, 2]], [[2, 1], [3, 1]]]]))]
        counts, weights = logs([[1, 2], [3, 1]]), np.array([[0, 1], [0, 2]])
        got = tail.word_tail(blocks, counts, weights, 2, 3)
        expected = [tail.word_tail([block[i] for block in blocks], counts[i], weights[i], 2, 3)
                    for i in range(2)]
        np.testing.assert_allclose(got, expected, rtol=2e-15)

    def test_zero_word_removed_without_subtraction_or_underflow(self):
        # The nonzero coefficient is far below binary64's ordinary range.
        matrices = np.full((2, 2, 2), -np.inf)
        matrices[0, 0, 0] = 0.
        matrices[1, 0, 0] = -18000.
        result = tail.word_tail([matrices], [0., 0.], np.array([0, 1]), 18, 18)
        self.assertTrue(np.isfinite(result))
        self.assertAlmostEqual(float(result), -18000., places=10)
        empty = tail.word_tail([matrices[:1]], [0.], np.array([0]), 18, 1)
        self.assertTrue(np.isneginf(empty))

    def test_time_guard_and_validation(self):
        def stop():
            raise TimeoutError('tiny test')
        with self.assertRaises(TimeoutError):
            tail.word_tail([np.zeros((1, 2, 2))], [0.], np.array([0]), 1, 0, check_time=stop)
        for weights, groups in ((np.array([5]), 1), (np.array([1.5]), 1), (np.array([0]), 19)):
            with self.assertRaises(ValueError):
                tail.branch_polynomial(np.zeros((1, 2, 2)), [0.], weights, groups)
        with self.assertRaises(ValueError):
            tail.matrix_polynomial_product(np.zeros((577, 2, 2)), np.zeros((2, 2, 2)))
        with self.assertRaises(ValueError):
            tail.tail_from_polynomial(np.zeros((1, 2, 2)), -1)
        self.assertFalse(tail.metadata()['whole_code_certificate'])
        self.assertFalse(tail.metadata()['original_column_activity_floor_included'])
        self.assertFalse(tail.metadata()['rowset_root_sum_included'])


if __name__ == '__main__':
    unittest.main()
