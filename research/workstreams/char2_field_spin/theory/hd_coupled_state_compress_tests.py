from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main

import numpy as np

import hd_coupled_state_compress as compression


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def exact_compress(matrix, weights):
    return [[matrix[0][0], max(matrix[0][j+1]/weights[j] for j in range(4))],
        [sum(weights[i]*matrix[i+1][0] for i in range(4)),
         max(sum(weights[i]*matrix[i+1][j+1] for i in range(4))/weights[j]
             for j in range(4))]]


class CoupledCompressionTests(TestCase):
    def setUp(self):
        self.weights = [F(1, 15), F(2, 15), F(4, 15), F(8, 15)]
        self.matrices = [
            [[F((3*i+5*j+1)%11, 13) for j in range(5)] for i in range(5)],
            [[F((7*i+2*j+3)%13, 17) for j in range(5)] for i in range(5)]]

    def test_exact_coefficient_embedding_domination(self):
        embedding = [[F(1), 0, 0, 0, 0], [0]+self.weights]
        for matrix in self.matrices:
            reduced = exact_compress(matrix, self.weights)
            left, right = multiply(embedding, matrix), multiply(reduced, embedding)
            self.assertTrue(all(a <= b for la, lb in zip(left, right) for a, b in zip(la, lb)))
            self.assertEqual(left[0][0], right[0][0])
            self.assertEqual(left[1][0], right[1][0])

    def test_noncommuting_chronological_words_terminal_once(self):
        self.assertNotEqual(multiply(*self.matrices), multiply(*reversed(self.matrices)))
        reduced = [exact_compress(matrix, self.weights) for matrix in self.matrices]
        for word in product(range(2), repeat=4):
            full, small = [[F(1), 0, 0, 0, 0]], [[F(1), 0]]
            for index in word:
                full = multiply(full, self.matrices[index])
                small = multiply(small, reduced[index])
            self.assertLessEqual(sum(full[0]), sum(small[0]))

    def test_float_compression_matches_exact_and_preserves_zero(self):
        weights = np.log2(list(map(float, self.weights)))
        for matrix in self.matrices:
            with np.errstate(divide='ignore'):
                actual = compression.compress_matrix(np.log2(np.array(matrix, dtype=float)), weights)
            expected = np.array(exact_compress(matrix, self.weights), dtype=float)
            np.testing.assert_allclose(np.exp2(actual), expected, atol=1e-14)
        with np.errstate(divide='ignore'):
            zero = compression.compress_matrix(np.log2(np.eye(5)), weights)
        np.testing.assert_allclose(np.exp2(zero), np.eye(2), atol=1e-15)

    def test_normalization_is_common_and_scale_invariant(self):
        one, total = compression.normalize_live_weights((-3., -2., 0., -7.))
        two, total2 = compression.normalize_live_weights((17., 18., 20., 13.))
        np.testing.assert_allclose(one, two, atol=2e-15)
        self.assertAlmostEqual(total2-total, 20.)
        self.assertAlmostEqual(np.exp2(one).sum(), 1.)
        with np.errstate(divide='ignore'):
            matrix = np.log2(np.array(self.matrices[0], dtype=float))
        bank = {(m, f): matrix for m in range(3) for f in range(min(1, m//2)+1)}
        result = compression.compress_bank(bank, maximum=2)
        self.assertEqual(result['metadata']['terminal_vector'], [1, 1])
        self.assertEqual(result['metadata']['extra_terminal_factor_log2'], 0.)
        curves = compression.entrywise_curves(result['regions'], maximum=2)
        self.assertEqual(curves.shape, (4, 3))

    def test_arbitrary_point_family_remains_closed_under_mixtures(self):
        #P has two possible atoms. Positive output contributions produce an
        #arbitrary convex mixture, not a uniform reference or a density factor.
        weights = self.weights
        first = [F(1), F(0)]
        second = [F(0), F(1)]
        coefficients = [F(1, 7), F(2, 9)]
        total = sum(coefficients)
        point_measure = [(coefficients[0]*first[i]+coefficients[1]*second[i])/total
                         for i in range(2)]
        self.assertEqual(sum(point_measure), 1)
        self.assertNotEqual(point_measure, [F(1, 2)]*2)
        cap = total/weights[3]
        self.assertEqual([cap*weights[3]*p for p in point_measure], coefficients)

    def test_fail_closed_guards(self):
        for values in ((0., 0., 0.), (0., 0., 0., -np.inf)):
            with self.assertRaises(ValueError):
                compression.normalize_live_weights(values)
        with self.assertRaises(ValueError):
            compression.compress_matrix(np.zeros((5, 5)), np.zeros(4))
        with self.assertRaises(ValueError):
            compression.compress_bank({}, maximum=2)


if __name__ == '__main__':
    main()
