from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest import TestCase, main

import numpy as np

import hd_two_stream_mean_gate as gate
from hd_two_sum_gf16_mean_column import column_polynomial


def matmul(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


class TwoStreamMeanTests(TestCase):
    def test_common_potential_bounds_noncommuting_words(self):
        matrices = [
            [[F(1, 2), F(1, 10), 0, F(1, 30)],
             [F(1, 100), F(1, 3), F(1, 7), 0],
             [0, F(1, 7), F(1, 3), F(1, 11)],
             [F(1, 91), 0, F(1, 13), F(1, 4)]],
            [[F(2, 3), 0, F(1, 9), 0],
             [F(1, 97), F(1, 8), 0, F(1, 9)],
             [F(1, 101), F(1, 6), F(1, 5), 0],
             [0, F(1, 67), F(1, 43), F(1, 7)]]]
        self.assertNotEqual(matmul(*matrices), matmul(*reversed(matrices)))
        potential = [F(1), F(1, 4), F(1, 4), F(1, 16)]
        caps = [max(sum(row[j]*potential[j] for j in range(4))/potential[i]
                    for i, row in enumerate(matrix)) for matrix in matrices]
        for indices in product(range(2), repeat=4):
            row, bound = [[F(1), 0, 0, 0]], F(16)
            for index in indices:
                row = matmul(row, matrices[index])
                bound *= caps[index]
            self.assertLessEqual(sum(row[0]), bound)

    def test_potential_shape_and_once_only_endpoint(self):
        np.testing.assert_equal(gate.product_potential(1, -3), [0, -3])
        np.testing.assert_equal(gate.product_potential(2, -3), [0, -3, -3, -6])
        np.testing.assert_equal(gate.product_potential(2, -3, lumped=True), [0, -3, -6])

    def test_max_over_histograms_is_after_row_sum(self):
        with np.errstate(divide='ignore'):
            identity = np.log2(np.eye(2))
            left = np.log2([[.8, .1], [.2, .3]])
            right = np.log2([[.1, .8], [.2, .3]])
        bank = {(0, 0): identity, (1, 0): left, (2, 0): left, (2, 1): right}
        curve, contacts = gate.potential_curve(bank, 2, np.zeros(2))
        self.assertAlmostEqual(curve[2], log2(.9))
        self.assertEqual(contacts[2], [dict(full_pairs=0, initial_state=0),
                                      dict(full_pairs=1, initial_state=0)])
        self.assertLess(curve[2], np.logaddexp2.reduce(np.maximum(left, right)[0]))

    def test_probability_and_finite_mass_secants(self):
        for intercept, slope in ((-3., -2.), (2., -.5)):
            curve = intercept+slope*np.arange(19)
            curve[0] = 0.
            measure = gate.secant_measure(curve)
            self.assertAlmostEqual(measure['log2_atom'], slope)
            self.assertAlmostEqual(measure['log2_nonzero_weight'], intercept)
            self.assertEqual(measure['log2_total_mass'], max(0., intercept))
            self.assertLessEqual(measure['log2_probability'], 0.)
            self.assertEqual(measure['zero_moment_relaxed'], intercept > 0)

    def test_finite_measure_column_contraction_matches_labelled_masks(self):
        # Total massA=3, with normalized atoms0 andB=2 equally likely.
        mass, probability, atom = F(3), F(1, 2), F(2)
        weights, groups, regions = [F(1), F(2), F(3)], 2, 2
        exact = [F(0)]*3
        for latent in product((0, atom), repeat=regions):
            latent_mass = mass**regions*probability**regions
            for masks in product(range(4), repeat=groups):
                term, degree = latent_mass, 0
                for mask in masks:
                    weight = mask.bit_count()
                    degree += min(weight, 1)
                    term *= weights[weight]
                    for i in range(regions):
                        if mask >> i & 1:
                            term *= latent[i]
                exact[degree] += term
        actual = column_polynomial(np.log2(list(map(float, weights))), groups,
            log2(probability), log2(atom), arity=1)+regions*log2(mass)
        np.testing.assert_allclose(actual, list(map(log2, exact)), atol=2e-14)

    def test_exact_symmetric_quotient_preserves_ordered_products(self):
        matrices = [np.array([[3, 2, 2, 1], [4, 5, 6, 7],
                              [4, 6, 5, 7], [8, 9, 9, 2]], dtype=object),
                    np.array([[1, 3, 3, 2], [2, 4, 1, 5],
                              [2, 1, 4, 5], [6, 7, 7, 8]], dtype=object)]
        def lump(matrix):
            return np.array([[matrix[i, 0], matrix[i, 1]+matrix[i, 2], matrix[i, 3]]
                             for i in (0, 1, 3)], dtype=object)
        self.assertFalse(np.array_equal(matrices[0]@matrices[1], matrices[1]@matrices[0]))
        for indices in product(range(2), repeat=4):
            full, reduced = np.eye(4, dtype=object), np.eye(3, dtype=object)
            for index in indices:
                full = full@matrices[index]
                reduced = reduced@lump(matrices[index])
            self.assertEqual(sum(full[0]), sum(reduced[0]))
        v4, v3 = np.array([1., .25, .25, .0625]), np.array([1., .25, .0625])
        for matrix in matrices:
            self.assertEqual(max(matrix@v4/v4), max(lump(matrix)@v3/v3))

    def test_fail_closed_guards(self):
        with self.assertRaises(ValueError):
            gate.product_potential(2, 1.)
        with self.assertRaises(ValueError):
            gate.secant_measure(np.zeros(10))
        with self.assertRaises(ValueError):
            gate.potential_curve({}, 2, np.zeros(4))
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        bank = {(m, f): np.zeros((2, 2)) for m in range(3) for f in range(min(1, m//2)+1)}
        with self.assertRaises(ArithmeticError):
            gate.potential_curve(bank, 2, np.zeros(2))


if __name__ == '__main__':
    main()
