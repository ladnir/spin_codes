from fractions import Fraction as F
from itertools import product
from math import log2
from unittest import TestCase, main

import numpy as np

import hd_two_sum_matrix_moment as gate


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def reference_column(weights, groups, arity, p, b0, bx, a01, a10, a11):
    size = len(weights)-1
    result = [[[F(0) for _ in range(groups*arity+1)] for _ in range(2)] for _ in range(2)]
    for masks in product(range(1 << size), repeat=groups):
        value = [[F(1), F(0)], [F(0), F(1)]]
        mass, degree = F(1), 0
        for mask in masks:
            mass *= weights[mask.bit_count()]
            degree += min(mask.bit_count(), arity)
        for region in range(size):
            m = sum(mask >> region & 1 for mask in masks)
            regional = [[(1-p if m == 0 else 0)+p*b0**m, a01*bx**m],
                        [a10*bx**m, a11]]
            value = multiply(value, regional)
        for i, j in product(range(2), repeat=2):
            result[i][j][degree] += mass*value[i][j]
    return result


class MatrixMomentTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2/3', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom_cross',
            'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))
        self.weights = tuple(map(F, (1, 2, 3, 5)))

    def test_all_labelled_masks_and_noncommuting_paths(self):
        exact = reference_column(self.weights, 2, 2, *self.constants)
        actual, cost = gate.column_matrix(np.log2(list(map(float, self.weights))), 2,
                                         self.measure, arity=2)
        np.testing.assert_allclose(2**actual, np.asarray(exact, dtype=float), rtol=3e-14)
        self.assertGreater(cost['positive_path_atom_terms'], 0)
        p, b0, bx, a01, a10, a11 = self.constants
        first = [[F(1), a01], [a10, a11]]
        second = [[p*b0, a01*bx], [a10*bx, a11]]
        self.assertNotEqual(multiply(first, second), multiply(second, first))

    def test_scalar_evaluation_is_the_same_positive_polynomial(self):
        logs = np.log2(list(map(float, self.weights)))
        polynomial, _ = gate.column_matrix(logs, 2, self.measure, arity=2)
        for marker in (0., 1.5):
            scalar, _ = gate.column_matrix(logs, 2, self.measure, arity=2, log_marker=marker)
            expected = np.logaddexp2.reduce(polynomial+marker*np.arange(5), axis=2)
            np.testing.assert_allclose(scalar, expected, atol=3e-14)

    def test_probability_one_zero_atom_boundary(self):
        constants = (F(1), *self.constants[1:])
        measure = dict(self.measure, log2_probability=0.)
        exact = reference_column(self.weights, 1, 2, *constants)
        actual, _ = gate.column_matrix(np.log2(list(map(float, self.weights))), 1,
                                      measure, arity=2)
        np.testing.assert_allclose(2**actual, np.asarray(exact, dtype=float), rtol=2e-14)

    def test_matrix_polynomial_product_and_power_keep_order(self):
        left = np.array([[[1., 2.], [3., 1.]], [[2., 4.], [1., 5.]]])
        right = np.array([[[4., 1.], [2., 3.]], [[1., 6.], [2., 1.]]])
        expected = np.zeros((2, 2, 3))
        for d, e in product(range(2), repeat=2):
            expected[:, :, d+e] += left[:, :, d]@right[:, :, e]
        actual = gate.matrix_poly_product(np.log2(left), np.log2(right))
        np.testing.assert_allclose(2**actual, expected, rtol=2e-14)
        self.assertFalse(np.array_equal(left[:, :, 0]@right[:, :, 0],
                                       right[:, :, 0]@left[:, :, 0]))
        squared = gate.matrix_poly_product(np.log2(left), np.log2(left))
        cubed = gate.matrix_poly_product(squared, np.log2(left))
        np.testing.assert_allclose(gate.matrix_poly_power(np.log2(left), 3), cubed, atol=2e-14)

    def test_frozen_entry_slope_covers_every_occupancy(self):
        curves = np.array([[0., -6., -8., -10.], [-np.inf, -2., -5., -4.],
                           [-np.inf, -9., -10., -11.], [-4., -3., -5., -4.]])
        measure = gate.entry_measure(curves, -4., -2., cross_endpoints=(1, 3))
        self.assertEqual(measure['log2_atom_cross'], -1.)
        self.assertEqual(measure['log2_scale01'], -1.)
        self.assertEqual(measure['log2_scale10'], -8.)
        self.assertEqual(measure['log2_scale11'], -3.)
        self.assertEqual(measure['contacts'][0], [0, 1, 2, 3])

    def test_chernoff_marker_and_direction(self):
        word = np.log2([1., 3., 2., 1.])
        marker = gate.adjacent_floor_marker(word, 2)
        self.assertEqual(marker, 1.)
        tail = np.logaddexp2.reduce(word[2:])
        bound = np.logaddexp2.reduce(word+marker*np.arange(4))-2*marker
        self.assertGreaterEqual(bound, tail)
        self.assertEqual(gate.adjacent_floor_marker(np.array([0., 1., 2.]), 1), 0.)

    def test_guards(self):
        with self.assertRaises(ValueError):
            gate.run(mode='adaptive')
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        with self.assertRaises(ValueError):
            gate.adjacent_floor_marker([0.], 0)
        with self.assertRaises(ValueError):
            gate.column_matrix([0., 1.], 2, self.measure, arity=1, log_marker=-1.)
        with self.assertRaises(ValueError):
            gate.matrix_poly_power(np.zeros((2, 2, 2)), 0)


if __name__ == '__main__':
    main()
