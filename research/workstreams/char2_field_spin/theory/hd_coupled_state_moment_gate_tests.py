from fractions import Fraction as F
from itertools import product
from math import log2
from unittest import TestCase, main

import numpy as np

import hd_coupled_state_moment_gate as gate
import hd_two_sum_matrix_moment as old


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def reference_column(weights, groups, arity, p, b0, b01, b10, a01, a10, a11):
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
            regional = [[(1-p if m == 0 else 0)+p*b0**m, a01*b01**m],
                        [a10*b10**m, a11]]
            value = multiply(value, regional)
        for i, j in product(range(2), repeat=2):
            result[i][j][degree] += mass*value[i][j]
    return result


class CoupledMomentTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2/3', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))
        self.weights = tuple(map(F, (1, 2, 3, 5)))

    def test_separate_cross_atoms_match_exact_labelled_noncommuting_words(self):
        expected = reference_column(self.weights, 2, 2, *self.constants)
        actual, costs = gate.column_matrix(np.log2(list(map(float, self.weights))), 2,
            self.measure, arity=2)
        np.testing.assert_allclose(2**actual, np.asarray(expected, dtype=float), rtol=4e-14)
        self.assertTrue(costs['separate_cross_counts_in_cache'])
        p, b0, b01, b10, a01, a10, a11 = self.constants
        first, second = [[F(1), a01], [a10, a11]], [[p*b0, a01*b01], [a10*b10, a11]]
        self.assertNotEqual(multiply(first, second), multiply(second, first))

    def test_equal_cross_slopes_recover_pinned_old_engine(self):
        measure = dict(self.measure, log2_atom_cross=self.measure['log2_atom01'])
        converted = gate.from_shared_cross(measure)
        weights = np.log2(list(map(float, self.weights)))
        expected, _ = old.column_matrix(weights, 2, measure, arity=2)
        actual, _ = gate.column_matrix(weights, 2, converted, arity=2)
        np.testing.assert_array_equal(actual, expected)

    def test_cache_must_distinguish_direction(self):
        weights = np.log2([1., 2., 3., 5.])
        left = gate.group_polynomial(weights, (1, 1, 0, 1), 1, self.measure, 2)
        right = gate.group_polynomial(weights, (1, 0, 1, 1), 1, self.measure, 2)
        self.assertFalse(np.allclose(left, right))

    def test_entry_secants_and_all_occupancy_intercepts(self):
        curves = np.array([[0., -6., -8., -10.], [-np.inf, -2., -2., -4.],
                           [-np.inf, -9., -10., -15.], [-4., -3., -5., -4.]])
        measure = gate.entry_measure(curves, zero_endpoints=(1, 3), cross_endpoints=(1, 3))
        self.assertEqual(measure['log2_atom00'], -2.)
        self.assertEqual(measure['log2_probability'], -4.)
        self.assertEqual(measure['log2_atom01'], -1.)
        self.assertEqual(measure['log2_atom10'], -3.)
        self.assertEqual(measure['log2_scale01'], 0.)
        self.assertEqual(measure['log2_scale10'], -4.)
        self.assertEqual(measure['log2_scale11'], -3.)
        self.assertEqual(measure['contacts'][1], [2])
        self.assertGreaterEqual(measure['minimum_slack'], -1e-14)

    def test_probability_one_and_terminal_vector(self):
        constants = (F(1), *self.constants[1:])
        measure = dict(self.measure, log2_probability=0.)
        expected = reference_column(self.weights, 1, 2, *constants)
        actual, _ = gate.column_matrix(np.log2(list(map(float, self.weights))), 1, measure, arity=2)
        np.testing.assert_allclose(2**actual, np.asarray(expected, dtype=float), rtol=3e-14)
        word = old.matrix_poly_power(actual, 2)
        # The normalized live cone ends against ones, without an extra potential.
        scalar = np.exp2(actual).sum(axis=2)
        self.assertAlmostEqual(float(np.exp2(word[0]).sum()), float((scalar@scalar)[0].sum()), places=10)

    def test_nonprobability_00_is_reported_not_clipped(self):
        curves = np.array([[0., -1., -1., -5.], [-np.inf, -2., -3., -4.],
                           [-np.inf, -3., -4., -5.], [-1., -2., -3., -4.]])
        with self.assertRaises(gate.MomentFitError):
            gate.entry_measure(curves, zero_endpoints=(1, 3), cross_endpoints=(1, 3))

    def test_guards(self):
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        with self.assertRaises(ValueError):
            gate.column_matrix([0., 1.], 1, dict(self.measure, log2_probability=1.), arity=1)
        with self.assertRaises(ValueError):
            gate.group_polynomial([0., 1.], (1, 0, 0, 0), 2, self.measure, 1)


if __name__ == '__main__':
    main()
