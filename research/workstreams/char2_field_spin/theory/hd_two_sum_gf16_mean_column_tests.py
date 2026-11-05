"""Small independent positive contractions, no regional bank/proof run."""
from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np

from hd_two_sum_closed_kernel import branch_counts, uniform_mask_count
from hd_two_sum_gf16_mean_column import (column_polynomial, common_matrix,
    raw_holder_weights, state_factor, two_atom_moments, whole_word_screen)


def multiply(left, right):
    out = [F(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] += a*b
    return out


def matrix_multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def exact_column(weights, groups, p, atom, arity):
    size = len(weights)-1
    out = [F(0)]*(arity*groups+1)
    for masks in product(range(1 << size), repeat=groups):
        counts = [sum(mask >> i & 1 for mask in masks) for i in range(size)]
        degree = sum(min(mask.bit_count(), arity) for mask in masks)
        mass = F(1)
        for mask in masks:
            mass *= weights[mask.bit_count()]
        for count in counts:
            mass *= p*atom**count if count else 1
        out[degree] += mass
    return out


class MeanColumnTests(unittest.TestCase):
    def test_latent_contraction_against_all_labelled_masks(self):
        weights = tuple(map(F, (2, 3, 5, 7)))
        for p in (F(2, 7), F(1)):
            expected = exact_column(weights, 2, p, F(3, 5), 2)
            actual = column_polynomial(np.log2(list(map(float, weights))), 2,
                log2(p), log2(F(3, 5)), arity=2)
            np.testing.assert_allclose(actual, [log2(x) for x in expected], atol=2e-13)

    def test_full_state_average_and_chronological_power(self):
        weights, p, atom = tuple(map(F, (1, 2, 4))), F(3, 8), F(2, 3)
        matrix = [[F(2), F(1, 5)], [F(1, 7), F(3, 4)]]
        base = exact_column(weights, 2, p, atom, 1)
        # Independently multiply R(m)=M_m*C for each labelled mask array.
        coefficients = [[[F(0) for _ in range(2)] for _ in range(2)] for _ in base]
        for masks in product(range(4), repeat=2):
            degree = sum(min(mask.bit_count(), 1) for mask in masks)
            product_matrix = [[F(1), F(0)], [F(0), F(1)]]
            mass = weights[masks[0].bit_count()]*weights[masks[1].bit_count()]
            for coordinate in range(2):
                count = sum(mask >> coordinate & 1 for mask in masks)
                moment = p*atom**count if count else F(1)
                regional = [[entry*moment for entry in row] for row in matrix]
                product_matrix = matrix_multiply(product_matrix, regional)
            for i, j in product(range(2), repeat=2):
                coefficients[degree][i][j] += mass*product_matrix[i][j]
        squared = matrix_multiply(matrix, matrix)
        for degree, value in enumerate(base):
            self.assertEqual(coefficients[degree], [[value*x for x in row] for row in squared])
        matrix_power = [[F(1), F(0)], [F(0), F(1)]]
        for _ in range(6):
            matrix_power = matrix_multiply(matrix_power, matrix)
        actual = state_factor(np.log2(np.asarray(matrix, dtype=float)), 6)
        self.assertAlmostEqual(actual, log2(sum(matrix_power[0])), places=13)

    def test_common_matrix_dominates_every_entry(self):
        moments = two_atom_moments(3, -2., -1.5)
        matrix = np.array([[0., -3.], [-4., -2.]])
        regions = moments[:, None, None]+matrix
        regions[1, 0, 1] = -np.inf
        got, reproduced = common_matrix(regions, -2., -1.5)
        np.testing.assert_array_equal(got, matrix)
        np.testing.assert_array_equal(reproduced, moments)
        self.assertTrue(np.all(regions <= moments[:, None, None]+got))

    def test_half_arity_branch_separation(self):
        for group in (8, 16):
            arity = group//2
            for inside in range(arity+1):
                quiet, hole, full = branch_counts(arity, inside, group)
                if hole:
                    self.assertGreaterEqual(group-arity-1+inside, arity+1)
                if full:
                    self.assertGreaterEqual(group-arity+inside, arity+1)
                if quiet:
                    self.assertLessEqual(inside, arity)
            quiet = raw_holder_weights(quiet=True, group=group, arity=arity)
            full = raw_holder_weights(group=group, arity=arity)
            np.testing.assert_array_equal(quiet[:arity+1], full[:arity+1])
            self.assertTrue(np.all(np.isneginf(quiet[arity+1:])))

    def test_fixed_rowset_holder_control(self):
        group, arity, columns = 8, 4, 2
        weights = raw_holder_weights(group=group, arity=arity, columns=columns)
        for left, right in ((0, 0), (15, 240), (31, 127), (85, 255)):
            exact = sum(uniform_mask_count(a, left, group)*uniform_mask_count(a, right, group)
                        for a in range(1 << group) if a.bit_count() == arity)
            upper_log = weights[left.bit_count()]+weights[right.bit_count()]
            if exact:
                self.assertLessEqual(log2(exact), upper_log+1e-12)

    def test_normalization_once_and_no_rowset_factor(self):
        weights = np.log2([1., 2., 3.])
        matrix = np.log2([[2., 1.], [1., 3.]])
        args = dict(log_source_count=7., mark=.9, limit=11, columns=3,
                    arity=1, dimension=1, minimum=2)
        actual = whole_word_screen(weights, matrix, 2, -1., -2., **args)
        raw = exact_column(tuple(map(F, (1, 2, 3))), 2, F(1, 2), F(1, 4), 1)
        word = multiply(multiply(raw, raw), raw)
        expected = (7-11*log2(.9)-3*log2(255)+log2(sum(word[4:]))
                    +state_factor(matrix, 6))
        self.assertAlmostEqual(actual['log2_screen'], expected, places=12)
        self.assertFalse(actual['extra_rowset_factor'])

    def test_guards_and_zero_state(self):
        with self.assertRaises(ValueError):
            two_atom_moments(3, 1., -1.)
        with self.assertRaises(ValueError):
            raw_holder_weights(group=16, arity=7)
        with self.assertRaises(ValueError):
            common_matrix(np.zeros((3, 2)), -1., -1.)
        self.assertEqual(state_factor(np.full((2, 2), -np.inf), 0), 0.)
        self.assertEqual(state_factor(np.full((2, 2), -np.inf), 1), -np.inf)


if __name__ == '__main__':
    unittest.main()
