"""Exact toy checks for positive nonempty-column transfer coefficients."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

import numpy as np

from hd_xor_dense_column_floor import column_polynomial, matrix_product, split_column, tail_log_mass


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


class DenseColumnFloorTests(unittest.TestCase):
    def test_positive_flag_split_matches_all_labelled_occupancy_patterns(self):
        matrices = [
            [[[1, 0], [0, F(1, 4)]], [[F(1, 8), F(1, 2)], [F(1, 16), F(1, 4)]]],
            [[[F(1, 2), F(1, 4)], [0, F(1, 8)]], [[F(1, 4), F(1, 8)], [F(1, 2), F(1, 16)]]],
            [[[F(1, 4), 0], [F(1, 8), F(1, 2)]], [[F(1, 16), F(1, 8)], [F(1, 4), F(1, 2)]]],
        ]
        eta = F(1, 4)
        exact = [np.zeros((2, 2), dtype=object), np.zeros((2, 2), dtype=object)]
        for pattern in product(range(2), repeat=3):
            value = [[F(1), F(0)], [F(0), F(1)]]
            for position, occupied in enumerate(pattern):
                value = multiply(value, matrices[position][occupied])
            probability = eta**sum(pattern)*(1-eta)**(3-sum(pattern))
            exact[int(any(pattern))] += probability*np.array(value, dtype=object)
        with np.errstate(divide='ignore'):
            actual = split_column(np.log2(np.array(matrices, dtype=float)), float(eta))
        for index in range(2):
            np.testing.assert_allclose(np.exp2(actual[index]), np.array(exact[index], dtype=float), rtol=1e-14)

    def test_noncommuting_column_coefficients_match_word_enumeration(self):
        empty = [[F(1, 2), F(1, 4)], [F(0), F(1, 8)]]
        occupied = [[F(1, 16), F(1, 8)], [F(1, 4), F(1, 2)]]
        exact = np.zeros((5, 2, 2), dtype=object)
        for word in product(range(2), repeat=4):
            value = [[F(1), F(0)], [F(0), F(1)]]
            for nonempty in word:
                value = multiply(value, occupied if nonempty else empty)
            exact[sum(word)] += np.array(value, dtype=object)
        with np.errstate(divide='ignore'):
            actual = column_polynomial(np.log2(np.array(empty, dtype=float)),
                                       np.log2(np.array(occupied, dtype=float)), 4)
        np.testing.assert_allclose(np.exp2(actual), np.array(exact, dtype=float), rtol=1e-14)
        expected_tail = sum(exact[degree, 0, column] for degree in range(2, 5) for column in range(2))
        self.assertAlmostEqual(2**tail_log_mass(actual, 2), float(expected_tail), places=15)

    def test_boolean_column_mark_is_not_a_nonzero_entry_counter(self):
        with np.errstate(divide='ignore'):
            identity = np.log2(np.eye(2))
        regions = np.array([[identity, identity]]*3)
        empty, occupied = split_column(regions, .5)
        np.testing.assert_allclose(np.exp2(empty), np.eye(2)/8)
        np.testing.assert_allclose(np.exp2(occupied), np.eye(2)*7/8)
        # Three occupied positions still contribute degree one, not three.
        polynomial = column_polynomial(empty, occupied, 3)
        for degree in range(4):
            expected = comb(3, degree)*(7/8)**degree*(1/8)**(3-degree)
            self.assertAlmostEqual(2**polynomial[degree, 0, 0], expected, places=14)

    def test_zero_columns_are_removed_by_the_tail_not_subtracted(self):
        with np.errstate(divide='ignore'):
            identity = np.log2(np.eye(2))
        empty, occupied = split_column(np.array([[identity, identity]]*8), 0.)
        polynomial = column_polynomial(empty, occupied, 32)
        self.assertEqual(tail_log_mass(polynomial, 0), 0.)
        self.assertEqual(tail_log_mass(polynomial, 17), -np.inf)
        empty, occupied = split_column(np.array([[identity, identity]]*8), 1.)
        polynomial = column_polynomial(empty, occupied, 32)
        self.assertEqual(tail_log_mass(polynomial, 17), 0.)
        self.assertTrue(np.all(np.isneginf(polynomial[:32])))

    def test_toy_mds_words_and_invertible_outer_preserve_nonempty_column_floor(self):
        def gf4_multiply(a, b):
            value = 0
            while b:
                if b & 1:
                    value ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return value
        words = [(a, a ^ b, a ^ gf4_multiply(2, b))
                 for a, b in product(range(4), repeat=2) if a or b]
        minimum = 3
        for rows in product(words, repeat=3):
            occupied = 0
            for column in range(3):
                original = [rows[row][column] for row in range(3)]+[0]
                scaled = [gf4_multiply(value, 1+(column+row) % 3) for row, value in enumerate(original)]
                total = scaled[0] ^ scaled[1] ^ scaled[2] ^ scaled[3]
                output = [value ^ total for value in scaled]
                self.assertEqual(any(original), any(output))
                occupied += bool(any(output))
            self.assertGreaterEqual(occupied, 2)
            minimum = min(minimum, occupied)
        self.assertEqual(minimum, 2)

    def test_invalid_shapes_and_floor_are_rejected(self):
        with self.assertRaises(ValueError):
            split_column(np.zeros((2, 2, 2)), .5)
        with self.assertRaises(ValueError):
            split_column(np.zeros((8, 2, 2, 2)), 1.01)
        with self.assertRaises(ValueError):
            tail_log_mass(np.zeros((33, 2, 2)), 33)


if __name__ == '__main__':
    unittest.main()
