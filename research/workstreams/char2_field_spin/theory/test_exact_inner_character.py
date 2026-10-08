"""Independent exhaustive checks for the joint local marked transition."""

from fractions import Fraction
import itertools
import math
import unittest
from unittest.mock import patch

import numpy as np

from exact_inner_character import (aggregate, build_histogram, conditional_transfers,
                                   conditional_weight_transfer, field_tables,
                                   moment, integer_enumerator, positive_aggregate)


def brute(q, state, width, v, z, weight_filter=None):
    product, _ = field_tables(q)
    a = []
    for point in range(width):
        row = [1]
        for j in range(1, state):
            row.append(int(product[row[-1], point]))
        a.append(row)
    result = [[Fraction(0), Fraction(0)] for _ in range(2)]
    for r in itertools.product(range(q), repeat=state):
        nonzero = int(any(r))
        factor = Fraction(1, q**state-1) if nonzero else Fraction(1)
        emitted = []
        for lane in a:
            y = 0
            for coefficient, value in zip(lane, r):
                y ^= int(product[coefficient, value])
            emitted.append(y)
        for x in itertools.product(range(q), repeat=width):
            weight = sum(value != 0 for value in x)
            if weight_filter is not None and weight != weight_filter:
                continue
            following = list(r)
            for j, value in enumerate(x):
                for k in range(state):
                    following[k] ^= int(product[a[j][k], value])
            output_weight = sum(value != mask for value, mask in zip(x, emitted))
            result[nonzero][int(any(following))] += factor*v**weight*z**output_weight
    if weight_filter is not None:
        mass = math.comb(width, weight_filter)*(q-1)**weight_filter
        result = [[entry/mass for entry in row] for row in result]
    return result


class CharacterMomentTests(unittest.TestCase):
    def test_positive_enumerator_agrees_exactly(self):
        for q, state, width in ((4, 2, 4), (16, 2, 16), (16, 3, 16)):
            hist = build_histogram(q, state, width)
            polys = integer_enumerator(q, state, width)
            self.assertTrue(all(c > 0 for row in polys for poly in row for c in poly.values()))
            for v, z in ((Fraction(1, 11), Fraction(7, 13)), (Fraction(2), Fraction(1, 9))):
                self.assertEqual(positive_aggregate(v, z, state, q, width), moment(v, z, hist))
            for w in range(width+1):
                for row in (0, 1):
                    count = sum(c for poly in polys[row] for (i, _), c in poly.items() if i == w)
                    expected = math.comb(width, w)*(q-1)**w*(q**state-1 if row else 1)
                    self.assertEqual(count, expected)

    def test_gf16_three_state_coefficient_normalizations(self):
        q, state, width = 16, 3, 16
        polys = integer_enumerator(q, state, width)
        # Since C has full row rank, every syndrome has q^(width-state)
        # preimages.  Returned row1 counts one syndrome for each R!=0.
        self.assertEqual(sum(polys[0][0].values()), q**(width-state))
        self.assertEqual(sum(polys[1][0].values()), (q**state-1)*q**(width-state))
        self.assertEqual(polys[0][0].get((0, 0)), 1)
        self.assertFalse(any(w == 0 for w, _ in polys[0][1]))
        self.assertFalse(any(w == 0 for w, _ in polys[1][0]))
        # Summing over every input also gives every output vector for each R.
        for row in (0, 1):
            for h in range(width+1):
                count = sum(c for poly in polys[row] for (_, j), c in poly.items() if j == h)
                self.assertEqual(count, math.comb(width, h)*(q-1)**h*(q**state-1 if row else 1))

    def test_conditional_weights_against_exhaustive_inputs(self):
        z = Fraction(3, 7)
        for state in (1, 2):
            transfers = conditional_transfers(z, state, 4, 4)
            for w, transfer in enumerate(transfers):
                with self.subTest(state=state, weight=w):
                    self.assertEqual(transfer, brute(4, state, 4, Fraction(1), z, w))
                    self.assertEqual(transfer, conditional_weight_transfer(w, z, state, 4, 4))

    def test_conditional_weights_reconstruct_gf16_aggregate(self):
        q, state, width = 16, 3, 16
        v, z = Fraction(2, 9), Fraction(4, 7)
        transfers = conditional_transfers(z, state, q, width)
        self.assertTrue(all(isinstance(entry, Fraction)
                            for transfer in transfers for row in transfer for entry in row))
        rebuilt = [[sum(math.comb(width, w)*(q-1)**w*v**w*transfer[row][column]
                        for w, transfer in enumerate(transfers))
                    for column in (0, 1)] for row in (0, 1)]
        self.assertEqual(rebuilt, positive_aggregate(v, z, state, q, width))
        for transfer in conditional_transfers(Fraction(1), state, q, width):
            self.assertEqual([sum(row) for row in transfer], [1, 1])

    def test_aggregate_uses_positive_path(self):
        expected = np.asarray(positive_aggregate(Fraction(1, 10), Fraction(1, 3)), dtype=float)
        with patch('exact_inner_character.moment', side_effect=AssertionError('signed path used')):
            np.testing.assert_allclose(aggregate(.1, 1/3), expected, rtol=2e-15, atol=0)

    def test_positive_float_and_conditional_transfers(self):
        for v, z in ((Fraction(1, 10**9), Fraction(99, 100)),
                     (Fraction(1, 5), Fraction(1, 10)),
                     (Fraction(3, 2), Fraction(1, 2))):
            exact = np.asarray(positive_aggregate(v, z), dtype=float)
            floating = np.asarray(positive_aggregate(float(v), float(z)))
            np.testing.assert_allclose(floating, exact, rtol=3e-15, atol=0)
        exact = np.asarray(conditional_transfers(Fraction(1, 7)), dtype=float)
        floating = np.asarray(conditional_transfers(1/7))
        np.testing.assert_allclose(floating, exact, rtol=3e-15, atol=0)

    def test_exact_gf4_against_all_states_and_inputs(self):
        for state in (1, 2):
            hist = build_histogram(4, state, 4)
            for v, z in ((Fraction(1, 7), Fraction(2, 3)),
                         (Fraction(3, 2), Fraction(1, 5)),
                         (Fraction(1), Fraction(1))):
                self.assertEqual(moment(v, z, hist), brute(4, state, 4, v, z))

    def test_gf16_short_direct_comparison(self):
        hist = build_histogram(16, 1, 2)
        v, z = Fraction(2, 9), Fraction(3, 7)
        self.assertEqual(moment(v, z, hist), brute(16, 1, 2, v, z))

    def test_uniform_input_preserves_uniform_output(self):
        for q, state, width in ((4, 2, 4), (16, 1, 16), (16, 2, 16)):
            hist = build_histogram(q, state, width)
            z = Fraction(2, 7)
            matrix = moment(Fraction(1), z, hist)
            expected = (1+(q-1)*z)**width
            self.assertEqual(sum(matrix[0]), expected)
            self.assertEqual(sum(matrix[1]), expected)

    def test_zero_input_and_unit_mark(self):
        hist = build_histogram(16, 2, 16)
        matrix = moment(Fraction(0), Fraction(1, 3), hist)
        self.assertEqual(matrix[0], [1, 0])
        self.assertEqual(matrix[1][0], 0)
        for v in (Fraction(1, 9), Fraction(2)):
            matrix = moment(v, Fraction(1), hist)
            mass = (1+15*v)**16
            self.assertEqual(sum(matrix[0]), mass)
            self.assertEqual(sum(matrix[1]), mass)

    def test_float_agrees_with_exact_and_preserves_row_sums(self):
        hist = build_histogram(16, 2, 16)
        for v, z in ((Fraction(1, 1000), Fraction(99, 100)),
                     (Fraction(1, 5), Fraction(1, 10)),
                     (Fraction(3, 2), Fraction(1, 2))):
            exact = np.array(moment(v, z, hist), dtype=float)
            floating = np.array(moment(float(v), float(z), hist))
            np.testing.assert_allclose(floating, exact, rtol=2e-11, atol=1e-18)

    def test_field_and_geometry_guards(self):
        for args in ((8, 1, 8), (16, 0, 16), (16, 4, 16), (4, 1, 5), (16, 2.5, 16),
                     (np.int64(16), np.int64(16), np.int64(16))):
            with self.assertRaises(ValueError):
                build_histogram(*args)
        hist = build_histogram(4, 1, 4)
        for v, z in ((-1, .5), (1, 0), (1, 2), (math.nan, .5),
                     (math.inf, .5), (1, math.nan), (1, math.inf)):
            with self.assertRaises(ValueError):
                moment(v, z, hist)
            with self.assertRaises(ValueError):
                positive_aggregate(v, z, 1, 4, 4)
        for z in (math.nan, math.inf, 0, 2):
            with self.assertRaises(ValueError):
                conditional_transfers(z, 1, 4, 4)
        for w in (-1, 5, 1.5):
            with self.assertRaises(ValueError):
                conditional_weight_transfer(w, .5, 1, 4, 4)
        with self.assertRaises(ArithmeticError):
            positive_aggregate(1e100, .5, 1, 4, 4)


if __name__ == '__main__':
    unittest.main()
