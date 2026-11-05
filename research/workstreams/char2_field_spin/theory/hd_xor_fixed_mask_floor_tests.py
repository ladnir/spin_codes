"""Independent finite-field law and RS measure normalization checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np

from hd_xor_column_probe import mds_prescribed_support, xor_column
from hd_xor_dense_column_floor import column_polynomial
from hd_xor_fixed_mask_floor import beta_log2, comparison_counts, fixed_column_split, mask_operators


def multiply4(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


class FixedMaskFloorTests(unittest.TestCase):
    def test_uniform_input_including_zero_and_output_diagonals_exactly(self):
        q, group = 4, 4
        for active_mask in range(1, 16):
            active = [i for i in range(group) if active_mask >> i & 1]
            actual = Counter()
            for values in product(range(q), repeat=len(active)):
                original = [0]*group
                for index, value in zip(active, values):
                    original[index] = value
                mixed = xor_column(original)
                for diagonal in product(range(1, q), repeat=group):
                    actual[tuple(multiply4(x, d) for x, d in zip(mixed, diagonal))] += 1
            counts = comparison_counts(active_mask, group, q)
            total = q**len(active)*(q-1)**group
            for target in product(range(q), repeat=group):
                mask = sum(bool(value) << i for i, value in enumerate(target))
                expected = F(counts[mask], q**len(active)*(q-1)**mask.bit_count())
                self.assertEqual(F(actual[target], total), expected)
            self.assertEqual(F(actual[(0,)*group], total), F(1, q**len(active)))
        self.assertEqual(comparison_counts(15, group, q), tuple(3**mask.bit_count() for mask in range(16)))

    def test_beta_dominates_actual_nonzero_rs_rows_on_every_toy_support(self):
        q, n, k = 4, 3, 2
        for active in range(1, 5):
            beta = F(q**(n*active), (q-1)**((n-k)*active))
            comparison_point_mass = beta/F(q**(n*active))
            self.assertEqual(comparison_point_mass, F(1, (q-1)**((n-k)*active)))
            self.assertAlmostEqual(beta_log2(active, n, k, q), log2(beta), places=12)
            for supports in product(range(1 << n), repeat=active):
                actual = F(1)
                for mask in supports:
                    weight = mask.bit_count()
                    # Active rows exclude the zero word, even though the
                    # uniform comparison assigns positive mass to that row.
                    actual *= F(mds_prescribed_support(q, n, k, weight), (q-1)**weight) if weight else 0
                self.assertLessEqual(actual, comparison_point_mass)

    def test_mask_table_and_fixed_split_preserve_labelled_chronology(self):
        tables = np.array([[[[1., 0.], [.25, .125]], [[.125, .5], [.25, .5]]],
                           [[[.5, .25], [.125, .5]], [[.25, .125], [.5, .25]]],
                           [[[.5, 0.], [.25, .5]], [[.125, .25], [.5, .125]]],
                           [[[.25, .5], [0., .5]], [[.5, .125], [.25, .125]]]])
        with np.errstate(divide='ignore'):
            operators = mask_operators(np.log2(tables))
        for mask in range(16):
            expected = np.eye(2)
            for position in range(4):
                expected = expected @ tables[position, (mask >> position) & 1]
            np.testing.assert_allclose(np.exp2(operators[mask]), expected, rtol=1e-14)
        counts = comparison_counts(7, 4, 4)
        empty, nonempty = fixed_column_split(operators, 7, 4, 4)
        expected = sum(counts[mask]*np.exp2(operators[mask])/64 for mask in range(1, 16))
        np.testing.assert_allclose(np.exp2(empty), np.exp2(operators[0])/64)
        np.testing.assert_allclose(np.exp2(nonempty), expected, rtol=1e-14)

    def test_fixed_row_sets_are_not_averaged_before_powering(self):
        # A deliberately labelled operator bank makes premature averaging
        # visible. Both sets have the same arity and beta, but different laws.
        operators = np.log2(np.array([[[1.+mask**3, .25], [.125, .5+mask**3/2]] for mask in range(16)]))
        splits = [fixed_column_split(operators, mask, 4, 4) for mask in (7, 14)]
        correct = sum(np.exp2(column_polynomial(*split, columns=3)) for split in splits)
        average = [np.logaddexp2(splits[0][i], splits[1][i])-1 for i in range(2)]
        resampled = 2*np.exp2(column_polynomial(*average, columns=3))
        self.assertGreater(float(np.max(np.abs(correct-resampled))), 1.)

    def test_identity_operators_have_exact_probability_normalization(self):
        with np.errstate(divide='ignore'):
            identity = np.log2(np.eye(2))
        operators = np.tile(identity, (256, 1, 1))
        for mask in (7, 31, 127, 255):
            empty, nonempty = fixed_column_split(operators, mask)
            np.testing.assert_allclose(np.exp2(np.logaddexp2(empty, nonempty)), np.eye(2), rtol=1e-14)
            self.assertAlmostEqual(2**empty[0, 0], 256**(-mask.bit_count()), places=30)


if __name__ == '__main__':
    unittest.main()
