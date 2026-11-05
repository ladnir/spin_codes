"""Exact toy checks for the offset DP and two independent support floors."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_geometry_strata import category_resources, pair_resource_float
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_floor import (offset_resource_vectors, mixed_mask_vectors,
    fixed_mask_branch_caps, column_log_coefficients, log_convolve, log_powers,
    bivariate_log_tail)


def multiply(left, right):
    return tuple(sum(left[2*i+k]*right[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


class MixedMaskFloorTests(unittest.TestCase):
    def test_all_labelled_offsets_match_separate_resource_dps(self):
        tables = np.array([[[[1+i+m, .25+m], [.125+i, .5+m]] for m in range(4)]
                           for i in range(4)])
        logs, logp = np.log2(tables), -3.
        actual = offset_resource_vectors(logs, logp, 0, 2)
        self.assertEqual(actual.shape, (16, 9, 2))
        for mask in range(16):
            shifted = np.array([logs[i, ((mask >> i) & 1):((mask >> i) & 1)+3]
                                for i in range(4)])
            expected = pair_resource_float(shifted, logp, 0, 2, 4, length=4)
            np.testing.assert_allclose(actual[mask], expected, rtol=0, atol=1e-13)
        # Exact noncommuting word products must lie below the relaxed DP.
        exact_tables = [[tuple(F(str(value)) for value in matrix.flat) for matrix in table]
                        for table in tables]
        for mask, sparse in product(range(16), product(range(3), repeat=4)):
            value = (F(1), F(0), F(0), F(1))
            for i, occupancy in enumerate(sparse):
                value = multiply(value, exact_tables[i][occupancy+((mask >> i) & 1)])
            vector = (value[0]+value[1]/8, value[2]+value[3]/8)
            for state in range(2):
                self.assertGreaterEqual(actual[mask, sum(sparse), state]+1e-12, log2(vector[state]))

    def test_cancellation_mixture_matches_explicit_mask_dp(self):
        logs = np.log2(np.array([[[[1+i+m, .5], [.25, 2+i+m]] for m in range(4)]
                                 for i in range(4)]))
        # Full eight-region comparison checks every actual resource target.
        logs = np.tile(logs, (2, 1, 1, 1))
        values = mixed_mask_vectors(logs, -2., groups=2, mu=2)
        for index, (single, double) in enumerate(categories(2)):
            for mask in (0, 1, 37, 128, 255):
                shifted = np.array([logs[i, ((mask >> i) & 1):((mask >> i) & 1)+3]
                                    for i in range(8)])
                expected = np.full(2, -np.inf)
                for cancelled in range(double+1):
                    lo, hi, limit, target = category_resources(single, double, cancelled, 2)
                    envelope = pair_resource_float(shifted, -2., lo, hi, limit)
                    weight = log2(comb(double, cancelled))+(double-cancelled)*log2(254)-double*log2(255)
                    expected = np.logaddexp2(expected, envelope[target]+weight)
                np.testing.assert_allclose(values[index, mask], expected, rtol=0, atol=1e-13)

    def test_dense_weights_precede_potential_quotient(self):
        vector_values = np.array([[[1.+mask**2, .25+((mask+7) % 16)**2]
                                    for mask in range(16)]])
        logp, potential = -2., F(1, 4)
        empty, nonempty = fixed_mask_branch_caps(np.log2(vector_values), 7, logp, 4, 4)
        counts = comparison_counts(7, 4, 4)
        exact = [sum(F(counts[mask], 64)*F(str(vector_values[0, mask, state]))
                     for mask in range(1, 16)) for state in range(2)]
        self.assertAlmostEqual(nonempty[0], log2(max(exact[0], exact[1]/potential)), places=12)
        expected_empty = max(F(str(vector_values[0, 0, 0])), F(str(vector_values[0, 0, 1]))/potential)/64
        self.assertAlmostEqual(empty[0], log2(expected_empty), places=12)
        premature = sum(F(counts[mask], 64)*max(F(str(vector_values[0, mask, 0])),
                         F(str(vector_values[0, mask, 1]))/potential) for mask in range(1, 16))
        self.assertLess(max(exact[0], exact[1]/potential), premature)

    def test_bivariate_tail_equals_exact_positive_word_enumeration(self):
        empty, nonempty = [F(2, 3), F(1, 5), F(3, 7)], [F(1, 11), F(2, 5)]
        exact = {}
        choices = [(flag, degree, value) for flag, poly in enumerate((empty, nonempty))
                   for degree, value in enumerate(poly)]
        for word in product(choices, repeat=3):
            x, y, mass = 0, 0, F(1)
            for flag, degree, value in word:
                x += degree
                y += flag
                mass *= value
            exact[x, y] = exact.get((x, y), 0)+mass
        for x_floor, y_floor in product(range(5), range(4)):
            value, by_y = bivariate_log_tail(np.log2(list(map(float, empty))),
                np.log2(list(map(float, nonempty))), 3, x_floor, y_floor)
            expected = sum(mass for (x, y), mass in exact.items() if x >= x_floor and y >= y_floor)
            if expected:
                self.assertAlmostEqual(value, log2(expected), places=12)
            else:
                self.assertEqual(value, -np.inf)
            for row in by_y:
                expected = sum(mass for (x, y), mass in exact.items()
                               if x >= x_floor and y == row['nonempty_dense_columns'])
                if expected:
                    self.assertAlmostEqual(row['log2_mass'], log2(expected), places=12)
                else:
                    self.assertEqual(row['log2_mass'], -np.inf)

    def test_cached_empty_scale_and_two_floor_normalization(self):
        q, a, groups, columns = 4, 3, 1, 3
        zero = column_log_coefficients(np.zeros(len(categories(groups))), groups, q)
        empty_probability = F(1, q**a)
        one = zero+log2(1-empty_probability)
        empty = zero+log2(empty_probability)
        powers = log_powers(zero, columns)
        value, _ = bivariate_log_tail(empty, one, columns, 4, 2,
                                      empty_powers=powers, empty_shift=log2(empty_probability))
        direct, _ = bivariate_log_tail(empty, one, columns, 4, 2)
        self.assertAlmostEqual(value, direct, places=12)
        support = sum(comb(6, degree)*(q-1)**degree for degree in range(4, 7))
        dense = sum(F(comb(3, degree))*(1-empty_probability)**degree*empty_probability**(3-degree)
                    for degree in range(2, 4))
        self.assertAlmostEqual(value, log2(support*dense), places=12)
        # Sparse shortened-support normalizer is separate from dense beta.
        normalized = F(support, (q-1)**(2*groups*(3-2)))*dense
        self.assertAlmostEqual(value-2*log2(q-1), log2(normalized), places=12)

    def test_scope_and_guard(self):
        regions = np.zeros((8, 18, 2, 2))
        with self.assertRaises(ValueError):
            mixed_mask_vectors(regions, -2, groups=16, mu=14)
        with self.assertRaises(ValueError):
            offset_resource_vectors(regions, -2, 0, 17)
        with self.assertRaises(ValueError):
            bivariate_log_tail([0], [0], columns=2, y_floor=3)


if __name__ == '__main__':
    unittest.main()
