"""Independent checks for nonzero amplitudes and regional scalar potentials."""
from fractions import Fraction
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_grouped_subfield_screen import density_bits
from hd_randomized_outer_state_profile import multiply, power
from hd_xor_column_screen import one_group_region_words
from hd_xor_multigroup import (exact_group_moment, homogeneous_bound,
    log_polynomial_product, nonzero_coefficients, nonzero_kernel_probability,
    nonzero_matrices, potential_envelope, profile_arrays, region_operators,
    directed_one_row_receipt, replay_one_row_128, region_operators_arb_float)


class XorMultigroupTests(unittest.TestCase):
    def test_small_field_kernel_counts(self):
        # Independent exhaustive GF4 Vandermonde syndromes, polynomial x^2+x+1.
        def mul(a, b):
            out = 0
            while b:
                if b & 1:
                    out ^= a
                a <<= 1
                if a & 4:
                    a ^= 7
                b >>= 1
            return out
        for rank in (1, 2):
            for occupancy in range(5):
                quiet = 0
                for amplitudes in product(range(1, 4), repeat=occupancy):
                    syndrome = [0]*rank
                    for node, value in enumerate(amplitudes):
                        factor = 1
                        for check in range(rank):
                            syndrome[check] ^= mul(value, factor)
                            factor = mul(factor, node)
                    quiet += not any(syndrome)
                self.assertEqual(nonzero_kernel_probability(occupancy, rank, 4),
                                 Fraction(quiet, 3**occupancy))

    def test_feedback_floor_and_first_cancellation(self):
        self.assertEqual(nonzero_kernel_probability(0), 1)
        for occupancy in range(1, 5):
            self.assertEqual(nonzero_kernel_probability(occupancy), 0)
        self.assertEqual(nonzero_kernel_probability(5), Fraction(1, 255**4))
        for occupancy in range(5, 17):
            probability = nonzero_kernel_probability(occupancy)
            self.assertGreater(probability, 0)
            self.assertLessEqual(probability, Fraction(256, 255)**occupancy/Fraction(256)**4)

    def test_first_row_mass_and_multiblock_cancellation(self):
        z = .9375
        occupancies, _ = profile_arrays()
        degrees, _, matrices = nonzero_matrices(z)
        b = 2**density_bits(256, 16)
        np.testing.assert_allclose(matrices[:, 0, 0]+matrices[:, 0, 1]/b,
                                   z**degrees, rtol=2e-14)
        lookup = {tuple(occ): index for index, occ in enumerate(occupancies)}
        for profile in ((0, 0, 0, 1), (0, 0, 0, 4), (1, 1, 1, 1), (0, 0, 1, 5)):
            self.assertEqual(matrices[lookup[profile], 0, 0], 0.)
        expected = z**10*float(nonzero_kernel_probability(5))**2
        self.assertAlmostEqual(matrices[lookup[(0, 0, 5, 5)], 0, 0]/expected, 1.)

    def test_profile_coefficients_count_all_subsets(self):
        degrees, counts, matrices = nonzero_matrices(1.)
        b = 2**density_bits(256, 16)
        coefficients = np.exp2(nonzero_coefficients(1.))
        for degree in range(65):
            self.assertEqual(sum(int(c) for c in counts[degrees == degree]), comb(64, degree))
            self.assertAlmostEqual((coefficients[degree, 0, 0]+coefficients[degree, 0, 1]/b)
                                   /comb(64, degree), 1., places=12)

    def test_log_matrix_polynomial_noncommuting(self):
        left = np.array([[[1, 2], [0, 3]], [[2, 0], [1, 1]]], float)
        right = np.array([[[0, 1], [2, 1]], [[3, 1], [1, 0]], [[1, 0], [0, 2]]], float)
        expected = np.zeros((4, 2, 2))
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                expected[i+j] += a @ b
        with np.errstate(divide='ignore'):
            actual = np.exp2(log_polynomial_product(np.log2(left), np.log2(right), 3))
        np.testing.assert_allclose(actual, expected, rtol=1e-14)

    def test_region_one_agrees_with_existing_forced_nonzero(self):
        z = .995
        regions = region_operators(z, 1)
        word = np.array(multiply(tuple(map(tuple, regions[1])),
                                 power(tuple(map(tuple, regions[0])), 7)))
        previous = np.array(one_group_region_words(z)[1])
        finite = np.isfinite(word) & np.isfinite(previous)
        self.assertTrue(np.array_equal(np.isfinite(word), np.isfinite(previous)))
        self.assertLessEqual(float(np.max(word[finite]-previous[finite])), 1e-10)
        np.testing.assert_allclose(word[0, 1], previous[0, 1], rtol=1e-13)

    def test_arb_proposal_engine_matches_log_convolution(self):
        actual = region_operators_arb_float(.983, 32)
        expected = region_operators(.983, 32)
        self.assertTrue(np.array_equal(np.isfinite(actual), np.isfinite(expected)))
        finite = np.isfinite(actual)
        np.testing.assert_allclose(actual[finite], expected[finite], rtol=2e-12, atol=1e-10)

    def test_potential_handles_correlated_chronological_paths(self):
        # Every sequence is allowed; no independence assumption enters.
        matrices = np.array([[[1, 0], [0, .2]], [[.1, .2], [.01, .3]],
                             [[.03, .06], [.004, .1]]])
        with np.errstate(divide='ignore'):
            mark, _ = potential_envelope(np.log2(matrices), -2.)
        endpoint_cost = 4.
        for occupancies in product(range(3), repeat=4):
            run = np.eye(2)
            for occupancy in occupancies:
                run = run @ matrices[occupancy]
            self.assertLessEqual(run[0].sum(), endpoint_cost*2**(mark*sum(occupancies))*(1+1e-13))

    def test_outer_moments_include_nonzero_rows(self):
        self.assertAlmostEqual(exact_group_moment(1, 0), log2(256**16-1), places=12)
        self.assertAlmostEqual(exact_group_moment(2, 0), 2*log2(256**16-1), places=12)

    def test_fixed_multigroup_witnesses(self):
        # Fixed witnesses, not optimizer regression: include group and row choices.
        for groups, z, potential, target in ((2, .99285, -62.72, -880),
                                             (4, .99214, -62.76, -1010),
                                             (16, .98883, -62.9, -1540),
                                             (128, .96777, -63.53, -3500)):
            value, _, _ = homogeneous_bound(groups, 1, z, potential)
            self.assertLess(value, target)

    def test_directed_receipt_rejects_false_potential(self):
        with self.assertRaisesRegex(ArithmeticError, 'potential check failed'):
            directed_one_row_receipt(2, 2, '19857/20000', '-251/4', '-100')

    def test_saved_directed_coverage(self):
        from flint import arb
        receipt = replay_one_row_128()
        self.assertEqual(receipt['active_group_interval'], [1, 128])
        self.assertEqual(receipt['interval_count'], 29)
        self.assertTrue(arb(receipt['failure_log2_upper']) < -97)
        self.assertFalse(receipt['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
