from math import comb
import unittest

import numpy as np

from hd_shifted_general_screen import general_coefficients
from hd_shifted_good_pair_screen import good_pair_coefficients, good_step_coefficients, paired_frequency
from hd_shifted_swap_screen import (CONSTRUCTION, SwapPositiveWitness, bound,
    positive_pair_matrix, regional, swap_coefficients)
from hd_shifted_unoriented_screen import unsigned_coefficients, unsigned_histograms


class SwapScreenTests(unittest.TestCase):
    def test_coefficients_are_fair_orientation_average_and_below_old_envelope(self):
        original = general_coefficients(.5, .1, (0., 0., 0.))
        actual = swap_coefficients(.5, .1)
        expected = np.zeros_like(actual)
        for singleton in range(33):
            for first in range(singleton+1):
                expected[:, singleton, :] += original[:, first, singleton-first, :]/2**singleton
        np.testing.assert_allclose(actual, expected, rtol=3e-14, atol=1e-280)
        self.assertTrue(np.all(actual <= unsigned_coefficients(.5, .1)*(1+3e-14)+1e-280))

    def test_good_coefficients_unchanged(self):
        actual = swap_coefficients(.5, .1)
        good = good_step_coefficients(.5, .1).reshape((2, 33, 4))
        np.testing.assert_allclose(actual[:, :2, :].transpose(1, 2, 0), good,
                                   rtol=3e-14, atol=1e-280)

    def test_fast_positive_matches_full_polynomial(self):
        for odds in ((-4., -2.), (-12., -4.), (1., 2.)):
            dirty = swap_coefficients(.5, .1).copy()
            dirty[:, :2, :] = 0
            marker = np.exp2(np.arange(33)[:, None]*odds[0]+np.arange(33)[None, :]*odds[1])
            dirty = np.sum(dirty*marker, axis=(1, 2))
            def evaluate(coefficients):
                marker = np.exp2(np.arange(len(coefficients))[:, None]*odds[0]
                                 +np.arange(coefficients.shape[1])[None, :]*odds[1])
                return np.sum(coefficients*marker[:, :, None, None], axis=(0, 1)).ravel()
            expected = paired_frequency(dirty, evaluate(good_step_coefficients(.5, .1)),
                                         evaluate(good_pair_coefficients(.5, .1)))
            np.testing.assert_allclose(positive_pair_matrix(.5, .1, odds).ravel(), expected, rtol=5e-14)

    def test_alpha_one_reduces_to_squared_single_step(self):
        odds = (-3., -2.)
        coeff = swap_coefficients(.5, 1.)
        marker = np.exp2(np.arange(33)[:, None]*odds[0]+np.arange(33)[None, :]*odds[1])
        step = np.sum(coeff*marker, axis=(1, 2)).reshape((2, 2))
        np.testing.assert_allclose(positive_pair_matrix(.5, 1., odds), step@step, rtol=5e-13)

    def test_witness_preserves_single_matching_envelope_and_labels(self):
        profile, odds = (1024, 400, 400), ((-3.5, -2.3),)*2
        witness = SwapPositiveWitness.create(.5, .1, odds, profile)
        envelope, _ = witness.evaluate(tuple((x, x) for x in profile))
        actual = bound(profile, .5, .1, odds)
        self.assertLessEqual(actual['log2_upper_screen'], envelope+1e-7)
        self.assertEqual(witness.serializable()['construction'], CONSTRUCTION)
        self.assertTrue(witness.serializable()['pair_swaps'])

    def test_alias_below_positive_and_validates_inputs(self):
        histogram = unsigned_histograms(4096, 1024, 400, 400)[0]
        positive, _ = regional(.5, .1, histogram, (-3.5, -2.3), method='global')
        alias, _ = regional(.5, .1, histogram, (-3.5, -2.3), grid=64)
        self.assertTrue(np.all(alias <= positive+1e-7))
        with self.assertRaises(ValueError):
            regional(.5, .1, histogram, (-3.5, -2.3), grid=63)
        with self.assertRaises(ValueError):
            swap_coefficients(1., .1)

    def test_large_finite_boundary_markers_do_not_overflow_perron_scaling(self):
        value, _ = regional(1/32, .25, (0, 4096, 0), (13., -40.), grid=64)
        self.assertTrue(np.all(np.isfinite(value)))


if __name__ == '__main__':
    unittest.main()
