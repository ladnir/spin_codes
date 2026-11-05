from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_shifted_good_pair_screen import good_step_coefficients, matrix_product
from hd_shifted_unoriented_screen import (
    bound, global_box_bound, low_degree_frequency, positive_pair_matrix,
    unsigned_coefficients, unsigned_histograms,
)


class UnorientedScreenTests(unittest.TestCase):
    def test_unsigned_coefficients_dominate_each_oriented_average(self):
        from hd_shifted_general_screen import general_coefficients
        original = general_coefficients(.5, .1, (0., 0., 0.))
        unsigned = unsigned_coefficients(.5, .1)
        for singleton in range(33):
            for first in range(singleton+1):
                self.assertTrue(np.all(unsigned[:, singleton, :]*comb(singleton, first)
                                       >= original[:, first, singleton-first, :]*(1-2e-15)))
        good = good_step_coefficients(.5, .1).reshape((2, 33, 4))
        np.testing.assert_allclose(unsigned[:, 0, :].T, good[0], rtol=2e-12, atol=1e-280)
        np.testing.assert_allclose(unsigned[:, 1, :].T, good[1], rtol=2e-12, atol=1e-280)

    def test_constant_kernel_cancels_orientation_multiplicity(self):
        coefficients = np.zeros((4, 33, 33, 33))
        for first in range(33):
            for second in range(33-first):
                for full in range(33-first-second):
                    singleton = first+second
                    count = comb(32, full)*comb(32-full, singleton)*comb(singleton, first)
                    coefficients[:, first, second, full] = count
        unsigned_coefficients.cache_clear()
        with patch('hd_shifted_unoriented_screen.general_coefficients', return_value=coefficients):
            unsigned = unsigned_coefficients(.4, .2)
        unsigned_coefficients.cache_clear()
        for singleton in range(33):
            for full in range(33-singleton):
                expected = comb(32, full)*comb(32-full, singleton)
                np.testing.assert_allclose(unsigned[:, singleton, full], expected, rtol=2e-15)

    def test_unsigned_histograms_ignore_even_odd_split(self):
        from hd_shifted_general_count import matching_histograms
        for even in range(450, 543):
            oriented = matching_histograms(4096, even, 992-even, 450, 420)
            self.assertEqual(unsigned_histograms(4096, 992, 450, 420),
                             tuple((row[0], row[1]+row[2], row[3]) for row in oriented))

    def test_full_degree_64_is_folded(self):
        coefficients = np.zeros((3, 65, 2, 2))
        coefficients[2, 64] = ((1, 2), (3, 4))
        result = low_degree_frequency(coefficients, (-1., 0.), 64, (2, 64))
        roots = np.exp(-2j*np.pi*np.arange(64)/64)
        expected = np.array((1, 2, 3, 4))[:, None]*(roots[None, :]/2)**2
        np.testing.assert_allclose(result[:, :, 0], expected, atol=1e-14)
        np.testing.assert_allclose(result[:, :, 37], expected, atol=1e-14)

    def test_global_corner_bound_dominates_small_integer_domain(self):
        ranges = ((1023, 1025), (127, 129), (127, 129))
        odds = ((-1.4, -3.5), (-1.4, -3.5))
        envelope = global_box_bound(ranges, .5, .14, odds)
        for point in product(*(range(lo, hi+1) for lo, hi in ranges)):
            actual = bound(point, .5, .14, odds, method='global', matching_count=True)
            self.assertLessEqual(actual['log2_upper_screen'], envelope['class_log2_upper_screen']+1e-7)
        self.assertEqual(envelope['class_count'], 27)


if __name__ == '__main__':
    unittest.main()
