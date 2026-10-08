"""Independent source counts, exact geometry, and baseline agreement."""

from decimal import Decimal
from fractions import Fraction as F
import math
import unittest

import numpy as np

from gf16_sizes import class_upper as baseline_upper, screen_bound as baseline_screen
from outer_rate_screen import geometry, source_eta, exact_source, class_upper, screen_bound
from q1_batched_screen import fixed_support_multiplicity
from subfield_certificate import Interval


class OuterRateTests(unittest.TestCase):
    def test_natural_geometry_and_requested_subcode(self):
        expected = [(8, 1048576, 131072, 2097152, 0),
                    (7, 1048656, 149808, 2396928, 80),
                    (6, 1048608, 174768, 2796288, 32)]
        for dimension, k, rows, n, padding in expected:
            actual = geometry(dimension=dimension)
            self.assertEqual((actual['k'], actual['rows'], actual['n'], actual['padding']),
                             (k, rows, n, padding))
            self.assertEqual(rows % 16, 0)
            self.assertEqual(k, rows*dimension)
            self.assertEqual(actual['threshold'], (n+9)//10-1)
            self.assertEqual(k-padding, 1 << 20)
            self.assertLess(k-16*dimension, 1 << 20)
        for k, dimension in [(0, 7), (-1, 7), (1, 0), (1, 9)]:
            with self.assertRaises(ValueError):
                geometry(k, dimension)

    def test_mds_spectrum_total_and_minimum_support(self):
        for dimension in [6, 7, 8]:
            counts = [fixed_support_multiplicity(16, dimension, w, 16) for w in range(1, 17)]
            self.assertTrue(all(x == 0 for x in counts[:16-dimension]))
            self.assertEqual(counts[16-dimension], 15)
            self.assertEqual(sum(math.comb(16, w)*counts[w-1] for w in range(1, 17)),
                             16**dimension-1)

    def test_source_envelope_has_correct_codimension(self):
        for dimension in [6, 7, 8]:
            for rho in [F(4, 5), F(871, 1000), F(24, 25), F(2)]:
                eta = source_eta(dimension, rho)
                for weight in range(17-dimension, 17):
                    count = fixed_support_multiplicity(16, dimension, weight, 16)
                    self.assertLessEqual(count, eta*rho**weight*16**(weight-(16-dimension)))

    def test_directed_exact_source_against_rational_polynomial(self):
        # Use noncommuting matrices so that time-ordering remains tested.
        matrices = [np.array([[F(1, 7), F(2, 9)], [F(1, 11), F(2, 5)]], dtype=object),
                    np.array([[F(1, 17), F(1, 13)], [F(2, 23), F(1, 7)]], dtype=object)]
        coefficients = [np.array([F(1), F(0)], dtype=object)]
        for _ in range(16):
            following = [np.array([F(0), F(0)], dtype=object) for _ in range(len(coefficients)+1)]
            for weight, vector in enumerate(coefficients):
                following[weight] += vector @ matrices[0]
                following[weight+1] += vector @ matrices[1]
            coefficients = following
        upper = [[Interval.rational(x.numerator, x.denominator).hi for x in matrix.flat]
                 for matrix in matrices]
        for dimension in [6, 7, 8]:
            expected = sum(sum(coefficients[w])*fixed_support_multiplicity(16, dimension, w, 16)
                           for w in range(17-dimension, 17))
            actual = F(exact_source(upper, 1, dimension, directed=True))
            self.assertGreaterEqual(actual, expected)
            self.assertLess(actual-expected, expected*F(1, 10**70))

    def test_preserved_half_rate_bounds(self):
        for active in [1, 4, 8]:
            expected = baseline_screen(20, active, .999, state=4, source_rho=.871)
            actual = screen_bound(active, .999, dimension=8, state=4, rho=.871)
            self.assertAlmostEqual(actual, expected, places=6)
        expected = baseline_upper(20, 1, 999, 1000, state=4)
        actual = class_upper(1, F(999, 1000), dimension=8, state=4)
        self.assertEqual(Decimal(actual['failure_upper']), Decimal(expected['failure_upper']))
        self.assertFalse(actual['complete'])


if __name__ == '__main__':
    unittest.main()
