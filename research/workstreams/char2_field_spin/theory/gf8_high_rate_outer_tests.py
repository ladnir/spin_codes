"""Regressions for outer-length exploration; no timing assertions."""
from fractions import Fraction as F
from itertools import product
import math
import unittest

from gf8_high_rate_outer_screen import (RETURN_MINIMA, bulk_source_log,
    expected_window_dimension, middle_bound, profile, sparse_bound, window_screen, window_witness)
from gf8_packet_width_screen import bound as old_middle
from gf8_packet_sparse_screen import bound as old_sparse
from gf8_rank_history_screen import source_row_log


class HighRateOuterTests(unittest.TestCase):
    def test_middle_length128_regression(self):
        point = (-.2, -.7, math.log(12), -.1)
        g = profile(128, 83)
        weight = round(.1*g['n'])
        actual = middle_bound(point, weight, 128, 83)
        expected = old_middle(point, weight, dimension=83, minima=RETURN_MINIMA)
        self.assertEqual(actual, expected)

    def test_sparse_length128_regression(self):
        point = (math.log(640), math.log(19))
        self.assertEqual(sparse_bound(4, point, 128, 83),
                         old_sparse(4, point, dimension=83, minima=RETURN_MINIMA))

    def test_requested_rate_uses_padding(self):
        g = profile(256, 192)
        self.assertEqual(g['rows'], 5504)
        self.assertEqual(g['nominal_rate'], '3/4')
        self.assertEqual(g['requested_rate'], '32/43')
        self.assertEqual(g['padded_k'], 1056768)
        self.assertEqual(g['n'], 1409024)
        self.assertFalse(g['whole_code'])

    def test_bulk_omits_at_most_one_bit_per_source_row(self):
        for length, dimension in ((128, 96), (192, 144), (256, 192)):
            for alpha, tau, tilt in product((.03, .1, .14), (2., 16., 255.), (.1, 1., 10.)):
                logtau, logtilt = math.log(tau), math.log(tilt)
                exact = source_row_log(alpha, logtau, logtilt, length, dimension, False)
                bulk = bulk_source_log(alpha, logtau, logtilt, length, dimension)
                self.assertGreaterEqual(exact+1e-12, bulk)
                self.assertLessEqual(exact-bulk, math.log(2)+1e-12)

    def test_window_expectation_by_support_enumeration(self):
        b, dimension, rows = 5, 3, 16
        occupancy = F(3, 8)
        exact = sum(max(0, sum(bits)-(b-dimension))*occupancy**sum(bits)
                    *(1-occupancy)**(b-sum(bits)) for bits in product((0, 1), repeat=b))*rows
        self.assertEqual(expected_window_dimension(b, dimension, rows, occupancy), exact)
        for length, dimension in ((128, 96), (256, 192)):
            result = window_screen(length, dimension)
            self.assertEqual(result['window_width'] % 8, 0)
            self.assertLessEqual(result['total_window_size'], result['threshold'])
            self.assertEqual(result['actual_obstruction_threshold'], 4*length+1)

    def test_invalid_geometries(self):
        for length, dimension in ((256, 256), (257, 192), (128, 0)):
            with self.assertRaises(ValueError):
                profile(length, dimension)
        with self.assertRaises(ValueError):
            expected_window_dimension(128, 96, 10, F(11, 10))

    def test_short_outer_realized_obstruction(self):
        result = window_witness(32, 24, seed=1)
        self.assertEqual(result['realized_shortened_dimension'], 176)
        self.assertEqual(result['kernel_dimension_lower'], 48)
        self.assertTrue(result['actual_route_obstruction'])
        self.assertLessEqual(result['total_window_size'], result['threshold'])


if __name__ == '__main__':
    unittest.main()
