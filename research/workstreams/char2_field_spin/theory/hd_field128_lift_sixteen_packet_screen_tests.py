"""Tiny local/nested coefficient tests; no distance gate or benchmark."""

from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np

import hd_field128_lift_sixteen_packet_screen as sixteen
from hd_field128_lift_four_packet_screen_tests import multiply


eight, four, old = sixteen.eight, sixteen.four, sixteen.old


class SixteenPacketTests(unittest.TestCase):
    def test_small_raw_scope_has53_records_and_retains_tiny_field128_terms(self):
        raw = sixteen.small_octets(2.**-200, 1 << 128, 64, old.Budget(10))
        self.assertEqual(len(raw[0]), 53)
        self.assertTrue(np.all(raw[0] == 0))
        self.assertLessEqual(int(np.max(raw[1])), 2)
        self.assertLess(np.min(raw[3][np.isfinite(raw[3])]), -2048)
        composed = four.four_from_pairs(raw, .01, 4, old.Budget(10))
        self.assertTrue(np.all(np.isfinite(composed[:, 0, 1])))

    def test_raw_sixteen_product_matches_exact_ordered_toy(self):
        matrices = ((F(1), F(0), F(0), F(1, 2)), (F(1, 4), F(1, 3), F(1, 5), F(2, 3)))
        exact_records = []
        for chosen in product(range(2), repeat=8):
            if sum(chosen) > 2:
                continue
            matrix = (F(1), F(0), F(0), F(1))
            for index in chosen:
                matrix = multiply(matrix, matrices[index])
            exact_records.append((sum(chosen), 3**sum(chosen), matrix))
        records = (np.zeros(len(exact_records), dtype=int), np.array([row[0] for row in exact_records]),
                   np.log2([row[1] for row in exact_records]),
                   np.array([[old.log_fraction(x) for x in row[2]] for row in exact_records]).reshape((-1, 2, 2)))
        alpha, maximum = .01, 8
        actual = four.four_from_pairs(records, alpha, maximum, old.Budget(10))
        expected = old._empty(maximum)
        for f, count, a in exact_records:
            for g, other, b in exact_records:
                values = np.array([old.NEG if not x else log2(count*other)+alpha*old.log_fraction(x)
                                   for x in multiply(a, b)])
                expected[:, 0, f+g] = np.logaddexp2(expected[:, 0, f+g], values)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=5e-13)

    def test_nested_partition_retains_complete_eight_and_dominates_its_square(self):
        for alpha, small in ((.01, 0), (.01, 2), (.01, 3), (1., 2)):
            actual, baseline = sixteen.sixteen_operator(.5, alpha, 4, 256, old.Budget(10), small_full=small)
            expected_eight, _ = eight.eight_operator(.5, alpha, 4, 256, old.Budget(10))
            np.testing.assert_allclose(baseline, expected_eight, rtol=0, atol=3e-13)
            square = old.polynomial_product(expected_eight, expected_eight, 4, old.Budget(10))
            self.assertTrue(np.all(actual <= square+1e-12))
            if alpha == 1 or small == 0:
                np.testing.assert_allclose(actual, square, rtol=0, atol=3e-13)
            else:
                self.assertTrue(np.any(actual < square-1e-5))

    def test_small_region_improves_complete_eight_region(self):
        for alpha in (.01, 1.):
            actual = sixteen.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            expected = eight.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            for index in actual:
                a, b = np.array(actual[index]), np.array(expected[index])
                self.assertTrue(np.all(a <= b+1e-10))
                if alpha == 1:
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)

    def test_guards_and_source_closure(self):
        for small in (-1, 4):
            with self.assertRaises(ValueError):
                sixteen.regional_screen(.5, .1, 2, small_full=small)
        with self.assertRaises(ValueError):
            sixteen.regional_screen(.5, .1, 65)
        with self.assertRaises(TimeoutError):
            sixteen.regional_screen(.5, .1, 2, max_seconds=1e-12)
        pins = sixteen.source_pins()
        for name in ('hd_field128_lift_sixteen_packet_screen.py', 'hd_field128_lift_eight_packet_screen.py',
                     'hd_field128_lift_four_packet_screen.py', 'hd_field128_lift_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
