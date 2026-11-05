"""Tiny nested-partition tests; no distance gate or benchmark."""

from fractions import Fraction as F
import unittest

import numpy as np

import hd_field128_lift_eight_packet_screen as eight
import hd_field128_lift_four_packet_screen as four
from hd_field128_lift_four_packet_screen_tests import explicit


old = four.old


class EightPacketTests(unittest.TestCase):
    def test_raw_quartets_and_eighth_power_match_exact_ordered_toy(self):
        rows = ((0, 0, 1, (F(1), F(0), F(0), F(1, 2))),
                (0, 1, 3, (F(1, 4), F(1, 3), F(1, 5), F(2, 3))))
        s, f, counts = (np.array([row[i] for row in rows], dtype=int) for i in range(3))
        logs = np.array([[old.log_fraction(x) for x in row[3]] for row in rows]).reshape((-1, 2, 2))
        raw = (s, f, np.log2(counts), logs)
        raw = eight.combine_raw(eight.combine_raw(raw, 8, old.Budget(10)), 8, old.Budget(10))
        self.assertEqual(len(raw[0]), 16)
        for alpha in (.25, 1.):
            actual = four.four_from_pairs(raw, alpha, 8, old.Budget(10))
            expected = explicit(rows, 8, alpha, 8, lambda _: True)
            np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-13)

    def test_clean_geometry_and_raw_field128_exponent_range(self):
        raw = eight.clean_quartets(2.**-200, 1 << 128, 16, old.Budget(10))
        self.assertEqual(len(raw[0]), 256)
        self.assertTrue(np.all(raw[0] == 0))
        self.assertLess(np.min(raw[3][np.isfinite(raw[3])]), -2048)
        cap = four.four_from_pairs(raw, .01, 2, old.Budget(10))
        self.assertTrue(np.all(np.isfinite(cap[:, 0, 1])))
        self.assertTrue(np.all(np.isneginf(cap[:, 1:, :])))

    def test_nested_partition_retains_frozen_four_and_improves_its_square(self):
        maximum = 4
        for alpha in (.01, 1.):
            actual, baseline = eight.eight_operator(.5, alpha, maximum, 256, old.Budget(10))
            selected, rest, c4 = four.positive_parts(.5, alpha, maximum, 256, old.Budget(10))
            expected_four, _ = four.four_operator(selected, rest, c4, maximum, old.Budget(10))
            np.testing.assert_allclose(baseline, expected_four, rtol=0, atol=2e-13)
            expected_eight = old.polynomial_product(expected_four, expected_four, maximum, old.Budget(10))
            self.assertTrue(np.all(actual <= expected_eight+1e-12))
            if alpha == 1:
                np.testing.assert_allclose(actual, expected_eight, rtol=0, atol=2e-13)
            else:
                self.assertTrue(np.any(actual < expected_eight-1e-5))

    def test_small_region_is_dominated_by_complete_four_packet_region(self):
        for alpha in (.5, 1.):
            actual = eight.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            expected = four.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            for index in actual:
                a, b = np.array(actual[index]), np.array(expected[index])
                self.assertTrue(np.all(a <= b+1e-10))
                if alpha == 1:
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)

    def test_guards_and_complete_source_closure(self):
        with self.assertRaises(ValueError):
            eight.regional_screen(.5, .1, 65)
        with self.assertRaises(ValueError):
            eight.regional_screen(.5, .1, 2, max_v_full=3)
        with self.assertRaises(TimeoutError):
            eight.regional_screen(.5, .1, 2, max_seconds=1e-12)
        pins = eight.source_pins()
        for name in ('hd_field128_lift_eight_packet_screen.py', 'hd_field128_lift_four_packet_screen.py',
                     'hd_field128_lift_screen.py', 'hd_common_refresh_principal.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
