"""Tiny exact and coefficient tests; no construction-size distance gate."""

from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

import numpy as np

import hd_field128_lift_four_packet_screen as four
import hd_field128_lift_screen as old


def multiply(a, b):
    return tuple(sum(a[2*r+k]*b[2*k+c] for k in range(2)) for r in range(2) for c in range(2))


def explicit(rows, length, alpha, maximum, predicate):
    result = old._empty(maximum)
    identity = tuple(F(int(i in (0, 3))) for i in range(4))
    for chosen in product(rows, repeat=length):
        if not predicate(chosen):
            continue
        s, f = sum(row[0] for row in chosen), sum(row[1] for row in chosen)
        if s+2*f > maximum:
            continue
        count, matrix = 1, identity
        for _, _, n, raw in chosen:
            count *= n
            matrix = multiply(matrix, raw)
        values = np.array([old.NEG if not x else log2(count)+alpha*old.log_fraction(x) for x in matrix])
        result[:, s, f] = np.logaddexp2(result[:, s, f], values)
    return result


class FourPacketTests(unittest.TestCase):
    def toy(self):
        rows = ((0, 0, 1, (F(1), F(0), F(0), F(1, 2))),
                (1, 0, 3, (F(1, 4), F(1, 3), F(1, 5), F(2, 3))),
                (0, 1, 2, (F(1, 7), F(3, 5), F(2, 9), F(3, 8))))
        s, f, counts = (np.array([row[i] for row in rows], dtype=int) for i in range(3))
        logs = np.array([[old.log_fraction(x) for x in row[3]] for row in rows]).reshape((-1, 2, 2))
        return rows, s, f, counts, logs

    def test_positive_partition_matches_explicit_exact_ordered_products(self):
        rows, s, f, counts, logs = self.toy()
        for alpha in (.25, 1.):
            selected, rest, pairs = four.pair_partition(s, f, counts, logs, alpha, 6, old.Budget(10), 0)
            raw_four = four.four_from_pairs(pairs, alpha, 6, old.Budget(10))
            expected_vv = explicit(rows, 2, alpha, 6, lambda chosen: all(row[1] == 0 for row in chosen))
            expected_rest = explicit(rows, 2, alpha, 6, lambda chosen: any(row[1] != 0 for row in chosen))
            expected_four = explicit(rows, 4, alpha, 6, lambda chosen: all(row[1] == 0 for row in chosen))
            for actual, expected in ((selected, expected_vv), (rest, expected_rest), (raw_four, expected_four)):
                np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-13)

    def test_four_operator_is_positive_and_dominated_by_old_square(self):
        rows, s, f, counts, logs = self.toy()
        for alpha in (.25, 1.):
            selected, rest, pairs = four.pair_partition(s, f, counts, logs, alpha, 6, old.Budget(10), 0)
            c4 = four.four_from_pairs(pairs, alpha, 6, old.Budget(10))
            refined, old_two = four.four_operator(selected, rest, c4, 6, old.Budget(10))
            old_four = old.polynomial_product(old_two, old_two, 6, old.Budget(10))
            expected = c4.copy()
            for left, right in ((selected, rest), (rest, selected), (rest, rest)):
                expected = np.logaddexp2(expected, old.polynomial_product(left, right, 6, old.Budget(10)))
            np.testing.assert_allclose(refined, expected, rtol=0, atol=2e-13)
            self.assertTrue(np.all(refined <= old_four+1e-12))
            if alpha == 1:
                np.testing.assert_allclose(refined, old_four, rtol=0, atol=2e-13)
            else:
                self.assertTrue(np.any(refined < old_four-1e-5))
        self.assertNotEqual(multiply(rows[1][3], rows[2][3]), multiply(rows[2][3], rows[1][3]))

    def test_literal_two_packet_partition_matches_frozen_coefficients(self):
        maximum, z, alpha, q = 3, .5, .2, 256
        selected, rest, c4 = four.positive_parts(z, alpha, maximum, q, old.Budget(10), 2)
        actual = np.logaddexp2(selected, rest)
        single = old.single_coefficients(z, alpha, maximum, q, old.Budget(10))
        good, dirty = single.copy(), single.copy()
        good[:, 2:, :], dirty[:, :2, :] = old.NEG, old.NEG
        expected = old.good_pair_coefficients(z, alpha, maximum, q, old.Budget(10))
        for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
            expected = np.logaddexp2(expected, old.polynomial_product(left, right, maximum, old.Budget(10)))
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-13)
        refined, _ = four.four_operator(selected, rest, c4, maximum, old.Budget(10))
        expected_square = old.polynomial_product(expected, expected, maximum, old.Budget(10))
        self.assertTrue(np.all(refined <= expected_square+1e-12))

    def test_small_regional_prefix_improves_or_equals_frozen_two_packet(self):
        for alpha in (.5, 1.):
            actual = four.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            expected = old.regional_screen(.5, alpha, 3, q=256, max_seconds=10)
            for index in actual:
                new, previous = np.array(actual[index]), np.array(expected[index])
                self.assertTrue(np.all(new <= previous+1e-10))
                if alpha == 1:
                    np.testing.assert_allclose(new, previous, rtol=0, atol=1e-10)

    def test_geometry_scope_and_very_small_field128_raw_terms(self):
        s, f, counts, indices = old.good_geometry()
        keep = f <= 2
        self.assertEqual(int(keep.sum()), 10)
        # This is only ten local matrices and their tiny degree2 coefficients.
        logs = old.local_logs(2.**-200, 1 << 128, old.OCC[indices[keep]])
        selected, rest, pairs = four.pair_partition(s[keep], f[keep], counts[keep], logs, .01, 2, old.Budget(10), 2)
        c4 = four.four_from_pairs(pairs, .01, 2, old.Budget(10))
        self.assertTrue(np.all(np.isfinite(c4[:, 1, 0])))
        self.assertTrue(np.all(np.isneginf(rest)))
        self.assertLess(np.min(pairs[3][np.isfinite(pairs[3])]), -2048)

    def test_scope_guards_and_source_pins(self):
        for cutoff in (-1, 3):
            with self.assertRaises(ValueError):
                four.regional_screen(.5, .1, 2, max_v_full=cutoff)
        with self.assertRaises(ValueError):
            four.regional_screen(.5, .1, 65)
        with self.assertRaises(TimeoutError):
            four.regional_screen(.5, .1, 2, max_seconds=1e-12)
        pins = four.source_pins()
        self.assertIn('hd_field128_lift_four_packet_screen.py', pins)
        self.assertIn('hd_field128_lift_screen.py', pins)
        self.assertIn('hd_common_refresh_principal.py', pins)


if __name__ == '__main__':
    unittest.main()
