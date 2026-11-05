from fractions import Fraction as F
from itertools import product
from math import log2
import unittest

from hd_two_sum_fiber_activity import (compare, quiet_word_logs,
    row_weight_coefficients, weighted_row_enumerator)
from hd_two_sum_syndrome_rank_tests import mul


class ActivityTests(unittest.TestCase):
    def test_exact_weighted_rs_row(self):
        points = (0, 1, 2)
        words = [tuple(a ^ mul(b, x) for x in points)
                 for a, b in product(range(4), repeat=2) if a or b]
        counts = tuple(sum(sum(x != 0 for x in c) == w for c in words)
                       for w in range(4))
        self.assertEqual(row_weight_coefficients(4, 3, 2), counts)
        markers = (F(1, 2), F(2), F(3, 7))
        direct = F(0)
        for word in words:
            term = F(1)
            for value, marker in zip(word, markers):
                if value:
                    term *= marker
            direct += term
        self.assertEqual(weighted_row_enumerator(markers, q=4, dimension=2), direct)

    def test_exact_two_column_repetition_quiet_count(self):
        # Three GF4 repetition rows: each of two active columns has3
        # nonzero solutions to its two moment checks; denominator3^3.
        for marker in (F(1), F(1, 2), F(2)):
            row = quiet_word_logs(3, log2(marker), q=4, length=2, dimension=1)
            self.assertAlmostEqual(row['log2_upper'], log2(F(1, 3)*marker**6), places=12)
            self.assertTrue(row['all_rows_seen'])

    def test_positive_fractional_count_toy(self):
        # Full row activity forces all3 rows active in both columns.
        row = quiet_word_logs(3, -2., .5, q=4, length=2, dimension=1)
        expected = (2*.5-3*.5)*log2(3)+6*(-2.)
        self.assertAlmostEqual(row['log2_upper'], expected, places=12)

    def test_alpha_one_ordering_control(self):
        row = compare(8, -.2, 1.)
        self.assertAlmostEqual(row['old_mean_then_alpha_log2'],
                               row['old_column_alpha_log2'], places=10)
        self.assertFalse(row['regional_envelope_recovered'])
        self.assertFalse(row['distance_certificate'])


if __name__ == '__main__':
    unittest.main()
