"""Numerical bookkeeping checks without regenerating the large state census."""
from fractions import Fraction as Q
import unittest

import tune


class BookkeepingTests(unittest.TestCase):
    def test_exact_map_and_saved_baseline(self):
        record, margins = tune.retained()
        self.assertEqual(record['map_record']['map_sha256'], tune.MAP_SHA256)
        self.assertAlmostEqual(margins[89], 44.841583, places=6)
        self.assertAlmostEqual(tune.union_margin(margins.values()), 42.090840060352, places=8)

    def test_log_union_retains_tiny_terms_without_underflow(self):
        self.assertAlmostEqual(tune.union_margin([60, 60, 2000]), 59)
        self.assertAlmostEqual(tune.union_margin([2000, 2000]), 1999)

    def test_dyadic_log_does_not_convert_to_float_probability(self):
        self.assertEqual(tune.endpoint_margin([1, -2000]), 2000)
        self.assertEqual(tune.endpoint_margin([1 << 255, -2255]), 2000)

    def test_search_budget_and_exact_rationals(self):
        self.assertEqual(tune.checked_tilts(['.085']), (Q(17, 200),))
        for values in ([], ['.085', '17/200'], ['0'], [str(i) for i in range(1, 22)]):
            with self.assertRaises(ValueError):
                tune.checked_tilts(values)


if __name__ == '__main__':
    unittest.main()
