"""Small exhaustive checks and independent fixed-threshold route witnesses."""
from fractions import Fraction
from itertools import combinations, product
from math import comb, log
import unittest

import route_counts


class RouteCountsTests(unittest.TestCase):
    def test_small_subset_counts(self):
        for slots, windows in ((6, 2), (8, 4), (9, 3)):
            for q in range(slots+1):
                expected = [0]*(slots//windows+1)
                for selected in combinations(range(slots), q):
                    expected[len({x//windows for x in selected})] += 1
                self.assertEqual(route_counts.counts(q, slots, windows), tuple(expected))

    def test_small_tail_chernoff_and_union(self):
        slots, windows, regions = 6, 2, 2
        for q in range(slots+1):
            positions = list(combinations(range(slots), q))
            histogram = [0]*(regions*slots//windows+1)
            for selected in product(positions, repeat=regions):
                occupied = sum(len({x//windows for x in region}) for region in selected)
                histogram[occupied] += 1
            for h in range(len(histogram)+1):
                union = min(Fraction(1), Fraction(comb(slots, q)*sum(histogram[:h]), len(positions)**regions))
                for x in (Fraction(1, 5), Fraction(1, 2), Fraction(4, 5), Fraction(1)):
                    bound = route_counts.rational_bad_route_bound(
                        q, h, x, regions=regions, slots=slots, windows=windows)
                    self.assertGreaterEqual(bound, union)

    def test_float_matches_rational(self):
        for q, h, x in ((64, 1057, Fraction(1, 5)), (119, 1483, Fraction(3, 20)),
                        (128, 1531, Fraction(1, 7))):
            bound = route_counts.rational_bad_route_bound(q, h, x)
            exact_log = (log(bound.denominator)-log(bound.numerator))/log(2)
            self.assertAlmostEqual(route_counts.bad_route_margin(q, h, -log(float(x))), exact_log, places=8)

    def test_selected_exact_sixty_bit_thresholds(self):
        for q, h, x in ((32, 643, '3/16'), (64, 1057, '1/5'), (119, 1483, '3/20'),
                        (128, 1531, '1/7'), (192, 1765, '1/15'), (256, 1888, '1/46')):
            witness = route_counts.rational_witness(q, h, x)
            self.assertTrue(witness['exact_integer_check_passed'])
            self.assertGreaterEqual(witness['certified_margin_bits_floor'], 60)
            self.assertFalse(witness['whole_code_certificate'])

    def test_bad_inputs(self):
        for arguments in ((1, 5, 2), (7, 6, 2), (-1, 6, 2), (True, 6, 2)):
            with self.assertRaises(ValueError):
                route_counts.counts(*arguments)
        for x in (0, -1, '3/2', .5):
            with self.assertRaises(ValueError):
                route_counts.rational_bad_route_bound(3, 3, x, regions=2, slots=6, windows=2)


if __name__ == '__main__':
    unittest.main()
