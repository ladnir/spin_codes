"""Exhaustive small capped-occupancy tests and selected exact witnesses."""
from fractions import Fraction
from itertools import combinations, product
from math import comb, log
import unittest

import cap_counts
import route_counts


def statistic(selected, slots, windows, cap):
    counts = [0]*(slots//windows)
    for slot in selected:
        counts[slot//windows] += 1
    return sum(min(x, cap) for x in counts)


class CapCountsTests(unittest.TestCase):
    def test_region_moments(self):
        slots, windows = 9, 3
        for q in range(slots+1):
            selected = list(combinations(range(slots), q))
            for cap in range(1, windows+1):
                for x in (Fraction(1, 3), Fraction(2, 5), Fraction(1)):
                    expected = sum(x**statistic(s, slots, windows, cap) for s in selected)/len(selected)
                    self.assertEqual(cap_counts.one_region_moment(q, x, cap=cap,
                        slots=slots, windows=windows), expected)

    def test_small_tail_union(self):
        slots, windows, regions = 6, 3, 2
        for q in range(slots+1):
            selected = list(combinations(range(slots), q))
            for cap in range(1, windows+1):
                histogram = [0]*(regions*(slots//windows)*cap+1)
                for regions_selected in product(selected, repeat=regions):
                    h = sum(statistic(s, slots, windows, cap) for s in regions_selected)
                    histogram[h] += 1
                for h in range(len(histogram)+1):
                    tail = min(Fraction(1), Fraction(comb(slots, q)*sum(histogram[:h]), len(selected)**regions))
                    bound = cap_counts.rational_bad_route_bound(q, h, '2/5', cap=cap,
                        regions=regions, slots=slots, windows=windows)
                    self.assertGreaterEqual(bound, tail)

    def test_cap_one_agrees_with_occupied_steps(self):
        for q, h, x in ((64, 1057, '1/5'), (119, 1483, '3/20'), (128, 1531, '1/7')):
            self.assertEqual(cap_counts.rational_bad_route_bound(q, h, x, cap=1),
                route_counts.rational_bad_route_bound(q, h, x))

    def test_exact_sixty_bit_witnesses_and_floating_check(self):
        for q, h, x in ((64, 1615, '3/14'), (119, 2559, '9/35'), (128, 2680, '9/35')):
            witness = cap_counts.rational_witness(q, h, x)
            self.assertTrue(witness['exact_integer_check_passed'])
            self.assertGreaterEqual(witness['certified_margin_bits_floor'], 60)
            bound = cap_counts.rational_bad_route_bound(q, h, x)
            margin = (log(bound.denominator)-log(bound.numerator))/log(2)
            self.assertAlmostEqual(cap_counts.bad_route_margin(q, h, -log(float(Fraction(x)))), margin, places=8)

    def test_invalid_geometry_after_cached_valid_call(self):
        cap_counts.one_region_moment(1, '1/2', cap=1, slots=6, windows=3)
        for q, cap in ((True, 1), (1, True), (-1, 1), (7, 1), (1, 0), (1, 4)):
            with self.assertRaises(ValueError):
                cap_counts.one_region_moment(q, '1/2', cap=cap, slots=6, windows=3)


if __name__ == '__main__':
    unittest.main()
