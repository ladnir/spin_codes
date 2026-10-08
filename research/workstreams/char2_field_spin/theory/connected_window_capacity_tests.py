"""Independent small checks of merged-window counts and interval replay."""

from collections import Counter
from fractions import Fraction as F
from itertools import product
import unittest

from connected_window_capacity import (evaluate_family, families, rational_moment,
                                       shape_count)
from subfield_certificate import Interval


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    elif parts > 1:
        for first in range(1, total-parts+2):
            for rest in compositions(total-first, parts-1):
                yield (first,)+rest


def actual_window_counts(regions, units, budget):
    choices = [None]+[(start, end) for start in range(units)
                     for end in range(start+1, units+1)]
    counts = Counter()
    for selected in product(choices, repeat=regions):
        used = [(region, item) for region, item in enumerate(selected) if item]
        if not used or sum(end-start for _, (start, end) in used) > budget:
            continue
        lengths = []
        previous = None
        for region, (start, end) in used:
            if previous and previous[0] == region-1 and previous[1] == units and start == 0:
                lengths[-1] += 1
            else:
                lengths.append(1)
            previous = (region, end)
        tail = 0
        if previous == (regions-1, units):
            tail = 1 if lengths[-1] == 1 else 2
        counts[(len(used), len(lengths), lengths.count(1), tail)] += 1
    return counts


class ConnectedCapacityTests(unittest.TestCase):
    def test_component_shape_formula(self):
        for used in range(1, 9):
            for count in range(1, used+1):
                shapes = list(compositions(used, count))
                for singles in range(count+1):
                    for tail in range(3):
                        actual = sum(shape.count(1) == singles and
                                     (tail == 0 or (shape[-1] == 1) == (tail == 1))
                                     for shape in shapes)
                        self.assertEqual(shape_count(used, count, singles, tail), actual)

    def test_exhaustive_small_geometry_count(self):
        for regions in range(1, 5):
            for units in (1, 2, 3):
                for budget in range(1, min(regions*units, 6)+1):
                    cover = {(entry['used'], entry['components'], entry['singles'],
                              entry['tail_kind']): entry['count']
                             for entry in families(regions, 0, units, budget)}
                    for key, actual in actual_window_counts(regions, units, budget).items():
                        self.assertIn(key, cover)
                        self.assertGreaterEqual(cover[key], actual)

    def test_interval_moment_contains_exact(self):
        for probability in (F(1, 7), F(1, 2), F(1)):
            for mark in (F(1), F(4, 3), F(13, 4)):
                for parity in (-1, 0, 2, 4):
                    exact = rational_moment(4, parity, probability, mark)**11/mark**17
                    pi = Interval.rational(probability.numerator, probability.denominator)
                    zi = Interval.rational(mark.numerator, mark.denominator)
                    interval = rational_moment(4, parity, pi, zi)**11/zi**17
                    self.assertLessEqual(F(interval.lo), exact)
                    self.assertGreaterEqual(F(interval.hi), exact)

    def test_directed_family_replay(self):
        for family in families(5, 1, 3, 5):
            entry, upper = evaluate_family(family, 3, 1, 2, 3, True)
            probability, mark = F(*entry['probability']), F(*entry['mark'])
            exact = rational_moment(family['variables'], 1-family['full'],
                                    probability, mark)**3/mark**entry['crossing']
            exact = min(F(1), exact)*family['count']
            self.assertGreaterEqual(F(upper.hi), exact)

    def test_joining_and_tail_reduce_capacity(self):
        entries = list(families(4, 1, 8, 12))
        chain = next(e for e in entries if e['used'] == 2 and e['components'] == 1
                     and e['tail_kind'] == 0)
        tail = next(e for e in entries if e['used'] == 2 and e['components'] == 1
                    and e['tail_kind'] == 2)
        self.assertEqual(chain['full'], 0)
        self.assertEqual(tail['full'], 1)
        self.assertEqual(evaluate_family(chain, 8, 1, 4, 8)[0]['crossing'], 5)
        self.assertEqual(evaluate_family(tail, 8, 1, 4, 8)[0]['crossing'], 1)


if __name__ == '__main__':
    unittest.main()
