"""Exact checks separating initialization scores from true source dimension."""

import unittest

import numpy as np

from repaired_geometry_windows import (GEOMETRIES, anchored_starts, inverse_positions,
                                       initialize, parity_ladder, screen)
from window_dimension_search import generate_route, count_hits, dimension, strict_window_geometry


class RepairedWindowTests(unittest.TestCase):
    def test_centered_alignment_and_guaranteed_positive_dimension(self):
        route = generate_route(32, 8, 519)
        inverse = inverse_positions(route)
        for row in (0, 7, 31):
            starts = anchored_starts(inverse, (row,), 8, 4)
            self.assertTrue(np.all(starts % 4 == 0))
            self.assertTrue(np.all((0 <= starts) & (starts+8 <= 32)))
            counts = count_hits(route, starts, 8)
            self.assertEqual(int(counts[row]), 8)
            self.assertGreaterEqual(dimension(counts, 5), 3)

    def test_multi_anchor_choices_are_exact_weighted_maxima(self):
        route = generate_route(32, 8, 914)
        inverse = inverse_positions(route)
        for anchors in ((0, 31), (0, 7, 31), (0, 7, 19, 31)):
            starts = anchored_starts(inverse, anchors, 8, 4)
            weights = [1 << (len(anchors)-1-i) for i in range(len(anchors))]
            for j, first in enumerate(starts):
                def score(start):
                    return sum(w for row, w in zip(anchors, weights) if start <= inverse[j, row] < start+8)
                self.assertEqual(score(first), max(score(start) for start in range(0, 25, 4)))
            self.assertEqual(int(count_hits(route, starts, 8)[anchors[0]]), 8)

    def test_real_geometry_strict_budgets_and_padding(self):
        expected = {'33/22': (4752, 157343, 18), '64/44': (2384, 152780, 41)}
        for name, geometry in GEOMETRIES.items():
            rows, regions, parity = (geometry[k] for k in ('rows', 'regions', 'parity'))
            width, threshold = strict_window_geometry(rows, regions)
            k = regions-parity
            padding_rows = (rows*k-2**20+k-1)//k
            self.assertEqual((width, threshold, padding_rows), expected[name])
            self.assertLessEqual(width*regions, threshold)
            self.assertLessEqual((rows-padding_rows)*k, 2**20)

    def test_ladder_reports_true_objective_separately(self):
        route = generate_route(24, 6, 132)
        initial = anchored_starts(inverse_positions(route), (5, 13), 6, 2)
        initial_true = dimension(count_hits(route, initial, 6), 3)
        result, stages = parity_ladder(route, 6, initial, 3, (0, 1, 2, 3),
                                       alignment=2, max_sweeps=8, seed=191)
        self.assertEqual(result.dimension, dimension(count_hits(route, result.starts, result.widths), 3))
        self.assertGreaterEqual(result.dimension, initial_true)
        self.assertTrue(all(stage['true_parity'] == 3 for stage in stages))
        self.assertTrue(all(stage['surrogate_dimension'] >= stage['true_dimension'] for stage in stages))
        self.assertEqual(stages[0]['surrogate_dimension'], 36)

    def test_64_region_initialization_avoids_zero_plateau(self):
        route = generate_route(256, 64, 21)
        starts = initialize(route, 16, 901, 4)
        self.assertEqual(len(starts), 8)
        for _, window_starts, anchors in starts:
            if anchors:
                counts = count_hits(route, window_starts, 16)
                self.assertEqual(int(counts[anchors[0]]), 64)
                self.assertGreaterEqual(dimension(counts, 20), 44)

    def test_small_full_screen_is_replayable(self):
        geometry = {'rows': 80, 'regions': 6, 'parity': 2}
        parameters = dict(exchanges=8, max_sweeps=4, alignment=2, requested_source=300)
        first = screen(geometry, 33, **parameters)
        repeated = screen(geometry, 33, **parameters)
        self.assertEqual(first, repeated)
        self.assertGreaterEqual(first['max_true_dimension_found'], first['true_dimension_before_exchange'])
        self.assertLessEqual(first['total_window_symbols'], first['bad_output_limit'])
        self.assertFalse(first['probability_claimed'])
        self.assertTrue(first['no_find_is_not_certificate'])

    def test_anchor_guards(self):
        inverse = inverse_positions(generate_route(16, 4, 1))
        for anchors, width, alignment in (((), 4, 2), ((1, 1), 4, 2), ((1,), 4, 0), ((1,), 3, 2)):
            with self.assertRaises(ValueError):
                anchored_starts(inverse, anchors, width, alignment)


if __name__ == '__main__':
    unittest.main()
