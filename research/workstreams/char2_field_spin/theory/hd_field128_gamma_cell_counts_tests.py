"""Tiny exact combinatorial oracles only; no full-size distance gates."""

from itertools import product
import unittest

import hd_field128_gamma_cell_counts as counts
from hd_shifted_pairs_count import exact_profiles, one_pairing_count, circular_run_count, profile_count_upper


class Budget:
    def __init__(self, limit=None):
        self.calls, self.limit = 0, limit

    def check(self):
        self.calls += 1
        if self.limit is not None and self.calls > self.limit:
            raise TimeoutError('mock count deadline')


def actual(cell, profiles):
    if cell.kind == 'interior':
        return sum(count for point, count in profiles.items() if cell.contains(point))
    even = cell.kind == 'even_face'
    total = 0
    for (r, a, b), count in profiles.items():
        if r % 2 == 0 and (a if even else b) == r//2:
            point = r//2, r-a-b
            if cell.contains(point):
                total += count
    return total


class CellCountTests(unittest.TestCase):
    def test_prefix_recurrences_match_independent_binomial_oracles(self):
        cache, budget = counts.PrefixCache(2), Budget()
        for pairs in range(1, 11):
            for weight in range(2*pairs+1):
                marginal, cycle = cache.get('interior', pairs, weight, budget)
                for full in range(-1, pairs+2):
                    self.assertEqual(marginal.interval(full, full), one_pairing_count(pairs, weight, full))
                for runs in range(-1, 2*pairs+2):
                    self.assertEqual(cycle.interval(runs, runs), circular_run_count(2*pairs, weight, runs))
                self.assertEqual(marginal.interval(-100, 100), cycle.interval(-100, 100))
                self.assertEqual(marginal.interval(3, 2), 0)
                self.assertLessEqual(len(cache), 2)

    def test_single_profile_cells_match_existing_interior_profile_upper(self):
        budget, cache = Budget(), counts.PrefixCache()
        for pairs in range(2, 7):
            for r in range(1, 2*pairs):
                for a in range(max(0, r-pairs), (r-1)//2+1):
                    for b in range(max(0, r-pairs), (r-1)//2+1):
                        cell = counts.Cell('interior', ((r, r), (a, a), (b, b)), pairs)
                        self.assertEqual(counts.support_count(cell, budget, cache=cache),
                                         profile_count_upper(pairs, r, a, b))

    def test_interior_projection_sums_bound_actual_supports_without_profile_factor(self):
        budget, cache = Budget(), counts.PrefixCache()
        for pairs in range(2, 7):
            profiles = exact_profiles(pairs)
            boxes = (((1, 2*pairs-1), (0, pairs), (0, pairs)),
                     ((1, pairs), (0, 1), (0, pairs)),
                     ((pairs, 2*pairs-1), (1, pairs), (0, pairs)),
                     ((2, 2*pairs-2), (0, 0), (1, pairs)),
                     ((pairs, pairs), (0, pairs), (0, pairs)))
            for ranges in boxes:
                cell = counts.Cell('interior', ranges, pairs)
                result = counts.count_diagnostics(cell, budget, cache=cache)
                expected = [0, 0, 0]
                # Tiny independent lattice projection sets, never used by production.
                for r in range(1, 2*pairs):
                    points = [(a, b) for a, b in product(range(pairs+1), repeat=2)
                              if cell.contains((r, a, b))]
                    expected[0] += sum(one_pairing_count(pairs, r, a) for a in {a for a, b in points})
                    expected[1] += sum(one_pairing_count(pairs, r, b) for b in {b for a, b in points})
                    expected[2] += sum(circular_run_count(2*pairs, r, k) for k in {r-a-b for a, b in points})
                self.assertEqual(list(result['projection_sums'].values()), expected)
                self.assertEqual(result['support_count_upper'], min(expected))
                self.assertGreaterEqual(result['support_count_upper'], actual(cell, profiles))
                self.assertEqual(result['profile_pairs_enumerated'], 0)
                self.assertFalse(result['exact_support_count'])

    def test_both_faces_are_exact_contracted_counts_and_split_additively(self):
        budget, cache = Budget(), counts.PrefixCache()
        for pairs in range(2, 8):
            profiles = exact_profiles(pairs)
            for kind in ('even_face', 'odd_face'):
                for ranges in (((1, pairs-1), (1, pairs)),
                               ((1, max(1, pairs//2)), (1, 1)),
                               ((max(1, pairs//2), pairs-1), (2, pairs))):
                    cell = counts.Cell(kind, ranges, pairs)
                    result = counts.count_diagnostics(cell, budget, cache=cache)
                    self.assertTrue(result['exact_support_count'])
                    self.assertEqual(result['support_count_upper'], actual(cell, profiles))
                    expected = sum(circular_run_count(pairs, u, k)
                                   for u in range(1, pairs) for k in range(1, pairs)
                                   if cell.contains((u, k)))
                    self.assertEqual(result['support_count_upper'], expected)
                if pairs > 2:
                    parent = counts.Cell(kind, ((1, pairs-1), (1, pairs)), pairs)
                    left, right = parent.split(0, pairs//2)
                    self.assertEqual(counts.support_count(parent, budget, cache=cache),
                                     counts.support_count(left, budget, cache=cache)
                                     +counts.support_count(right, budget, cache=cache))

    def test_cache_is_bounded_reuses_completed_tables_and_never_caches_timeout(self):
        cache, budget = counts.PrefixCache(2), Budget()
        first = cache.get('interior', 5, 2, budget)
        self.assertIs(first, cache.get('interior', 5, 2, budget))
        cache.get('face', 5, 2, budget)
        cache.get('interior', 5, 3, budget)
        self.assertEqual(len(cache), 2)
        self.assertIsNot(first, cache.get('interior', 5, 2, budget))
        cache.clear()
        with self.assertRaises(TimeoutError):
            cache.get('interior', 5, 2, Budget(1))
        self.assertEqual(len(cache), 0)
        cache.get('interior', 5, 2, budget)
        with self.assertRaises(TimeoutError):
            cache.get('interior', 5, 2, Budget(0))
        cell = counts.Cell('interior', ((1, 9), (0, 5), (0, 5)), 5)
        with self.assertRaises(TimeoutError):
            counts.support_count(cell, Budget(1), cache=cache)

    def test_empty_domains_guards_and_source_closure(self):
        budget = Budget()
        for cell in (counts.Cell('interior', ((1, 1), (1, 1), (0, 0)), 4),
                     counts.Cell('even_face', ((1, 1), (2, 2)), 4),
                     counts.Cell('odd_face', ((3, 3), (2, 3)), 4)):
            self.assertEqual(counts.support_count(cell, budget), 0)
        with self.assertRaises(ValueError):
            counts.support_count(object(), budget)
        for size in (0, 17, True):
            with self.assertRaises(ValueError):
                counts.PrefixCache(size)
        for kind, pairs, weight in (('unknown', 4, 1), ('face', 4, 5), ('interior', True, 1)):
            with self.assertRaises(ValueError):
                counts.PrefixCache().get(kind, pairs, weight, budget)
        pins = counts.source_pins()
        self.assertIn('hd_field128_gamma_cell_counts.py', pins)
        self.assertIn('hd_field128_gamma_cell_geometry.py', pins)


if __name__ == '__main__':
    unittest.main()
