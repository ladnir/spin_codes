from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from unittest.mock import patch
import unittest

import hd_lambda_renyi_orbits as orbit
import hd_lambda_renyi_support as support


class RenyiOrbitTests(unittest.TestCase):
    def test_explicit_symbolic_maps_and_group_closure(self):
        maps = orbit.orthogonal_maps()
        self.assertEqual(len(set(maps)), 48)
        for image in maps:
            self.assertEqual(sorted(image), list(range(16)))
            for i, j in product(range(16), repeat=2):
                old = int(i == j) ^ ((i & j).bit_count() & 1)
                moved = int(image[i] == image[j]) ^ ((image[i] & image[j]).bit_count() & 1)
                self.assertEqual(old, moved)
        known = set(maps)
        for left, right in product(maps, repeat=2):
            self.assertIn(tuple(left[right[i]] for i in range(16)), known)

    def test_all_complements_have_actual_map_once(self):
        for b in range(5):
            seen = set()
            for data in orbit.complement_orbits(b):
                representative = data['representative']
                self.assertEqual(representative, min(member for member, _ in data['members']))
                for target, image in data['members']:
                    self.assertNotIn(target, seen)
                    seen.add(target)
                    self.assertEqual(tuple(sorted(image[i] for i in representative)), target)
                    self.assertIn(image, orbit.orthogonal_maps())
            self.assertEqual(seen, set(combinations(range(16), b)))
            self.assertEqual(len(seen), comb(16, b))

    def test_tiny_restricted_code_permutation_equivalence(self):
        # Six physical coordinates keep exact counts tiny. Move both parity
        # labels and coordinate labels, then reverse their local ordering.
        source, coordinates = (0, 1, 3), (0, 1, 2, 4, 8, 15)
        rows = support.parity_binary_rows(source)
        local = tuple(sum(((row >> col) & 1) << j for j, col in enumerate(coordinates)) for row in rows)
        original = support.support_counts_from_binary(local, 6, q=4)
        for image in orbit.orthogonal_maps()[::7]:
            moved_rows = support.parity_binary_rows(tuple(image[i] for i in source))
            moved_cols = tuple(image[i] for i in coordinates)[::-1]
            packed = tuple(sum(((row >> col) & 1) << j for j, col in enumerate(moved_cols)) for row in moved_rows)
            moved = support.support_counts_from_binary(packed, 6, q=4)
            for s, count in enumerate(original['counts']):
                target = sum(((s >> i) & 1) << (5-i) for i in range(6))
                self.assertEqual(moved['counts'][target], count)
            self.assertEqual(orbit.exact_moments(original['counts'], original['codimension'], q=4),
                             orbit.exact_moments(moved['counts'], moved['codimension'], q=4))

    def test_generic_moments_direct_uniform_space(self):
        data = support.support_counts_from_binary((0b1011, 0b0101), 4, q=4)
        orders = (1, 2, 3, 4, 5, 8, 16)
        actual = orbit.exact_moments(data['counts'], data['codimension'], orders=orders, q=4)
        expected = {p: F(0) for p in orders}
        for y in product(range(4), repeat=4):
            mask = sum(bool(v) << i for i, v in enumerate(y))
            likelihood = F(4**data['codimension']*data['counts'][mask], 3**mask.bit_count())
            for p in orders:
                expected[p] += likelihood**p/4**4
        self.assertEqual(actual, expected)
        for p in range(1, 5):
            self.assertEqual(actual[p], support.density_moment_from_counts(data['counts'], data['codimension'], q=4, order=p))

    def test_inventory_reuses_once_and_retains_labelled_rows(self):
        calls = []
        def fake(complement, *, orders):
            calls.append(complement)
            return dict(moments={p: F(2+p, 2) for p in orders})
        with patch.object(orbit, 'representative_moments', side_effect=fake):
            data = orbit.inventory(complement_sizes=(0, 1, 2), orders=(3, 4), max_seconds=10)
        self.assertTrue(data['complete'])
        self.assertTrue(data['source_pins_unchanged'])
        self.assertEqual(len(calls), sum(len(orbit.complement_orbits(b)) for b in range(3)))
        self.assertEqual(len(data['rows']), 1+16+120)
        self.assertEqual(data['actual_support_tables'], len(calls))
        for row in data['rows']:
            p = row['representative_to_complement_permutation']
            self.assertEqual(sorted(p[i] for i in row['representative']), row['complement'])
            self.assertEqual(row['moments'], {'3': '5/2', '4': '3'})
        self.assertEqual(data['regional_provider_runs'], 0)
        self.assertFalse(data['whole_code_certificate'])

    def test_timeout_is_partial_without_rep_call(self):
        with patch.object(orbit, 'perf_counter', side_effect=(0., 2., 2.)), \
                patch.object(orbit, 'representative_moments') as call:
            data = orbit.inventory(complement_sizes=(4,), max_seconds=1)
        self.assertFalse(data['complete'])
        self.assertEqual(data['rows'], [])
        call.assert_not_called()

    def test_guards(self):
        for b in (-1, 5, True):
            with self.assertRaises(ValueError):
                orbit.complement_orbits(b)
        for orders in ((), (3, 3), (17,), (True,)):
            with self.assertRaises(ValueError):
                orbit.exact_moments((1, 3), 0, orders=orders, q=4)
        for sizes in ((), (4, 4), (5,), (True,)):
            with self.assertRaises(ValueError):
                orbit.inventory(complement_sizes=sizes, max_seconds=1)
        with self.assertRaises(ValueError):
            orbit.exact_moments((1, 2), 0, q=4)


if __name__ == '__main__':
    unittest.main()
