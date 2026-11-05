from fractions import Fraction as F
from itertools import combinations, product
from random import Random
from unittest.mock import patch
import unittest

import hd_lambda_renyi_middle as middle
import hd_lambda_renyi_orbits as orbit
import hd_lambda_renyi_support as frozen
from hd_lambda_renyi_support_tests import brute, rank4


class MiddleRenyiTests(unittest.TestCase):
    def test_frozen_small_checks_agree(self):
        rng = Random(4735)
        for b in range(5):
            for width in (1, 4, 8):
                for _ in range(3):
                    rows = tuple(rng.randrange(1 << width) for _ in range(b))
                    self.assertEqual(middle.subset_ranks(rows, width), frozen.subset_ranks(rows, width))
                    self.assertEqual(middle.support_counts_from_binary(rows, width, q=4),
                                     frozen.support_counts_from_binary(rows, width, q=4))

    def test_five_to_seven_checks_direct_GF4(self):
        rng = Random(59833)
        for b in (5, 6, 7):
            for width in (3, 6):
                rows = tuple(rng.randrange(1 << width) for _ in range(b))
                ranks = middle.subset_ranks(rows, width)
                for subset, actual in enumerate(ranks):
                    columns = tuple(i for i in range(width) if subset >> i & 1)
                    matrix = tuple(tuple(2 ^ ((row >> i) & 1) for i in columns) for row in rows)
                    self.assertEqual(actual, rank4(matrix, len(columns)))
                if width == 3:
                    self.assertEqual(middle.support_counts_from_binary(rows, width, q=4)['counts'], brute(rows, width))

    def test_full_parity_restriction_and_identity(self):
        for complement in (tuple(range(5)), (0, 1, 2, 4, 8, 15), tuple(range(7))):
            rows = middle.parity_binary_rows(complement)
            for o, row in zip(complement, rows):
                self.assertEqual(row, sum((int(o == i) ^ ((o & i).bit_count() & 1)) << i for i in range(16)))
            restricted = tuple(row & 15 for row in rows)
            data = middle.support_counts_from_binary(restricted, 4, q=4)
            self.assertEqual(data['counts'], brute(restricted, 4))
            moments = orbit.exact_moments(data['counts'], data['codimension'], orders=(2, 3, 4), q=4)
            for p in (2, 3, 4):
                self.assertEqual(moments[p], frozen.density_moment_from_counts(data['counts'], data['codimension'], order=p, q=4))

    def test_full_rank_shortcut_and_wide_basis(self):
        # Seven independent pivot slots occupy49 bits, beyond byte/state IDs.
        rows = tuple((1 << i) | (1 << (i+7)) for i in range(7))
        ranks = middle.subset_ranks(rows, 14)
        for subset in (0, 1, 127, 128, 255, 0x1555, 0x2aaa, 0x3fff):
            columns = tuple(i for i in range(14) if subset >> i & 1)
            matrix = tuple(tuple(2 ^ ((row >> i) & 1) for i in columns) for row in rows)
            self.assertEqual(ranks[subset], rank4(matrix, len(columns)))
        self.assertEqual(ranks[-1], 7)

    def test_all_labelled_orbits_with_witness(self):
        for b in (5, 6, 7):
            seen = set()
            for data in middle.complement_orbits(b):
                representative = data['representative']
                for complement, image in data['members']:
                    self.assertNotIn(complement, seen)
                    seen.add(complement)
                    self.assertEqual(tuple(sorted(image[i] for i in representative)), complement)
                    self.assertIn(image, orbit.orthogonal_maps())
            self.assertEqual(seen, set(combinations(range(16), b)))

    def test_mock_inventory_one_table_per_orbit(self):
        calls = []
        def fake(complement, *, orders):
            calls.append(complement)
            return dict(moments={p: F(p) for p in orders})
        with patch.object(middle, 'representative_moments', side_effect=fake):
            data = middle.inventory(max_seconds=10)
        self.assertTrue(data['complete'])
        self.assertEqual(data['complement_sizes'], [5])
        self.assertEqual(data['moment_orders'], [2, 3, 4])
        self.assertEqual(len(calls), len(middle.complement_orbits(5)))
        self.assertEqual(len(data['rows']), 4368)
        self.assertEqual(data['maximum_moments_by_arity']['11'], {'2': '2', '3': '3', '4': '4'})
        self.assertTrue(data['source_pins_unchanged'])
        self.assertFalse(data['whole_code_certificate'])

    def test_guards_and_timeout(self):
        for rows, width in (((1,)*8, 4), ((16,), 4), ((True,), 4), ((), 17)):
            with self.assertRaises(ValueError):
                middle.subset_ranks(rows, width)
        for complement in (tuple(range(4)), tuple(range(8)), (0, 1, 2, 3, 3), (0, 1, 2, 3, True)):
            with self.assertRaises(ValueError):
                middle.parity_binary_rows(complement)
        for kwargs in ({'complement_sizes': (6,)}, {'complement_sizes': (5, 5)},
                       {'max_complement_size': 8}, {'max_seconds': float('inf')}):
            with self.assertRaises(ValueError):
                middle.inventory(**kwargs)
        with patch.object(middle, 'perf_counter', side_effect=(0., 2., 2.)), \
                patch.object(middle, 'representative_moments') as call:
            data = middle.inventory(max_seconds=1)
        self.assertFalse(data['complete'])
        call.assert_not_called()


if __name__ == '__main__':
    unittest.main()
