"""Small algebra, exact restricted-field checks, and mocked inventories only."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from random import Random
import unittest
from unittest.mock import patch

import hd_lambda_extra_transvection as target
from hd_quaternary_outer_relations import multiply_gf4


def rank_gf4(rows, columns):
    rows = [list(row) for row in rows]
    pivots = 0
    for j in range(columns):
        i = next((i for i in range(pivots, len(rows)) if rows[i][j]), None)
        if i is None:
            continue
        rows[pivots], rows[i] = rows[i], rows[pivots]
        inverse = next(c for c in (1, 2, 3) if multiply_gf4(rows[pivots][j], c) == 1)
        rows[pivots] = [multiply_gf4(x, inverse) for x in rows[pivots]]
        for i in range(len(rows)):
            if i != pivots and rows[i][j]:
                coefficient = rows[i][j]
                rows[i] = [x ^ multiply_gf4(coefficient, y) for x, y in zip(rows[i], rows[pivots])]
        pivots += 1
    return pivots


def brute_counts(binary_rows, width):
    counts = [0]*(1 << width)
    for values in product(range(4), repeat=width):
        valid = True
        for row in binary_rows:
            syndrome = 0
            for j, x in enumerate(values):
                syndrome ^= multiply_gf4(2 ^ ((row >> j) & 1), x)
            valid &= syndrome == 0
        if valid:
            counts[sum(bool(x) << i for i, x in enumerate(values))] += 1
    return tuple(counts)


class ExtraTransvectionTests(unittest.TestCase):
    def test_new_matrix_is_involution_and_all_entries_nonzero(self):
        for q, lam in ((4, 2), (4, 3), (256, 0xbc)):
            for support in ((0, 1, 2, 3), (0, 1, 4, 5), (3, 6, 9, 12)):
                h = target.matrix(lam, q, support=support)
                self.assertEqual(h, tuple(zip(*h)))
                self.assertEqual(set(x for row in h for x in row), {lam, lam ^ 1})
                for i in range(16):
                    e = tuple(int(j == i) for j in range(16))
                    self.assertEqual(target.direct(target.direct(e, lam, q, support=support),
                                                   lam, q, support=support), e)
                v = tuple(int(i in support) for i in range(16))
                self.assertEqual(target.original.direct(v, lam, q), v)
                self.assertEqual(target.direct(v, lam, q, support=support), v)

    def test_default_fused_naive_and_direct_circuits_agree(self):
        rng = Random(100934)
        for q, lam in ((4, 2), (256, 0xbc)):
            vectors = [tuple(int(i == j) for i in range(16)) for j in range(16)]
            vectors += [tuple(rng.randrange(q) for _ in range(16)) for _ in range(12)]
            for x in vectors:
                direct = target.direct(x, lam, q)
                self.assertEqual(target.circuit(x, lam, q), direct)
                self.assertEqual(target.circuit(x, lam, q, fuse_default=False), direct)
                old = target.original.circuit(x, lam, q)
                parity = old[0] ^ old[1] ^ old[2] ^ old[3]
                self.assertEqual(parity, x[0] ^ x[1] ^ x[2] ^ x[3])
                self.assertEqual(tuple(y ^ (parity if i < 4 else 0) for i, y in enumerate(old)), direct)
        x = tuple(rng.randrange(256) for _ in range(16))
        support = (0, 1, 4, 5)
        self.assertEqual(target.circuit(x, support=support), target.direct(x, support=support))

    def test_exact_primitive_counts_are_not_native_claims(self):
        fused = target.primitive_counts()
        naive = target.primitive_counts(fuse_default=False)
        self.assertEqual((fused['binary_xors'], fused['fixed_products'], fused['extra_xors_over_Hlambda']), (61, 1, 4))
        self.assertEqual((naive['binary_xors'], naive['fixed_products'], naive['extra_xors_over_Hlambda']), (64, 1, 7))
        self.assertEqual(target.primitive_counts(support=(0, 1, 4, 5))['extra_xors_over_Hlambda'], 7)
        self.assertFalse(fused['native_performance_claim'])
        self.assertEqual(fused['construction'], target.CONSTRUCTION)

    def test_binary_parity_and_small_exact_support_moments(self):
        for complement in ((), (0,), (1, 4), tuple(range(5)), tuple(range(7))):
            full = target.parity_binary_rows(complement)
            h = target.matrix(2, 4)
            for o, row in zip(complement, full):
                self.assertEqual(tuple(2 ^ ((row >> i) & 1) for i in range(16)), h[o])
            # Restrict coordinates so all field vectors can be enumerated.
            width = 4
            restricted = tuple(row & ((1 << width)-1) for row in full)
            data = target.middle.support_counts_from_binary(restricted, width, q=4)
            self.assertEqual(data['counts'], brute_counts(restricted, width))
            for mask, rank in enumerate(data['subset_ranks']):
                columns = tuple(i for i in range(width) if mask >> i & 1)
                rows = tuple(tuple(2 ^ ((row >> i) & 1) for i in columns) for row in restricted)
                self.assertEqual(rank, rank_gf4(rows, len(columns)))
            moments = target.symmetry.exact_moments(data['counts'], data['codimension'], orders=(2, 3, 4), q=4)
            exact = {p: F(0) for p in moments}
            for y in product(range(4), repeat=width):
                mask = sum(bool(v) << i for i, v in enumerate(y))
                likelihood = F(4**data['codimension']*data['counts'][mask], 3**mask.bit_count())
                for p in exact:
                    exact[p] += likelihood**p/4**width
            self.assertEqual(moments, exact)

    def test_only_actual_stabilizer_preserves_variant(self):
        maps = target.stabilizer_maps()
        self.assertEqual(len(maps), 4)
        self.assertTrue(set(maps) < set(target.symmetry.orthogonal_maps()))
        h = target.matrix()
        for image in maps:
            self.assertEqual({image[i] for i in target.DEFAULT_SUPPORT}, set(target.DEFAULT_SUPPORT))
            self.assertTrue(all(h[image[i]][image[j]] == h[i][j] for i, j in product(range(16), repeat=2)))
        for left, right in product(maps, repeat=2):
            self.assertIn(tuple(left[right[i]] for i in range(16)), maps)
        excluded = next(image for image in target.symmetry.orthogonal_maps() if image not in maps)
        self.assertTrue(any(h[excluded[i]][excluded[j]] != h[i][j] for i, j in product(range(16), repeat=2)))

    def test_all_labelled_complements_retained_and_disjoint_parity_unchanged(self):
        for size in range(6):
            seen = set()
            for orbit in target.complement_orbits(size):
                for complement, image in orbit['members']:
                    self.assertNotIn(complement, seen)
                    seen.add(complement)
                    self.assertEqual(tuple(sorted(image[i] for i in orbit['representative'])), complement)
                    self.assertIn(image, target.stabilizer_maps())
            self.assertEqual(seen, set(combinations(range(16), size)))
        self.assertEqual(comb(12, 5), 792)
        for complement in combinations(range(4, 16), 5):
            self.assertEqual(target.parity_binary_rows(complement), target.middle.parity_binary_rows(complement))

    def test_mock_inventory_uses_new_subgroup_and_construction(self):
        calls = []
        def fake(complement, *, orders, support):
            calls.append(complement)
            return dict(moments={p: F(p) for p in orders})
        with patch.object(target, 'support_moments', side_effect=fake):
            data = target.inventory(complement_sizes=(0, 1, 2), max_seconds=10)
        self.assertTrue(data['complete'])
        self.assertEqual(data['subgroup_size'], 4)
        self.assertEqual(data['old_subgroup_size'], 48)
        self.assertEqual(data['construction'], target.CONSTRUCTION)
        self.assertEqual(len(data['rows']), 137)
        self.assertEqual(len(calls), sum(len(target.complement_orbits(b)) for b in range(3)))
        self.assertEqual(data['maximum_moments_by_arity']['14']['4'], '4')
        self.assertFalse(data['old_receipts_transferred'])
        self.assertEqual(data['regional_provider_runs'], 0)
        self.assertFalse(data['native_performance_claim'])

    def test_guards_timeout_and_pins(self):
        self.assertEqual(target.checked_dependencies()['hd_lambda_outer.py'], target.PINS['hd_lambda_outer.py'])
        for support in ((0, 1, 2), (0, 1, 2, 4), (0, 0, 1, 1), (0, 1, 2, 16), (False, 1, 2, 3)):
            with self.assertRaises(ValueError):
                target.matrix(support=support)
        for complement in (tuple(range(8)), (0, 0), (16,), (True,)):
            with self.assertRaises(ValueError):
                target.parity_binary_rows(complement)
        for kwargs in ({'complement_sizes': (6,)}, {'complement_sizes': (1, 1)}, {'max_seconds': float('nan')}):
            with self.assertRaises(ValueError):
                target.inventory(**kwargs)
        with patch.object(target, 'perf_counter', side_effect=(0., 2., 2.)), \
                patch.object(target, 'support_moments') as calculate:
            data = target.inventory(complement_sizes=(1,), max_seconds=1)
        calculate.assert_not_called()
        self.assertFalse(data['complete'])
        self.assertEqual(data['rows'], [])


if __name__ == '__main__':
    unittest.main()
