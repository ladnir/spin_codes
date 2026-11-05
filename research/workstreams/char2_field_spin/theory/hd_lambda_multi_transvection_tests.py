"""Tiny algebra and equivalence tests; no all-complement or candidate screen."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import random
import unittest
from unittest.mock import patch

import hd_lambda_multi_transvection as target


def apply(matrix, x, q):
    return tuple(_xor(target.original._multiply(a, b, q) for a, b in zip(row, x)) for row in matrix)


def _xor(values):
    result = 0
    for x in values:
        result ^= x
    return result


class MultiTransvectionTests(unittest.TestCase):
    def test_supports_and_independence_modulo_base(self):
        self.assertEqual(len(target.quiet_supports()), 140)
        candidates = target.concrete_candidates()
        self.assertEqual([target.syndrome_extension(s)['extended_binary_span_rank'] for s in candidates], [6, 7, 8, 8])
        self.assertEqual([target.syndrome_extension(s)['independent_modulo_RM1'] for s in candidates], [False, True, True, True])
        with self.assertRaises(ValueError):
            target.configuration(((0, 1, 2, 3), (0, 4, 8, 12)))
        for bad in (((0, 1, 2, 4),), ((0, 1, 2, 3),)*2, ((True, 1, 2, 3),)):
            with self.assertRaises(ValueError):
                target.configuration(bad)

    def test_matrix_circuit_involution_and_all_entries_nonzero(self):
        rng = random.Random(442)
        for supports in target.concrete_candidates():
            for q, lam in ((4, 2), (256, 188)):
                h = target.matrix(supports, lam, q)
                self.assertEqual(h, tuple(zip(*h)))
                self.assertEqual(set(x for row in h for x in row), {lam, lam ^ 1})
                vectors = [tuple(int(i == j) for i in range(16)) for j in range(16)]
                vectors += [tuple(rng.randrange(q) for _ in range(16)) for _ in range(2)]
                for x in vectors:
                    y = target.circuit(x, supports, lam, q)
                    self.assertEqual(y, apply(h, x, q))
                    self.assertEqual(target.circuit(y, supports, lam, q), x)

    def test_explicit_shared_circuit_costs(self):
        counts = [target.primitive_counts(s) for s in target.concrete_candidates()]
        self.assertEqual([c['extra_xors_over_Hlambda'] for c in counts], [8, 9, 14, 16])
        self.assertEqual([c['fixed_products'] for c in counts], [1]*4)
        self.assertTrue(all(not c['native_performance_claim'] for c in counts))

    def test_projective_binary_criterion_vs_actual_field_normalization(self):
        for q, lam in ((4, 2), (256, 188), (256, 2)):
            for height in range(1, 6):
                rays = []
                for pattern in range(1 << height):
                    vector = tuple(lam ^ ((pattern >> i) & 1) for i in range(height))
                    inv = target.original._inverse(vector[0], q)
                    rays.append(tuple(target.original._multiply(x, inv, q) for x in vector))
                for a, b in product(range(1 << height), repeat=2):
                    self.assertEqual(rays[a] == rays[b], target.projective_key(a, height) == target.projective_key(b, height))

    def test_rank_one_moment_part_against_exhaustive_toy_code(self):
        q, lam, width = 4, 2, 4
        for binary_rows in ((0b0101,), (0b0101, 0b0011), (0b0101, 0b0011, 0b0001)):
            h = tuple(tuple(lam ^ ((row >> j) & 1) for j in range(width)) for row in binary_rows)
            counts = Counter()
            for x in product(range(q), repeat=width):
                if not any(apply(h, x, q)):
                    counts[sum(bool(v) << i for i, v in enumerate(x))] += 1
            b = len(binary_rows)
            self.assertEqual(sum(counts.values()), q**(width-b))
            classes = target.projective_classes(binary_rows, width=width)
            sizes = tuple(sorted(len(c) for _, c in classes))
            rank_one_masks = set()
            for _, columns in classes:
                for size in range(1, len(columns)+1):
                    for subset in combinations(columns, size):
                        mask = sum(1 << i for i in subset)
                        rank_one_masks.add(mask)
                        self.assertEqual(counts[mask], target.zero_sum_nonzero_count(size, q))
            for p in (2, 3, 4):
                scale = F(q)**(p*b-width)
                exact = scale*sum(F(n**p, (q-1)**((p-1)*mask.bit_count())) for mask, n in counts.items())
                subset = scale*sum(F(n**p, (q-1)**((p-1)*mask.bit_count()))
                    for mask, n in counts.items() if mask == 0 or mask in rank_one_masks)
                bound = target.rank_one_moment_part(sizes, b, p, q=q, width=width)
                self.assertEqual(bound, subset)
                self.assertLessEqual(bound, exact)

    def test_all_pair_branches_vs_GF4_and_raw_mass(self):
        supports = target.concrete_candidates()[2]
        h = target.matrix(supports, 2, 4)
        profile = target.pair_profile(supports, lambda_value=2, q=4)
        self.assertEqual(profile['pair_count'], 120)
        self.assertGreater(profile['two_active_minimum_output_weight'], 0)
        for row in profile['pairs']:
            a, b = row['source']
            actual = Counter()
            for x, y in product(range(4), repeat=2):
                output = apply(h, tuple(x if i == a else y if i == b else 0 for i in range(16)), 4)
                actual[int(bool(x))+int(bool(y)), sum(bool(v) << i for i, v in enumerate(output))] += 1
            expected = Counter()
            for w, mask, count in row['branches']:
                expected[w, mask] += count
            self.assertEqual(actual, +expected)
        minima = []
        for s in target.concrete_candidates():
            profile = target.pair_profile(s)
            minima.append(profile['two_active_minimum_output_weight'])
            self.assertEqual(sum(profile['exceptional_ratio_histogram'].values()), 360)
            self.assertTrue(all(sum(c for _, _, c in row['branches']) == 256**2 for row in profile['pairs']))
        self.assertEqual(minima, [4, 4, 2, 4])

    def test_subset_screen_is_not_full_and_bound_direction_explicit(self):
        supports = target.concrete_candidates()[1]
        row = target.complement_screen(supports, complements=((0, 1, 2, 3, 4), (0, 7, 11, 13, 14)))
        self.assertEqual(row['checked_complements'], 2)
        self.assertTrue(row['complete_requested_scope'])
        self.assertFalse(row['complete_all4368'])
        self.assertTrue(row['bounds_are_lower_not_upper'])
        self.assertEqual(row['regional_provider_runs'], 0)
        self.assertEqual(row['full_support_tables'], 0)
        self.assertFalse(row['old_receipts_transferred'])
        with patch.object(target, 'perf_counter', side_effect=[0, 1, 1]):
            row = target.complement_screen(supports, complements=((0, 1, 2, 3, 4),), max_seconds=0.5)
        self.assertEqual(row['checked_complements'], 0)
        self.assertFalse(row['complete_requested_scope'])

    def test_search_guards_and_no_implicit_production(self):
        with self.assertRaises(ValueError):
            target.screen_candidates(target.concrete_candidates(), max_candidates=2)
        with self.assertRaises(ValueError):
            target.screen_candidates((), max_seconds=31)
        with patch.object(target, 'complement_screen') as screen:
            result = target.screen_candidates(())
        screen.assert_not_called()
        self.assertTrue(result['complete'])
        self.assertTrue(result['no_success_certificate'])


if __name__ == '__main__':
    unittest.main()
