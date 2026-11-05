from collections import Counter
from itertools import product
from math import comb
import unittest

import hd_lambda_arity5_relations as new
import hd_lambda_scheduled_relations as old
import hd_quaternary_outer_relations as base
from hd_pair_scheduled_outer import column_permutations


def gf4_multiply(a, b):
    # Independent polynomial arithmetic for the tiny exhaustive oracle.
    value = (a if b & 1 else 0) ^ ((a << 1) if b & 2 else 0)
    return value ^ (7 if value & 4 else 0)


class ArityFiveRelationsTests(unittest.TestCase):
    def test_category_count_and_exact_GF256_partition(self):
        data = new.category_data(5)
        self.assertEqual(len(data['kernel_rrefs']), 12278)
        self.assertEqual(Counter(map(len, data['kernel_rrefs'])),
                         {0: 1, 1: 341, 2: 5797, 3: 5797, 4: 341, 5: 1})
        self.assertEqual(sum(data['category_counts']), 256**5)
        self.assertEqual(data['category_counts'][0], 0)  # No GF4^5 injection into GF256.
        self.assertEqual(sum(c == 0 for c in data['category_counts']), 1)
        self.assertEqual(len(set(data['kernel_bitsets'])), 12278)
        by_mask = Counter()
        for rows, bits, mask, weight, count in zip(data['kernel_rrefs'], data['kernel_bitsets'],
                data['category_input_masks'], data['category_input_weights'], data['category_counts']):
            self.assertEqual(bits.bit_count(), 4**len(rows))
            self.assertEqual(weight, mask.bit_count())
            self.assertIs(type(count), int)
            by_mask[mask] += count
        self.assertEqual(dict(by_mask), {s: 255**s.bit_count() for s in range(32)})
        self.assertEqual(data['category_counts'][-1], 1)
        self.assertEqual(data['category_input_masks'][-1], 0)

    def test_GF4_exhaustive_one_fixed_A_and_column(self):
        source = (14, 0, 9, 3, 5)  # Deliberately not sorted.
        schedule, column = 'input-inverse-output-xor', 11
        before, after = column_permutations(schedule, column)
        counts = Counter()
        for x in product(range(4), repeat=5):
            mask = 0
            for row in range(16):
                value = 0
                for i, amplitude in zip(source, x):
                    p = before[i]
                    h = 2 ^ int(row == p) ^ ((row & p).bit_count() & 1)
                    value ^= gf4_multiply(h, amplitude)
                mask |= bool(value) << after[row]
            counts[sum(v != 0 for v in x), mask] += 1
        self.assertEqual(new.joint_counts(source, schedule, column, q=4), counts)
        self.assertEqual(sum(counts.values()), 1024)
        data = new.column_data(source, schedule, column, q=4)
        self.assertEqual(data['original_source'], source)
        self.assertEqual(len(data['branches']), 12278)
        self.assertEqual(sum(c != 0 for c in data['category_counts']), 342)
        self.assertFalse(data['metadata']['same_relation_kernel_across_columns_required'])
        self.assertFalse(data['metadata']['whole_code_certificate'])

    def test_lowarity_agrees_with_frozen_provider(self):
        for arity in range(5):
            for q in (4, 16, 256):
                self.assertEqual(new.category_data(arity, q), base.category_data(arity, q))
        for arity in range(1, 5):
            source = (15, 2, 9, 0)[:arity]
            for schedule in new.SCHEDULES:
                for column in (0, 7, 15):
                    for lam in (0xbc, 0xbd):
                        self.assertEqual(new.column_branches(source, schedule, column, lambda_value=lam),
                                         old.column_branches(source, schedule, column, lambda_value=lam))

    def test_original_input_category_identity_across_schedule(self):
        source = (2, 6, 11, 1, 15)
        left = new.column_data(source, 'static', 0)
        right = new.column_data(source, 'input-inverse-output-xor', 15)
        for name in ('category_counts', 'category_input_masks', 'category_input_weights', 'kernel_rrefs'):
            self.assertEqual(left[name], right[name])
        self.assertNotEqual(left['category_output_masks'], right['category_output_masks'])
        for q in (4, 16, 256):
            data = new.category_data(5, q)
            by_weight = Counter()
            for w, c in zip(data['category_input_weights'], data['category_counts']):
                by_weight[w] += c
            self.assertEqual(dict(by_weight), {w: comb(5, w)*(q-1)**w for w in range(6)})

    def test_physical_output_permutation(self):
        source = (0, 2, 4, 7, 10)
        plain = new.column_masks(source, 'input-inverse', 9)
        moved = new.column_masks(source, 'input-inverse-output-xor', 9)
        self.assertEqual(moved, tuple(sum(((m >> i) & 1) << (i ^ 9) for i in range(16)) for m in plain))

    def test_bounds_do_not_mutate_frozen_provider(self):
        self.assertEqual(base.MAX_ARITY, 4)
        with self.assertRaises(ValueError):
            base.category_data(5)
        for arity in (-1, 6, True):
            with self.assertRaises(ValueError):
                new.category_data(arity)
        for source in ((), (1, 1), tuple(range(6)), (16,), (True,)):
            with self.assertRaises(ValueError):
                new.column_masks(source, 'static', 0)
        with self.assertRaises(ValueError):
            new.column_masks(tuple(range(5)), 'static', 0, lambda_value=2)
        with self.assertRaises(ValueError):
            new.category_data(5, 8)
        with self.assertRaises(ValueError):
            new.masks_from_forms((1024,), 5)


if __name__ == '__main__':
    unittest.main()
