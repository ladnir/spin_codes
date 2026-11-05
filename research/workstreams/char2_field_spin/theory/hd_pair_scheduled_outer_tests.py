"""Tiny exact controls for deterministic scheduled Hbin branch preparation."""
from collections import Counter
from itertools import combinations, product
import unittest

import numpy as np

import hd_pair_scheduled_outer as scheduled
from hd_pair_state_control import branch_data as static_branch_data


def reference_multiply(a, b):
    value = 0
    for bit in range(4):
        if (b >> bit) & 1:
            value ^= a << bit
    for bit in range(6, 3, -1):
        if (value >> bit) & 1:
            value ^= 0x13 << (bit-4)
    return value


def reference_permutations(name, column):
    inverse = [0]+[next(j for j in range(1, 16) if reference_multiply(i, j) == 1)
                  for i in range(1, 16)]
    before, after = list(range(16)), list(range(16))
    if name == 'input-xor':
        before = [i ^ column for i in range(16)]
    if name in ('input-inverse', 'input-inverse-output-xor'):
        before = [inverse[i ^ column] for i in range(16)]
    if name in ('output-xor', 'input-inverse-output-xor'):
        after = [i ^ column for i in range(16)]
    return before, after


def reference_counts(source, name, column, q=4):
    """Exhaust field inputs through an independently formed dense matrix."""
    before, after = reference_permutations(name, column)
    generator = [[int(row == before[i]) ^ 1 ^ ((row & before[i]).bit_count() & 1)
                  for i in source] for row in range(16)]
    counts = Counter()
    for values in product(range(q), repeat=len(source)):
        mask = 0
        for row, coefficients in enumerate(generator):
            value = 0
            for coefficient, x in zip(coefficients, values):
                if coefficient:
                    value ^= x
            if value:
                mask |= 1 << after[row]
        counts[sum(value != 0 for value in values), mask] += 1
    return counts


class ScheduledOuterTests(unittest.TestCase):
    def test_gf16_inverse_and_all_multiplications(self):
        for a, b in product(range(16), repeat=2):
            self.assertEqual(scheduled.gf16_multiply(a, b), reference_multiply(a, b))
        inverse = scheduled.inverse_table()
        self.assertEqual(inverse[0], 0)
        self.assertEqual(inverse[2], 9)
        self.assertEqual(set(inverse), set(range(16)))
        for value in range(1, 16):
            self.assertEqual(reference_multiply(value, inverse[value]), 1)

    def test_all_prescribed_permutation_directions(self):
        for name, column in product(scheduled.SCHEDULES, range(16)):
            p, q = scheduled.column_permutations(name, column)
            expected_p, expected_q = reference_permutations(name, column)
            self.assertEqual(p, tuple(expected_p))
            self.assertEqual(q, tuple(expected_q))
            self.assertEqual(set(p), set(range(16)))
            self.assertEqual(set(q), set(range(16)))
        # Inverse(x xor j) is not inverse(x) xor j; bind the intended order.
        p, _ = scheduled.column_permutations('input-inverse', 7)
        self.assertNotEqual(p, tuple(x ^ 7 for x in scheduled.inverse_table()))

    def test_support_permutation(self):
        for name, column in product(scheduled.SCHEDULES, range(16)):
            for permutation in scheduled.column_permutations(name, column):
                inverse = tuple(permutation.index(i) for i in range(16))
                for mask in (0, 65535, 0x1234, 0x891a, *[1 << i for i in range(16)]):
                    result = scheduled.permute_mask(mask, permutation)
                    expected = sum(1 << permutation[i] for i in range(16) if mask >> i & 1)
                    self.assertEqual(result, expected)
                    self.assertEqual(result.bit_count(), mask.bit_count())
                    self.assertEqual(scheduled.permute_mask(result, inverse), mask)

    def test_gf4_all_pairs_selected_columns_all_schedules(self):
        # All120 original pairs; three columns include asymmetric inversions
        # and the last column. This is local exact algebra, not a proof gate.
        for source, name, column in product(combinations(range(16), 2), scheduled.SCHEDULES, (0, 7, 15)):
            got = Counter()
            for weight, mask, count in scheduled.column_branches(source, name, column, q=4):
                got[weight, mask] += count
            self.assertEqual(got, reference_counts(source, name, column))
            self.assertEqual(sum(got.values()), 16)

    def test_gf4_all_singletons_columns_and_schedules(self):
        for source, name, column in product(combinations(range(16), 1), scheduled.SCHEDULES, range(16)):
            got = Counter({(weight, mask): count for weight, mask, count in
                           scheduled.column_branches(source, name, column, q=4)})
            self.assertEqual(got, reference_counts(source, name, column))

    def test_static_exactly_matches_existing_branch_arrays(self):
        for arity in (1, 2):
            got = scheduled.branch_data(arity)
            old = static_branch_data(arity, mixer='Hbin')
            self.assertEqual(got['sources'], old['sources'])
            self.assertEqual(got['branch_count'], old['branch_count'])
            for key in ('log_counts', 'input_weights', 'activities'):
                np.testing.assert_array_equal(got[key], old[key])
            for block in range(8):
                np.testing.assert_array_equal(got['first'][block], old['first'])
                np.testing.assert_array_equal(got['second'][block], old['second'])
            self.assertTrue(np.all(got['counts'].sum(axis=1) == 256**(2*arity)))

    def test_nonstatic_fixed_source_block_shapes_and_weights(self):
        data = scheduled.branch_data(2, 'input-inverse-output-xor', q=4)
        self.assertEqual(data['first'].shape, (8, 120*25))
        self.assertEqual(data['second'].shape, (8, 120*25))
        self.assertEqual(data['input_weights'].shape, (120, 25))
        self.assertEqual(data['sources'], tuple(combinations(range(16), 2)))
        self.assertEqual(data['deltas'], tuple(range(8, 16)))
        self.assertTrue(np.all(data['counts'].sum(axis=1) == 256))
        for source_index in (0, 37, 119):
            source = data['sources'][source_index]
            for block in range(8):
                left = scheduled.column_branches(source, 'input-inverse-output-xor', 2*block, q=4)
                right = scheduled.column_branches(source, 'input-inverse-output-xor', 2*block+1, q=4)
                for index, (a, b) in enumerate(product(left, right)):
                    flat = source_index*25+index
                    self.assertEqual(int(data['first'][block, flat]), a[1])
                    self.assertEqual(int(data['second'][block, flat]), b[1])
                    self.assertEqual(int(data['input_weights'][source_index, index]), a[0]+b[0])
                    self.assertEqual(int(data['counts'][source_index, index]), a[2]*b[2])
        self.assertFalse(np.array_equal(data['first'][0], data['second'][0]))
        self.assertFalse(np.array_equal(data['first'][0], data['first'][1]))

    def test_category_identity_follows_original_not_sorted_mapped_rows(self):
        source = (0, 1)
        p, q = scheduled.column_permutations('input-inverse-output-xor', 3)
        self.assertGreater(p[source[0]], p[source[1]])
        rows = scheduled.column_branches(source, 'input-inverse-output-xor', 3)
        only_first = reference_counts(source[:1], 'input-inverse-output-xor', 3)
        self.assertIn((1, rows[1][1]), only_first)
        self.assertEqual(tuple(weight for weight, _, _ in rows), (0, 1, 1, 2, 2))

    def test_validation_and_scope(self):
        for arity in (0, 3, True):
            with self.assertRaises(ValueError):
                scheduled.branch_data(arity)
        for schedule, column in (('random', 0), ('static', 16), ('static', -1)):
            with self.assertRaises(ValueError):
                scheduled.column_permutations(schedule, column)
        for source in ((1, 0), (1, 1), (-1,), (16,)):
            with self.assertRaises(ValueError):
                scheduled.column_branches(source, 'static', 0)
        with self.assertRaises(ValueError):
            scheduled.permute_mask(1, (0,)*16)
        info = scheduled.metadata('input-inverse-output-xor')
        self.assertTrue(info['original_A_fixed_across_all_columns'])
        self.assertFalse(info['rowset_resampling'])
        self.assertFalse(info['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
