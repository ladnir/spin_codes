"""Tiny max-plus/exhaustive checks only; no saved production diagnostic run."""
from itertools import product
from math import log2
import unittest

import numpy as np

import hd_lambda_double_method_path as target
from hd_pair_identity_kernel import mask_table_from_entries


class MethodPathTests(unittest.TestCase):
    def test_maxplus_matches_exhaustive_ordered_term_and_row_counts(self):
        columns = []
        for j in range(3):
            columns.append([dict(kind=str(mask), original_mask=mask,
                norm_log2=np.log2(np.array(((2+j+mask, 3+mask), (5+j, 1+j+2*mask)), dtype=float)))
                for mask in (0, 2, 3)])
        best = -np.inf
        for choices in product(range(3), repeat=3):
            masks = [columns[j][c]['original_mask'] for j, c in enumerate(choices)]
            if any(sum(bool(mask & bit) for mask in masks) < 2 for bit in (1, 2)): continue
            for ends in product((0, 1), repeat=3):
                states = (0, *ends)
                value = sum(2*columns[j][c]['norm_log2'][states[j], states[j+1]] for j, c in enumerate(choices))
                best = max(best, value)
        actual = target.best_repeated_branch_term(columns, groups=2, floor=2, row_classes=(1, 2))
        self.assertAlmostEqual(actual['term_log2'], best)
        self.assertGreaterEqual(min(actual['original_row_activity_counts'].values()), 2)
        self.assertTrue(actual['identical_branch_in_all_replicas'])

    def test_same_branch_term_is_lower_than_full_expansion(self):
        # Every replica may instead choose its own branch; sum all such terms.
        columns = [[dict(kind=str(mask), original_mask=mask,
            norm_log2=np.log2(np.array(((1+mask+j, 2+mask), (3+j, 4+mask)), dtype=float)))
            for mask in (0, 2, 3)] for j in range(2)]
        total = 0.
        for selections in product(range(3), repeat=4):
            activity = [0, 0]
            for index, selected in enumerate(selections):
                mask = columns[index//2][selected]['original_mask']
                activity[0] += bool(mask & 1)
                activity[1] += bool(mask & 2)
            if min(activity) < 2: continue
            for middle, end in product((0, 1), repeat=2):
                states = (0, middle, end)
                value = 1.
                for j in range(2):
                    for k in range(2):
                        value *= 2**columns[j][selections[2*j+k]]['norm_log2'][states[j], states[j+1]]
                total += value
        witness = target.best_repeated_branch_term(columns, groups=2, floor=1, row_classes=(1, 2))
        self.assertLessEqual(witness['term_log2'], log2(total))

    def test_original_support_mapping_not_just_column_index(self):
        for j in range(16):
            before, after = target.previous.base.column_permutations(target.previous.base.fiber.SCHEDULE, j)
            for original in (0, 63, 60):
                physical = sum(1 << after[before[i]] for i in range(6) if original & (1 << i))
                self.assertEqual(target.original_mask(physical, j), original)

    def test_mask_moment_retains_noncommuting_order_and_exact_zeros(self):
        zero = np.array(((0., log2(2)), (-np.inf, log2(3))))
        nonzero = np.array(((log2(5), -np.inf), (log2(7), 0.)))
        table = mask_table_from_entries(zero, nonzero, rows=4)
        for mask in range(16):
            self.assertTrue(np.allclose(target.mask_moment(mask, zero, nonzero, rows=4), table[mask]))
        self.assertFalse(np.allclose(table[1], table[8]))

    def test_reject_partial_counter_classes_and_infeasible_paths(self):
        partial = [[dict(kind='bad', original_mask=1, norm_log2=np.zeros((2, 2)))]]
        with self.assertRaises(ValueError):
            target.best_repeated_branch_term(partial, floor=1, row_classes=(3, 12))
        empty = [[dict(kind='zero', original_mask=0, norm_log2=np.zeros((2, 2)))]]
        witness = target.best_repeated_branch_term(empty, floor=1, row_classes=(3, 12))
        self.assertFalse(witness['found'])


if __name__ == '__main__':
    unittest.main()
