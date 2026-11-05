"""Fixed-class assignment and per-mask search-budget tests, no proof rerun."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

import numpy as np

from hd_xor_mixed_mask_refine import (ALL_MASKS, INHERITED, A8_SEARCH,
    combine_fixed_banks, search_one)
from hd_xor_mixed_mask_search import BOUNDS
from hd_xor_mixed_mask_witnesses import checked_classes, PROPOSAL_SHA256


class MixedMaskRefineTests(unittest.TestCase):
    def test_frozen_assignment_matches_pinned_complete_proposal(self):
        path = INHERITED.with_name('hd_xor_mixed_mask_8_0_mu16_refine64.json')
        self.assertEqual(sha256(path.read_bytes()).hexdigest(), PROPOSAL_SHA256)
        proposal = json.loads(path.read_text())['combined']
        assigned = {}
        for label, z, logp, masks in checked_classes():
            self.assertEqual(F(proposal['witnesses'][label]['mark']), z)
            self.assertEqual(F(proposal['witnesses'][label]['log2_potential']), logp)
            for mask in masks:
                assigned[mask] = label
        self.assertEqual(len(assigned), 219)
        self.assertEqual(assigned, {row['active_row_mask']: row['witness_label'] for row in proposal['entries']})
        self.assertLess(proposal['failure_log2_screen'], -100)

    def test_actual_survivors_and_classwise_assignment(self):
        first = json.loads(INHERITED.read_text())
        second = json.loads(A8_SEARCH.read_text())['selected_all_dense_row_sets']
        result = combine_fixed_banks(dict(inherited=first, a8_fit=second))
        self.assertEqual(len(result['entries']), 219)
        self.assertEqual([row['active_row_mask'] for row in result['entries']], list(ALL_MASKS))
        for x, y, chosen in zip(first['entries'], second['entries'], result['entries']):
            self.assertEqual(chosen['failure_log2_screen'], min(x['failure_log2_screen'], y['failure_log2_screen']))
        self.assertAlmostEqual(result['failure_log2_screen'], 38.301069654511835, places=9)
        selected = [max((row for row in result['entries'] if row['active_constituents'] == arity),
                        key=lambda row: row['failure_log2_screen'])['active_row_mask'] for arity in (3, 4)]
        self.assertEqual(selected, [7, 240])
        self.assertTrue(result['witness_choice_is_fixed_by_active_row_mask'])
        self.assertFalse(result['setup_dependent_selection'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_wrong_scope_and_duplicate_class_rejected(self):
        bank = json.loads(INHERITED.read_text())
        for key, value in (('source_matching_counts', [0, 8]), ('mu_range', [15, 16]),
                           ('sparse_support_floor', 543), ('dense_nonempty_column_floor', 16),
                           ('sparse_support_normalizer_exponent', 511), ('directed', True)):
            bad = deepcopy(bank)
            bad[key] = value
            with self.assertRaises(ValueError):
                combine_fixed_banks(dict(bad=bad))
        bad = deepcopy(bank)
        bad['entries'][0] = bad['entries'][1]
        with self.assertRaises(ValueError):
            combine_fixed_banks(dict(bad=bad))

    def test_each_search_has32_evaluation_budget_and_one_fixed_A(self):
        calls = []

        def fake(mark, potential, active_masks):
            z, p = float(F(mark)), float(F(potential))
            calls.append((z, p, tuple(active_masks)))
            value = 100*((z-.9962)/.002)**2+((p+65)/5)**2
            return dict(failure_log2_screen=value, elapsed_seconds=0.,
                        entries=[dict(largest_dense_degree=32, largest_sparse_degree_at_peak=1000)])

        initial = dict(mark='996598217/1000000000', log2_potential='-135/2')
        with patch('hd_xor_mixed_mask_refine.evaluate_at', side_effect=fake):
            result = search_one(7, initial, -1.)
        self.assertLessEqual(len(calls), 32)
        self.assertEqual(len(calls), result['evaluations'])
        self.assertTrue(all(masks == (7,) for _, _, masks in calls))
        for z, p, _ in calls:
            self.assertTrue(all(low <= value <= high for value, (low, high) in zip((z, p), BOUNDS)))
        self.assertEqual(result['active_row_mask'], 7)
        self.assertEqual(result['active_constituents'], 3)
        self.assertFalse(result['directed'])


if __name__ == '__main__':
    unittest.main()
