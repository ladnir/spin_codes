"""Lightweight scope, count-transfer and bounded-search orchestration tests."""
from copy import deepcopy
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from hd_xor_mixed_all_geometry_refine import (ALL_MASKS, PENALTY,
    load_banks, all_geometry_summary, run)
from hd_xor_mixed_mask_refine import combine_fixed_banks


class AllGeometryRefineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.banks = load_banks()

    def test_count_transfer_preserves_fixed_assignment_and_all219_masks(self):
        before = deepcopy(self.banks)
        common = combine_fixed_banks(self.banks)
        result = all_geometry_summary(self.banks)
        self.assertEqual(self.banks, before)
        self.assertEqual([row['active_row_mask'] for row in result['entries']], list(ALL_MASKS))
        self.assertEqual(result['sparse_pair_assignments'], 28**16)
        self.assertEqual(result['source_matching_counts'], [8, 0])
        self.assertEqual(result['cohort_counts'], [0, 16, 1])
        for old, new in zip(common['entries'], result['entries']):
            self.assertEqual(old['witness_label'], new['witness_label'])
            self.assertEqual(new['failure_log2_screen'], old['failure_log2_screen']+15*log2(28))
        self.assertAlmostEqual(result['failure_log2_screen'], -37.75362667993719, places=10)
        self.assertFalse(result['directed'])
        self.assertTrue(result['proposal_requires_fresh_directed_replay'])

    def test_wrong_count_scope_and_arity_are_rejected(self):
        for key, value in (('cohort_counts', [1, 15, 1]), ('two_row_geometries', 28**16),
                           ('field_order', 255), ('forbidden_weight', 734004), ('total_columns', 31)):
            banks = deepcopy(self.banks)
            banks['inherited'][key] = value
            with self.assertRaises(ValueError):
                all_geometry_summary(banks)
        banks = deepcopy(self.banks)
        for bank in banks.values():
            bank['entries'][0]['active_constituents'] = 8
        with self.assertRaises(ValueError):
            all_geometry_summary(banks)

    def test_success_stops_after_one_search_and_one_complete_bank(self):
        baseline = all_geometry_summary(self.banks)
        initial = baseline['witnesses']['a3_mask7_fit']
        trial = dict(mark=initial['mark'], log2_potential=initial['log2_potential'],
                     failure_log2_screen=-200.)
        search = dict(evaluations=32, best=trial, proposals=[trial]*32)
        bank = deepcopy(self.banks['inherited'])
        bank['mark'], bank['log2_potential'] = trial['mark'], trial['log2_potential']
        for row in bank['entries']:
            row['failure_log2_screen'] = -200.
        bank['failure_log2_screen'] = float(np.logaddexp2.reduce([-200.]*219))
        with patch('hd_xor_mixed_all_geometry_refine.load_banks', return_value=deepcopy(self.banks)), \
                patch('hd_xor_mixed_all_geometry_refine.search_one', return_value=search) as fit, \
                patch('hd_xor_mixed_all_geometry_refine.evaluate_at', return_value=bank) as evaluate:
            result = run()
        self.assertEqual(fit.call_count, 1)
        self.assertEqual(fit.call_args.args[:2], (195, initial))
        self.assertEqual(evaluate.call_count, 1)
        self.assertEqual(evaluate.call_args.args[2], ALL_MASKS)
        self.assertEqual(result['search_evaluations'], 32)
        self.assertEqual(result['full_bank_evaluations'], 1)
        self.assertLess(result['combined']['failure_log2_screen'], -60.)
        self.assertFalse(result['combined']['setup_dependent_selection'])
        self.assertEqual(result['combined']['exact_mass_multiplier'], 28**15)

    def test_saved_trials_and_complete_assignment_recombine_without_operators(self):
        path = Path(__file__).with_name('hd_xor_mixed_all_geometry_8_0_refine64.json')
        data = path.read_bytes()
        self.assertEqual(sha256(data).hexdigest(),
                         '44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6')
        proposal = json.loads(data)
        self.assertEqual(proposal['search_evaluations'], 64)
        self.assertEqual(proposal['full_bank_evaluations'], 2)
        self.assertEqual([row['active_row_mask'] for row in proposal['searches']], [195, 224])
        for search in proposal['searches']:
            self.assertEqual(search['evaluations'], 32)
            self.assertEqual([row['evaluation'] for row in search['proposals']], list(range(1, 33)))
            self.assertEqual(search['boundary_coordinates'], [])
            self.assertTrue(all(row['active_row_mask'] == search['active_row_mask']
                                for row in search['proposals']))
        banks = dict(self.banks, **proposal['added_banks'])
        self.assertEqual(all_geometry_summary(banks), proposal['combined'])
        self.assertLess(proposal['combined']['failure_log2_screen'], -60.)
        self.assertEqual(sum(row['witness_label'].startswith('all_geometry_')
                             for row in proposal['combined']['entries']), 41)


if __name__ == '__main__':
    unittest.main()
