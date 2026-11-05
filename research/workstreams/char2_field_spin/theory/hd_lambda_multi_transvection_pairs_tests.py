"""Tiny equivalence tests only; the full9730-configuration screen is not run."""
from unittest.mock import patch
import unittest

import hd_lambda_multi_transvection_pairs as target


class PairScreenTests(unittest.TestCase):
    def test_fast_columns_are_full_matrix_transpose(self):
        for supports in target.core.concrete_candidates():
            rows = target.core.binary_rows(supports)
            self.assertEqual(target.binary_columns(supports),
                tuple(sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(16)))

    def test_all_example_masks_match_GF4_and_GF256_exceptional_ratios(self):
        for supports in target.core.concrete_candidates():
            columns = target.binary_columns(supports)
            for q, lam in ((4, 2), (256, 188)):
                slow = target.core.pair_profile(supports, lambda_value=lam, q=q)
                for row in slow['pairs']:
                    a, b = row['source']
                    self.assertEqual(target.exceptional_masks(columns[a], columns[b]),
                                     tuple(mask for _, mask, _ in row['branches'][3:6]))
                self.assertEqual(target.exact_pair_minimum(supports)['minimum_output_weight'],
                                 slow['two_active_minimum_output_weight'])

    def test_examples_have_expected_costs_and_minima(self):
        receipt = target.example_receipt()
        self.assertEqual([row['circuit']['extra_xors_over_Hlambda'] for row in receipt['rows']], [8, 9, 14, 16])
        self.assertEqual([row['pair_profile']['minimum_output_weight'] for row in receipt['rows']], [4, 4, 2, 4])
        self.assertFalse(receipt['complete_family_search'])
        self.assertEqual(receipt['regional_provider_runs'], 0)

    def test_budget_exit_does_not_claim_coverage(self):
        with patch.object(target, 'perf_counter', side_effect=[0, 1, 1]):
            result = target.complete_pair_screen(max_seconds=0.5)
        self.assertFalse(result['complete'])
        self.assertEqual(result['considered_pairs'], 0)
        self.assertEqual(result['survivor_count'], 0)
        self.assertIsNotNone(result['stop_reason'])

    def test_guards(self):
        for budget in (0, -1, True, 31, float('inf')):
            with self.assertRaises(ValueError):
                target.complete_pair_screen(max_seconds=budget)
        for a, b in ((-1, 0), (0, 65536), (True, 1)):
            with self.assertRaises(ValueError):
                target.exceptional_masks(a, b)


if __name__ == '__main__':
    unittest.main()
