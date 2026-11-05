"""Cheap completeness, normalization and deadline controls for the one-bank screen."""
from contextlib import ExitStack
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_affine_sparse_all_sources as target


class AllAffineSourcesTests(unittest.TestCase):
    def mocks(self, stack):
        counts = {(a, b): 1 for a in range(9) for b in range(9)}
        counts[0, 0] = comb(8192, 17)-80
        stack.enter_context(patch.object(target, 'matching_joint_counts', return_value=counts))
        raw = stack.enter_context(patch.object(target, 'pair_regions_screen', return_value=object()))
        ordered = stack.enter_context(patch.object(target, 'ordered_region_tables', return_value=None))
        action = stack.enter_context(patch.object(target.affine, 'column_action_float', return_value=np.array([0., -63.5])))
        polynomial = stack.enter_context(patch.object(target, 'sparse_log_count', return_value=-10000.))
        return counts, raw, ordered, action, polynomial

    def test_full_partition_one_bank_and_once_only_class_counts(self):
        with ExitStack() as stack:
            counts, raw, ordered, action, polynomial = self.mocks(stack)
            result = target.screen()
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(ordered.call_count, 81)
        self.assertEqual(action.call_count, 81*18*19//2)
        self.assertEqual(polynomial.call_count, 81*18)
        self.assertTrue(result['complete_source_coverage'])
        self.assertEqual(result['selected_cells'], 1458)
        self.assertEqual(result['completed_source_supports'], comb(8192, 17))
        self.assertEqual(result['mark'], '996248217/1000000000')
        self.assertEqual(result['log2_potential'], '-127/2')
        for source in result['sources']:
            key = tuple(source['source_matching_counts'])
            for n1, row in enumerate(source['entries']):
                n2 = 17-n1
                self.assertEqual(row['source_supports'], counts[key])
                self.assertEqual(row['source_colorings'], comb(17, n1))
                self.assertEqual(row['singleton_row_choices'], 8**n1)
                self.assertEqual(row['pair_row_choices'], 28**n2)
                self.assertEqual(row['sparse_support_floor'], 17*(n1+2*n2))
                self.assertEqual(row['sparse_normalizer_exponent'], 16*(n1+2*n2))
        common = -log2(255)-target.LIMIT*log2(float(target.MARK))-float(target.LOG_POTENTIAL)-10000.
        expected = log2(comb(8192, 17))+17*log2(36)+common
        self.assertAlmostEqual(result['failure_log2_screen'], expected, places=9)
        self.assertTrue(result['regional_matching_count_envelopes_used'])
        for key in ('column_occupancy_resource_maximum_used', 'pair_multiplicity_strata_used',
                    'directed', 'old_certificates_transferred', 'whole_code_certificate'):
            self.assertFalse(result[key])

    def test_deadline_before_raw_bank_has_no_complete_union(self):
        with ExitStack() as stack:
            _, raw, _, _, _ = self.mocks(stack)
            stack.enter_context(patch.object(target, 'perf_counter', side_effect=[0., 121., 122.]))
            result = target.screen()
        self.assertTrue(result['budget_exhausted'])
        self.assertFalse(result['complete_source_coverage'])
        self.assertIsNone(result['failure_log2_screen'])
        self.assertEqual(result['sources'], [])
        raw.assert_not_called()

    def test_partial_source_does_not_count_as_covered(self):
        with ExitStack() as stack:
            _, _, _, action, _ = self.mocks(stack)
            action.side_effect = TimeoutError('bounded toy interruption')
            result = target.screen()
        self.assertTrue(result['budget_exhausted'])
        self.assertEqual(result['source_conditioned_cap_banks_attempted'], 1)
        self.assertEqual(result['source_classes'], 0)
        self.assertEqual(result['selected_cells'], 0)
        self.assertIsNone(result['failure_log2_screen'])

    def test_source_partition_fails_closed_before_expensive_work(self):
        for counts in ({(0, 0): comb(8192, 17)},
                       {(a, b): 1 for a in range(9) for b in range(9)}):
            with patch.object(target, 'matching_joint_counts', return_value=counts), \
                    patch.object(target, 'pair_regions_screen') as raw:
                with self.assertRaises(ArithmeticError):
                    target.screen()
                raw.assert_not_called()

    def test_budget_cannot_expand_authorized_scope(self):
        for budget in (0, -1, 121, True, float('nan')):
            with self.assertRaises(ValueError):
                target.screen(budget_seconds=budget)


if __name__ == '__main__':
    unittest.main()
