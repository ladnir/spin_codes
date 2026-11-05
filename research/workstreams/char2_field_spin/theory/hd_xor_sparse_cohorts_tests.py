"""Cheap exact count/scope and shared-bank controls; no proof sweep in tests."""
from fractions import Fraction as F
from math import comb, log2
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_xor_sparse_cohorts as gate
from hd_xor_geometry_strata import category_resources
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import sparse_log_count


class SparseCohortsTests(unittest.TestCase):
    def test_bins_cover_each_multiplicity_once(self):
        for n2 in range(18):
            expected = list(range(1, n2+1)) if n2 else [0]
            self.assertEqual([mu for low, high in gate.cohort_bins(n2) for mu in range(low, high+1)], expected)
            self.assertEqual(sum(count for _, _, count in gate.geometry_bins(n2)), 28**n2)
        for bad in (-1, 18, True, 2.5):
            with self.assertRaises(ValueError):
                gate.cohort_bins(bad)

    def test_cohort_counts_partition_all_one_or_two_row_choices(self):
        actual = sum(comb(17, n1)*8**n1*sum(count for _, _, count in gate.geometry_bins(17-n1))
                     for n1 in range(18))
        self.assertEqual(actual, 36**17)

    def test_class_factors_and_floor(self):
        caps = {key: -sum(key) for key in categories(17)}
        with patch.object(gate, 'sparse_log_count', return_value=123.) as polynomial:
            row = gate.class_row(1, (9, 14), 19, 23, caps, 'test', F(249, 250), F(-63))
        polynomial.assert_called_once_with(caps, 1, 16)
        expected = log2(23)+log2(17)+3+log2(19)-log2(255)-gate.LIMIT*log2(249/250)+63+123
        self.assertAlmostEqual(row['failure_log2_screen'], expected)
        self.assertEqual((row['sparse_support_floor'], row['sparse_support_normalizer_exponent']), (561, 528))
        self.assertEqual((row['cohort_colorings'], row['one_row_choices']), (17, 8))

    def test_pure_one_ignores_all_double_caps(self):
        first = {key: -.1*key[0]-.2*key[1] for key in categories(17)}
        second = {key: value+(1000 if key[1] else 0) for key, value in first.items()}
        self.assertEqual(sparse_log_count(first, 17, 0), sparse_log_count(second, 17, 0))

    def test_canonical_upper_caps_equal_truncated_caps_on_used_categories(self):
        # A possible future replay saving, not used by the frozen control:
        # only d<=n2 is selected by a cohort's support polynomial.
        for n2 in range(18):
            for low, high in gate.cohort_bins(n2):
                canonical = next(upper for lower, upper in gate.CANONICAL_BINS if lower <= max(1, high) <= upper)
                for double in range(n2+1):
                    for single in range(18-double):
                        for cancelled in range(double+1):
                            self.assertEqual(category_resources(single, double, cancelled, max(1, high)),
                                             category_resources(single, double, cancelled, canonical))

    def test_two_banks_reuse_each_cap_and_choose_without_bank_union(self):
        evaluators = []
        class FakeCaps:
            def __init__(self, tables, potential, groups):
                self.calls = []
                self.bank = len(evaluators)
                self._envelope = SimpleNamespace(cache_info=lambda: SimpleNamespace(currsize=0))
                evaluators.append(self)
            def __call__(self, mu):
                self.calls.append(mu)
                return {'value': self.bank*100+mu}
        with patch.object(gate, 'matching_joint_counts', return_value={(8, 0): 23}) as counts, \
                patch.object(gate, 'pair_regions_screen', return_value=None) as raw, \
                patch.object(gate, 'ordered_region_tables', return_value=None), \
                patch.object(gate, 'MultiplicityCapsFloat', FakeCaps), \
                patch.object(gate, 'sparse_log_count', side_effect=lambda caps, n1, n2: caps['value']+n1):
            result = gate.run()
        self.assertEqual(counts.call_count, 1)
        self.assertEqual(raw.call_count, 2)
        self.assertEqual(len(evaluators), 2)
        self.assertTrue(all(e.calls == list(range(1, 18)) for e in evaluators))
        self.assertEqual(result['cohort_count'], 18)
        self.assertEqual([row['cohort_counts'] for row in result['cohorts']], [[n, 17-n, 0] for n in range(18)])
        self.assertTrue(result['complete_sparse_cohorts_for_selected_source'])
        self.assertFalse(result['complete_source_coverage'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['fitted_new_witnesses'], 0)
        self.assertEqual(len(result['banks'][0]['entries']), 33)
        for cohort in result['cohorts']:
            for chosen in cohort['entries']:
                alternatives = [row['failure_log2_screen'] for bank in result['banks'] for row in bank['entries']
                    if row['cohort_counts'] == chosen['cohort_counts'] and row['mu_range'] == chosen['mu_range']]
                self.assertEqual(chosen['failure_log2_screen'], min(alternatives))
        self.assertEqual(result['failure_log2_screen'], float(np.logaddexp2.reduce([
            row['failure_log2_screen'] for row in result['cohorts']])))

    def test_all_sources_checks_joint_partition_before_building_tables(self):
        with patch.object(gate, 'matching_joint_counts', return_value={(8, 0): 23}), \
                patch.object(gate, 'pair_regions_screen') as raw:
            with self.assertRaises(ArithmeticError):
                gate.run(all_sources=True)
            raw.assert_not_called()
        for invalid in (1, None, 'all'):
            with self.assertRaises(ValueError):
                gate.run(all_sources=invalid)

    def test_all_sources_are_enumerated_and_sum_cohort_masses(self):
        counts = {(a, b): 1 for a in range(9) for b in range(9)}
        counts[8, 0] = comb(gate.REGION_LENGTH, 17)-80
        def fake(source, count, partitions, raw_banks, **kwargs):
            cohorts = []
            for n1 in range(18):
                value = -1000+sum(source)+n1
                rows = [dict(source_matching_counts=list(source), cohort_counts=[n1, 17-n1, 0],
                    mu_range=[low, high], failure_log2_screen=value-10)
                    for low, high, _ in partitions[17-n1]]
                cohorts.append(dict(cohort_counts=[n1, 17-n1, 0], entries=rows, failure_log2_screen=value))
            return dict(source_matching_counts=list(source), cohorts=cohorts,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts])),
                failed_cohorts_at_40_bits=[])
        with patch.object(gate, 'matching_joint_counts', return_value=counts), \
                patch.object(gate, 'pair_regions_screen', return_value=None) as raw, \
                patch.object(gate, '_source_screen', side_effect=fake) as source_screen:
            result = gate.run(all_sources=True)
        self.assertEqual(raw.call_count, 2)
        self.assertEqual(source_screen.call_count, 81)
        self.assertEqual({call.args[0] for call in source_screen.call_args_list}, set(counts))
        self.assertEqual(result['selected_cells'], 81*33)
        self.assertEqual(result['source_supports'], comb(gate.REGION_LENGTH, 17))
        self.assertTrue(result['complete_source_coverage'])
        self.assertTrue(result['complete_sparse_cohort_coverage'])
        self.assertFalse(result['whole_code_certificate'])
        for n1, cohort in enumerate(result['cohorts']):
            expected = float(np.logaddexp2.reduce([-1000+a+b+n1 for a, b in counts]))
            self.assertEqual(cohort['failure_log2_screen'], expected)
        self.assertEqual(result['failed_cells_at_40_bits'], [])


if __name__ == '__main__':
    unittest.main()
