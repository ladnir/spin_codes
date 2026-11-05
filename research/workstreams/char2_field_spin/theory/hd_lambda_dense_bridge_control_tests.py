"""Tiny directed scalar and mocked-controller tests; no real regional bank."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_lambda_dense_bridge_control as target


def toy_terms(*unused, **kwargs):
    return [(even, odd, count, arb(count)*arb(2)**-60)
            for (even, odd), count in target.ARITHMETIC.profile_counts(18).items()]


class DenseBridgeTests(unittest.TestCase):
    def test_directed_density_threshold_and_previous_failure(self):
        for precision in (256, 384):
            with ctx.workprec(precision):
                step = target.ARITHMETIC.uniform.rational(F(8, 5))*target.ARITHMETIC.log2(arb(255))
                result = target.directed_cutoffs(arb(-60)+2*step-arb(1)/4)
                self.assertEqual(result['minimal_certified_total_arity'], 218)
                rows = result['cutoff_ceilings']
                self.assertEqual(len(rows), 73)
                self.assertGreater(F(rows[1]['failure_log2_upper']), -60)
                self.assertLessEqual(F(rows[2]['failure_log2_upper']), -60)
                self.assertEqual([row['total_arity_minimum'] for row in rows], list(range(216, 289)))
                self.assertTrue(all(F(left['failure_log2_upper']) > F(right['failure_log2_upper'])
                                    for left, right in zip(rows, rows[1:])))
                self.assertEqual(target.directed_cutoffs(arb(-60))['minimal_certified_total_arity'], 216)
                none = target.directed_cutoffs(arb(10000))
                self.assertIsNone(none['minimal_certified_total_arity'])
                self.assertTrue(none['no_cutoff_in_requested_range'])

    def test_fresh_union_normalization_and_no_second_source_factor(self):
        with ctx.workprec(256):
            terms = toy_terms()
            result = target.summarize_terms(target.MARKS[0], terms)
            expected = target.ARITHMETIC.high_prefix(18, 216, target.MARKS[0], F(1, 5))
            expected += target.ARITHMETIC.log2(arb(sum(row[2] for row in terms)))-60
            self.assertEqual(result['failure_log2_upper'],
                             str(target.ARITHMETIC.uniform.rational_ceiling(expected)))
            self.assertEqual(result['profile_count'], 99)
            self.assertEqual(result['cutoff_ceilings'][0]['failure_log2_upper'],
                             result['failure_log2_upper'])
            # The same bank yields all cutoffs; the density change is outside
            # the profile sum and does not introduce a row-set multiplier.
            direct247 = target.ARITHMETIC.high_prefix(18, 247, target.MARKS[0], F(1, 5))
            direct247 += target.ARITHMETIC.log2(arb(sum(row[2] for row in terms)))-60
            self.assertEqual(result['cutoff_ceilings'][31]['failure_log2_upper'],
                             str(target.ARITHMETIC.uniform.rational_ceiling(direct247)))

    def test_three_serial_fresh_calls_and_construction_binding(self):
        calls = []
        def bank(mark, alpha, maximum, bits):
            calls.append((mark, alpha, maximum, bits))
            return {'mock': True}
        with patch.object(target.ARITHMETIC.uniform, 'regional_directed', side_effect=bank), \
             patch.object(target.ARITHMETIC, 'source_profile_terms', side_effect=toy_terms):
            result = target.receipt(256)
        self.assertEqual(calls, [(mark, F(1, 5), 18, 60) for mark in target.MARKS])
        self.assertTrue(result['complete'])
        self.assertEqual(result['completed_controls'], 3)
        self.assertEqual(result['attempted_fresh_regional_banks'], 3)
        self.assertEqual(result['total_arity_minimum'], 216)
        self.assertEqual(result['lambda_aes_byte'], 0xbc)
        self.assertEqual(result['construction'], target.binding.SCOPE['construction'])
        self.assertEqual(result['source_count_method'], 'profile_count_upper')
        self.assertEqual(result['rowset_count_factor'], 1)
        self.assertTrue(result['source_support_union_outside_alpha'])
        self.assertFalse(result['saved_receipts_used'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertTrue(result['no_extra_union_over_marks'])
        best = min(F(row['failure_log2_upper']) for row in result['controls'])
        self.assertEqual(F(result['best_completed_cutoff216_log2_upper']), best)

    def test_profile_incompleteness_and_count_tampering_rejected(self):
        with ctx.workprec(256):
            terms = toy_terms()
            for bad in (terms[:-1], terms+[terms[0]],
                        [(terms[0][0], terms[0][1], terms[0][2]+1, terms[0][3])]+terms[1:]):
                with self.assertRaises(ValueError):
                    target.summarize_terms(target.MARKS[0], bad)
            with self.assertRaises(ValueError):
                target.summarize_terms(F(99, 100), terms)

    def test_timeout_preserves_only_completed_controls(self):
        with patch.object(target, 'checked_dependencies', return_value={'mock': 'pin',
                 'hd_lambda_dense_bridge_control.py': 'self'}), \
             patch.object(target.ARITHMETIC.uniform, 'regional_directed',
                          side_effect=[{}, TimeoutError('toy budget')]), \
             patch.object(target.ARITHMETIC, 'source_profile_terms', side_effect=toy_terms):
            result = target.receipt(256)
        self.assertFalse(result['complete'])
        self.assertEqual(result['completed_controls'], 1)
        self.assertEqual(result['attempted_fresh_regional_banks'], 2)
        self.assertEqual(result['stop_reason']['mark'], str(target.MARKS[1]))
        self.assertFalse(result['all_three_marks_tested'])
        self.assertFalse(result['whole_code_certificate'])

    def test_pin_changes_and_bad_resource_parameters_fail(self):
        for kwargs in ({'precision': True}, {'precision': 512},
                       {'max_seconds_per_mark': True}, {'max_seconds_per_mark': 61},
                       {'max_seconds_per_mark': float('nan')}):
            with self.assertRaises(ValueError):
                target.receipt(**kwargs)
        with patch.object(target, 'checked_dependencies', side_effect=[{'pin': 'old'}, {'pin': 'new'}]), \
             patch.object(target.ARITHMETIC.uniform, 'regional_directed', return_value={}), \
             patch.object(target.ARITHMETIC, 'source_profile_terms', side_effect=toy_terms):
            with self.assertRaises(ValueError):
                target.receipt(256)


if __name__ == '__main__':
    unittest.main()
