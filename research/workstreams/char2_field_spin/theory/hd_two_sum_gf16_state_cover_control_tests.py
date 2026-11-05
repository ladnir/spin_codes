"""Mocked six-control scope/budget tests; no coverage bank evaluation."""
from hashlib import sha256
import json
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_state_cover_control as target


def fixture():
    source = dict(construction=target.low.fractional.CONSTRUCTION, field_order=256,
        forbidden_weight=target.low.fractional.LIMIT, source_size=8192,
        directed=False, complete_requested_anchors=True, banks=[], witnesses=[])
    for r, _, mark in target.low.ANCHORS:
        counts = [{'full_max': f, 'count': f+1} for f in range(r//2+1)]
        counts[-1]['count'] += comb(8192, r)-sum(row['count'] for row in counts)
        templates = [dict(name=name, family='nonempty-column-floor' if name == 'old-all-zero'
            else 'pre-H-total-support-floor', markers=['0']*16) for name in target.TEMPLATES]
        bank = dict(active_groups=r, mark=str(mark), alpha=str(target.low.ALPHA), maximum=32,
            log2_slopes=[0.]*4, source_strata=counts, type_templates=templates,
            regions=[dict(occupancy=m, full_pairs=f, log2_matrix=[0., None, None, -1.])
                     for m in range(33) for f in range(m//2+1)])
        source['banks'].append(bank)
        source['witnesses'].extend(dict(active_groups=r, source_full_pair_max=1, template=metadata,
            failure_log2_screen=100+log2(2)-log2(comb(8192, r))) for metadata in templates)
    return source


class StateCoverControlTests(unittest.TestCase):
    def test_count_is_replaced_outside_alpha_and_scope_is_only_f1(self):
        before = dict(failure_log2_screen=100., alpha='1/5', all_source_supports=True)
        row = target.restrict_count(before, 18, 7)
        self.assertAlmostEqual(row['failure_log2_screen'], 100.+log2(7)-log2(comb(8192, 18)))
        self.assertEqual(before['failure_log2_screen'], 100.)
        self.assertFalse(row['all_source_supports'])
        self.assertFalse(row['source_partition_complete'])
        self.assertEqual(row['source_full_pair_max'], 1)
        self.assertTrue(row['source_count_outside_alpha'])

    def test_restricted_curves_use_only_saved_f_zero_and_one(self):
        bank = dict(active_groups=4, regions=[dict(occupancy=m, full_pairs=f,
            log2_matrix=[m+10*f]*4) for m in range(5) for f in range(m//2+1)])
        actual = target.restricted_curves(bank)
        np.testing.assert_array_equal(actual, [[0, 1, 12, 13, 14]]*4)

    def test_complete_mock_run_exactly_six_controls_and_preserves_fixed_scopes(self):
        source = fixture()
        templates = lambda bank: ((row, {}) for row in bank['type_templates'])
        with patch.object(target, 'load_source', return_value=source), \
                patch.object(target, 'fixed_templates', side_effect=templates), \
                patch.object(target.cover, 'prepared_envelope', return_value={}) as prepare, \
                patch.object(target.cover.emptycut, 'evaluate_prepared', return_value={'failure_log2_screen': 100.}), \
                patch.object(target.cover, 'evaluate_prepared', return_value={'failure_log2_screen': 90.}) as evaluate, \
                patch.object(target, 'perf_counter', return_value=0.):
            result = target.run(1)
        self.assertEqual(prepare.call_count, 6)
        self.assertEqual(evaluate.call_count, 6)
        self.assertTrue(result['complete_requested_controls'])
        self.assertEqual(result['regional_banks_rebuilt'], 0)
        self.assertFalse(result['all_source_supports'])
        self.assertFalse(result['continuous_weight_band_covered'])
        self.assertEqual(len(result['selected']), 3)
        for row in result['witnesses']:
            self.assertEqual(row['source_full_pair_max'], 1)
            self.assertEqual(row['source_supports'], 2)
            self.assertAlmostEqual(row['improvement_bits'], 10.)
            self.assertFalse(row['all_source_supports'])
        for index, (r, _, mark) in enumerate(target.low.ANCHORS):
            self.assertEqual(evaluate.call_args_list[2*index].args[1:], (r, mark, target.low.ALPHA))

    def test_timeout_retains_completed_rows_and_never_claims_completion(self):
        source = fixture()
        with patch.object(target, 'load_source', return_value=source), \
                patch.object(target, 'fixed_templates', side_effect=lambda bank: ((row, {}) for row in bank['type_templates'])), \
                patch.object(target.cover, 'prepared_envelope', return_value={}), \
                patch.object(target.cover.emptycut, 'evaluate_prepared', return_value={'failure_log2_screen': 100.}), \
                patch.object(target.cover, 'evaluate_prepared', side_effect=[
                    {'failure_log2_screen': 90.}, TimeoutError('fixture')]), \
                patch.object(target, 'perf_counter', return_value=0.):
            result = target.run(1)
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertEqual(len(result['witnesses']), 1)
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(result['in_progress'], {'active_groups': 18, 'template': 'joint-empty-floor1'})
        self.assertNotIn('selected_three_strata_log2_screen', result)

    def test_saved_baseline_disagreement_fails_before_coverage_evaluation(self):
        source = fixture()
        source['witnesses'][0]['failure_log2_screen'] += 1.
        with patch.object(target, 'load_source', return_value=source), \
                patch.object(target, 'fixed_templates', side_effect=lambda bank: ((row, {}) for row in bank['type_templates'])), \
                patch.object(target.cover, 'prepared_envelope', return_value={}), \
                patch.object(target.cover.emptycut, 'evaluate_prepared', return_value={'failure_log2_screen': 100.}), \
                patch.object(target.cover, 'evaluate_prepared') as evaluate:
            with self.assertRaises(ArithmeticError):
                target.run(1)
        evaluate.assert_not_called()

    def test_literal_markers_reused_without_new_template_fit(self):
        bank = fixture()['banks'][0]
        bank['type_templates'][1]['markers'] = ['1/3', '2/7']*8
        def weights(alpha, marker):
            self.assertEqual(alpha, target.low.ALPHA)
            return {a: np.full(17, marker) for a in range(1, 17)}
        with patch.object(target.low.fractional, 'fractional_dense_weights', return_value={}), \
                patch.object(target.cover.emptycut, 'fractional_cardinality_weights', return_value={}), \
                patch.object(target.joint, 'type_weights', side_effect=weights) as build:
            result = list(target.fixed_templates(bank))
        self.assertEqual(build.call_count, 2)
        self.assertEqual(result[1][0]['markers'], ['1/3', '2/7']*8)
        np.testing.assert_array_equal(result[1][1][1], np.full(17, 1/3))
        np.testing.assert_array_equal(result[1][1][2], np.full(17, 2/7))

    def test_pinned_source_schema_and_partition(self):
        source = fixture()
        raw = json.dumps(source).encode()
        with patch.object(target, 'PINS', {}), patch.object(target, 'BANK_SHA', sha256(raw).hexdigest()), \
                patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target.low, 'checked_dependencies', return_value={}):
            self.assertEqual(len(target.load_source()['banks']), 3)
        source['banks'][0]['source_strata'][1]['count'] += 1
        raw = json.dumps(source).encode()
        with patch.object(target, 'PINS', {}), patch.object(target, 'BANK_SHA', sha256(raw).hexdigest()), \
                patch.object(target.Path, 'read_bytes', return_value=raw), \
                patch.object(target.low, 'checked_dependencies', return_value={}):
            with self.assertRaises(ValueError):
                target.load_source()

    def test_resource_guard_runs_before_source_access(self):
        with patch.object(target, 'load_source') as load:
            for budget in (0, -1, 31, True, float('nan'), float('inf')):
                with self.assertRaises(ValueError):
                    target.run(budget)
        load.assert_not_called()


if __name__ == '__main__':
    unittest.main()
