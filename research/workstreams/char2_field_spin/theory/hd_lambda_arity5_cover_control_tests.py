"""Mocked controller/coverage checks only: no full row-set proof run."""
from math import log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_arity5_cover_control as cover


def record(source, marked=-30., exact=None):
    result = dict(original_source=list(source), original_mask=cover.source_mask(source),
                  marked_normalized_root_log2=marked,
                  normalized_fixed_A_root_log2=min(marked, exact) if exact is not None else marked)
    if exact is not None:
        result['exact_normalized_root_log2'] = exact
    return result


class ArityFiveCoverTests(unittest.TestCase):
    def test_complete_exact_inventory_and_root_sum(self):
        sources = cover.original_sources()
        self.assertEqual(len(sources), 4368)
        self.assertEqual(len({cover.source_mask(source) for source in sources}), 4368)
        records = [record(source) for source in sources]
        summary = cover.coverage_summary(records, 17.)
        self.assertTrue(summary['all_original_A_covered'])
        self.assertEqual(summary['missing_original_masks'], [])
        self.assertAlmostEqual(summary['full_selected_arity_bound_log2'],
                               17+18*(log2(4368)-30), places=10)

    def test_missing_duplicate_and_wrong_minimum_cannot_claim_complete(self):
        rows = [record(source) for source in cover.original_sources()]
        self.assertFalse(cover.coverage_summary(rows[:-1], 0.)['all_original_A_covered'])
        duplicate = cover.coverage_summary(rows[:-1]+[rows[0]], 0.)
        self.assertFalse(duplicate['all_original_A_covered'])
        self.assertTrue(duplicate['duplicate_original_masks'])
        self.assertIsNone(duplicate['covered_subset_bound_log2'])
        bad = dict(rows[0], exact_normalized_root_log2=-40.)
        invalid = cover.coverage_summary([bad]+rows[1:], 0.)
        self.assertFalse(invalid['all_original_A_covered'])
        self.assertTrue(invalid['invalid_records'])

    def test_batched_axes_fallback_cap_and_complete_tail_minimum(self):
        sources = cover.original_sources()[:3]
        prepare_calls, tail_shapes = [], []
        def prepare(source, *unused, **kwargs):
            prepare_calls.append(source)
            return np.full((16, 6, 2, 2), 200.), []
        def exact(source, *unused, **kwargs):
            return np.full((16, 6, 2, 2), 10.), {'toy': True}
        def tail(columns, arity, groups, **kwargs):
            columns = np.asarray(columns)
            tail_shapes.append(columns.shape)
            if columns.ndim == 5:
                return columns[0, :, 0, 0, 0]
            return columns[0, 0, 0, 0]
        with patch.object(cover, 'prepare_source', side_effect=prepare), \
             patch.object(cover, 'exact_columns', side_effect=exact), \
             patch.object(cover, 'column_deficit_tail', side_effect=tail):
            # Large prefix ensures every marked fixed-A bound requests fallback.
            result = cover._process_sources(sources, None, None, 10000., batch_size=2,
                fallback_limit=1, check_time=lambda: None)
        self.assertEqual(prepare_calls, list(sources))
        self.assertEqual(tail_shapes, [(16, 2, 6, 2, 2), (16, 6, 2, 2), (16, 1, 6, 2, 2)])
        self.assertEqual(len(result['exact_fallbacks']), 1)
        self.assertEqual(result['exact_fallback_cap_skipped'], 2)
        first = result['records'][0]
        self.assertEqual(first['selected_method'], 'exact-GF4')
        self.assertAlmostEqual(first['normalized_fixed_A_root_log2'], 10/18-40*log2(255))
        self.assertTrue(all(row['selected_method'] == 'marked-fibre' for row in result['records'][1:]))

    def test_timeout_during_exact_fallback_retains_valid_marked_batch(self):
        sources = cover.original_sources()[:2]
        with patch.object(cover, 'prepare_source', return_value=(np.zeros((16, 6, 2, 2)), [])), \
             patch.object(cover, 'column_deficit_tail', return_value=np.array((200., 210.))), \
             patch.object(cover, 'exact_columns', side_effect=TimeoutError('toy budget')):
            result = cover._process_sources(sources, None, None, 10000., batch_size=2,
                fallback_limit=128, check_time=lambda: None)
        self.assertEqual(len(result['records']), 2)
        self.assertEqual(result['exact_fallbacks'][0]['status'], 'timeout_marked_retained')
        self.assertEqual(result['stop_reason']['type'], 'TimeoutError')
        self.assertEqual(result['records'][0]['selected_method'], 'marked-fibre')
        self.assertNotIn('exact_normalized_root_log2', result['records'][0])

    def test_source_pin_change_withholds_full_result(self):
        sources = cover.original_sources()
        rows = [record(source) for source in sources]
        work = dict(records=rows, exact_fallbacks=[], exact_fallback_cap_skipped=0,
                    pending_preparations_without_tail=0, stop_reason=None)
        saved = dict(mark='9927/10000', source_count=123)
        with patch.object(cover, '_pin_files', side_effect=({'file': 'before'}, {'file': 'after'})), \
             patch.object(cover, '_input_hash', return_value=cover.INPUT_SHA), \
             patch.object(cover, 'load_measure', return_value=(saved, {}, 0.)), \
             patch.object(cover, 'moment_table', return_value=None), \
             patch.object(cover, '_process_sources', return_value=work):
            result = cover.run()
        self.assertTrue(result['coverage']['all_original_A_covered'])
        self.assertFalse(result['complete'])
        self.assertFalse(result['target_met'])
        self.assertFalse(result['source_pins_unchanged'])
        self.assertIsNone(result['coverage']['full_selected_arity_bound_log2'])

    def test_threshold_and_guard_scope(self):
        self.assertAlmostEqual(cover.FALLBACK_THRESHOLD+18*log2(4368), -60.)
        for kwargs in ({'batch_size': 64}, {'max_seconds': 361}, {'exact_fallback_limit': 129}):
            with self.assertRaises(ValueError):
                cover.run(**kwargs)
        with self.assertRaises(ValueError):
            cover.source_mask((0, 1, 2, 3, 3))
        self.assertIsNone(cover._jsonable(float('-inf')))


if __name__ == '__main__':
    unittest.main()
