"""Small schema/mock tests only: no real a4 table, tail, or batch run."""
from copy import deepcopy
from math import log2
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_a4_batches as target


PINS = {'reviewed.py': 'test-only'}
BANK = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
CONTEXT = target.common_context(BANK, {'test-measure': 1}, 17.)


def record(ordinal, root=2.):
    return dict(ordinal=ordinal, original_source=list(target.SOURCES[ordinal]), arity=4,
        status='complete', aggregate_original_input_floor=648,
        word_tail_log2=18*(root+32*log2(255)), normalized_group_root_log2=root,
        fixed_A_repeated18_bound_log2=17+18*root)


def receipt(ordinals, completed=None):
    ordinals = list(ordinals)
    completed = len(ordinals) if completed is None else completed
    return dict(schema=target.SCHEMA, kind='batch', **deepcopy(CONTEXT),
        requested_ordinals=ordinals, requested_original_sources=[list(target.SOURCES[i]) for i in ordinals],
        rows=[record(i) for i in ordinals[:completed]], completed_sources=completed,
        batch_complete=completed == len(ordinals), full_a4_class_covered=False, all_a4_bound_log2=None,
        reviewed_source_sha256=dict(PINS), source_pins_unchanged=True)


class BatchTests(unittest.TestCase):
    def test_interval_and_explicit_list_scope(self):
        self.assertEqual(target.select_ordinals(begin=0, end=384), tuple(range(384)))
        self.assertEqual(target.select_ordinals(ordinals=(2, 7, 1819)), (2, 7, 1819))
        for kwargs in (dict(begin=0, end=385), dict(begin=-1, end=2), dict(begin=1, end=1821),
                       dict(ordinals=()), dict(ordinals=(1, 1)), dict(ordinals=(2, 1)),
                       dict(ordinals=(False, 1)), dict(begin=0, end=1, ordinals=(0,))):
            with self.assertRaises(ValueError): target.select_ordinals(**kwargs)

    def test_mock_batch_keeps_original_labels_and_single_preparation(self):
        prepared = SimpleNamespace(groups=18, rows=16)
        seen = []
        def prepare(source, same_prepared, **kwargs):
            self.assertIs(same_prepared, prepared)
            seen.append(source)
            return np.zeros((8, 9, 2, 2))
        with patch.object(target, 'dependencies', return_value=PINS), \
             patch.object(target.core.frozen, 'load_measure', return_value=(BANK, {'test-measure': 1}, 17.)), \
             patch.object(target.core.ratio, 'prepare_measure', return_value=prepared) as once, \
             patch.object(target.core, 'prepare_source', side_effect=prepare), \
             patch.object(target.core.tail_backend, 'positive_deficit_tail', return_value=18*(2+32*log2(255))) as tail:
            data = target.run_batch(ordinals=(2, 7, 1819))
        self.assertEqual(seen, [target.SOURCES[i] for i in (2, 7, 1819)])
        self.assertEqual(once.call_count, 1)
        self.assertEqual(tail.call_args.args[1:], (18, 648))
        self.assertTrue(data['batch_complete'])
        self.assertFalse(data['contiguous_ordinals'])
        self.assertIsNone(data['all_a4_bound_log2'])
        self.assertAlmostEqual(data['rows'][0]['normalized_group_root_log2'], 2.)
        self.assertFalse(data['helper_result_rooted_again'])

    def test_timeout_prefix_and_resume_only_uncompleted_labels(self):
        def tail(*args, **kwargs):
            if tail.called: raise TimeoutError('mock budget')
            tail.called = True
            return 0.
        tail.called = False
        with patch.object(target, 'dependencies', return_value=PINS), \
             patch.object(target.core.frozen, 'load_measure', return_value=(BANK, {'test-measure': 1}, 17.)), \
             patch.object(target.core.ratio, 'prepare_measure', return_value=SimpleNamespace(groups=18, rows=16)), \
             patch.object(target.core, 'prepare_source', return_value=np.zeros((8, 9, 2, 2))), \
             patch.object(target.core.tail_backend, 'positive_deficit_tail', side_effect=tail):
            data = target.run_batch(begin=10, end=13)
        self.assertFalse(data['batch_complete'])
        self.assertEqual(data['completed_sources'], 1)
        self.assertEqual(data['timeout'], 'mock budget')
        assembled = target.assemble_receipts([data, receipt((11, 12))], CONTEXT, PINS)
        self.assertEqual(assembled['unique_original_rowsets'], 3)
        self.assertFalse(assembled['complete'])
        self.assertIsNone(assembled['all_a4_bound_log2'])

    def test_complete_1820_labels_prefix_once_without_extra_factor(self):
        batches = [receipt(range(begin, min(begin+384, 1820))) for begin in range(0, 1820, 384)]
        data = target.assemble_receipts(batches[::-1], CONTEXT, PINS)
        self.assertTrue(data['full_a4_class_covered'])
        self.assertEqual(data['unique_original_rowsets'], 1820)
        self.assertEqual([row['ordinal'] for row in data['rows']], list(range(1820)))
        self.assertEqual(data['missing_ordinals'], [])
        self.assertAlmostEqual(data['all_a4_bound_log2'], 17+18*(2+log2(1820)), places=10)
        self.assertFalse(data['extra_rowset_factor_applied'])

    def test_gap_and_completed_overlap(self):
        batches = [receipt(range(begin, min(begin+384, 1819))) for begin in range(0, 1819, 384)]
        data = target.assemble_receipts(batches, CONTEXT, PINS)
        self.assertFalse(data['complete'])
        self.assertEqual(data['missing_ordinals'], [1819])
        self.assertIsNone(data['all_a4_bound_log2'])
        with self.assertRaises(ValueError):
            target.assemble_receipts([receipt((0, 1)), receipt((1, 2))], CONTEXT, PINS)

    def test_reject_source_pins_measure_prefix_and_false_claim(self):
        changes = [dict(reviewed_source_sha256={}), dict(source_pins_unchanged=False),
                   dict(common_measure={}), dict(common_prefix_log2=18.), dict(source_count='8'),
                   dict(full_a4_class_covered=True), dict(all_a4_bound_log2=-100.),
                   dict(batch_complete=False), dict(helper_result_rooted_again=True)]
        for delta in changes:
            data = receipt((0, 1))
            data.update(delta)
            with self.assertRaises(ValueError): target.validate_batch(data, CONTEXT, PINS)

    def test_reject_label_permutation_bad_normalization_nonfinite_and_missing_row(self):
        for field, value in (('original_source', [0, 1, 3, 2]), ('aggregate_original_input_floor', 486),
                             ('normalized_group_root_log2', 0.), ('word_tail_log2', float('nan')),
                             ('fixed_A_repeated18_bound_log2', 1.), ('ordinal', True)):
            data = receipt((0, 1))
            data['rows'][0][field] = value
            with self.assertRaises(ValueError): target.validate_batch(data, CONTEXT, PINS)
        data = receipt((0, 1, 2))
        data['rows'] = [record(0), record(2)]
        data['completed_sources'], data['batch_complete'] = 2, False
        with self.assertRaises(ValueError): target.validate_batch(data, CONTEXT, PINS)

    def test_source_change_and_budget_limits_rejected(self):
        with self.assertRaises(ValueError): target.run_batch(begin=0, end=1, max_seconds=181)
        with patch.object(target, 'dependencies', side_effect=[PINS, {}]), \
             patch.object(target.core.frozen, 'load_measure', return_value=(BANK, {'test-measure': 1}, 17.)), \
             patch.object(target.core.ratio, 'prepare_measure', return_value=SimpleNamespace(groups=18, rows=16)), \
             patch.object(target.core, 'prepare_source', return_value=np.zeros((8, 9, 2, 2))), \
             patch.object(target.core.tail_backend, 'positive_deficit_tail', return_value=0.):
            with self.assertRaises(ValueError): target.run_batch(begin=0, end=1)


if __name__ == '__main__':
    unittest.main()
