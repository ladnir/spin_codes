"""Completeness and mocked all-label wiring; no actual all-a3 run."""
from math import log2
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_all_a3 as target


def row(source, root=2.):
    return dict(original_source=source, arity=3, status='complete', aggregate_original_input_floor=486,
        normalized_group_root_log2=root, word_tail_log2=18*(root+24*log2(255)))


class AllThreeTests(unittest.TestCase):
    def test_complete_unique560_and_one_prefix(self):
        data = target.coverage([row(source) for source in target.SOURCES], 17.)
        self.assertTrue(data['complete'])
        self.assertEqual(data['unique_original_rowsets'], 560)
        self.assertAlmostEqual(data['all_a3_bound_log2'], 17+18*(2+log2(560)), places=10)

    def test_missing_label_never_gets_full_class_bound(self):
        data = target.coverage([row(source) for source in target.SOURCES[:-1]], 17.)
        self.assertFalse(data['complete'])
        self.assertEqual(data['missing_original_rowsets'], 1)
        self.assertIsNone(data['all_a3_bound_log2'])

    def test_duplicates_foreign_labels_and_bad_normalization_rejected(self):
        for rows in ([row(target.SOURCES[0])]*2, [row((0, 1, 16))],
                     [dict(row((0, 1, 2)), normalized_group_root_log2=9.)]):
            with self.assertRaises(ValueError): target.coverage(rows, 17.)

    def test_mock_all560_one_catalog_and_one_table(self):
        bank = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        catalog = SimpleNamespace(arity=3, q=256, counts=range(66795), positive_mass_kernels=564357)
        progress, seen = [], []
        def prepare(source, *args, **kwargs):
            seen.append(source)
            return np.zeros((8, 7, 2, 2)), dict(merged_joint_terms_per_pair=[3]*8)
        with patch.object(target.core.frozen, 'load_measure', return_value=(bank, {}, 17.)), \
             patch.object(target.core.catalog_api, 'load_catalog', return_value=catalog) as load, \
             patch.object(target.core.frozen, 'mask_log_matrices', return_value=None) as table, \
             patch.object(target.core, 'prepare_source', side_effect=prepare), \
             patch.object(target.core.tail_backend, 'positive_deficit_tail', return_value=18*(2+24*log2(255))) as tail:
            data = target.run(progress=progress.append)
        self.assertEqual(tuple(seen), target.SOURCES)
        self.assertEqual(load.call_count, 1)
        self.assertEqual(table.call_count, 1)
        self.assertEqual(tail.call_args.args[1:], (18, 486))
        self.assertTrue(data['complete_all560_original_a3'])
        self.assertAlmostEqual(data['all_a3_bound_log2'], 17+18*(2+log2(560)), places=10)
        self.assertEqual([x['completed'] for x in progress], list(range(32, 560, 32))+[560])
        self.assertFalse(data['extra_rowset_factor_applied'])
        self.assertEqual(data['catalog_builds'], 0)

    def test_timeout_keeps_completed_rows_only(self):
        bank = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        catalog = SimpleNamespace(arity=3, q=256, counts=range(66795), positive_mass_kernels=564357)
        calls = 0
        def tail(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls > 2: raise TimeoutError('mock budget')
            return 0.
        with patch.object(target.core.frozen, 'load_measure', return_value=(bank, {}, 17.)), \
             patch.object(target.core.catalog_api, 'load_catalog', return_value=catalog), \
             patch.object(target.core.frozen, 'mask_log_matrices', return_value=None), \
             patch.object(target.core, 'prepare_source', return_value=(np.zeros((8, 7, 2, 2)), dict(merged_joint_terms_per_pair=[3]*8))), \
             patch.object(target.core.tail_backend, 'positive_deficit_tail', side_effect=tail):
            data = target.run()
        self.assertFalse(data['complete'])
        self.assertEqual(data['unique_original_rowsets'], 2)
        self.assertEqual(len(data['rows']), 2)
        self.assertIsNone(data['all_a3_bound_log2'])
        self.assertTrue(data['partial_rows_are_not_full_class_bound'])


if __name__ == '__main__':
    unittest.main()
