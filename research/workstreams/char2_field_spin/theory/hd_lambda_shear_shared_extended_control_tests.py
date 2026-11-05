"""Tiny tail/replay and mocked controller tests; no production catalogs/bank."""
from math import log2
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_shared_extended_control as control


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, dtype=float))


def direct_tail(columns, groups, threshold):
    word = np.eye(2)[None]
    for column in columns:
        powered = np.zeros((1, 2, 2))
        powered[0] = 1
        for _ in range(groups):
            new = np.zeros((len(powered)+len(column)-1, 2, 2))
            for i, left in enumerate(powered):
                for j, right in enumerate(column):
                    new[i+j] += left*right
            powered = new
        new = np.zeros((len(word)+len(powered)-1, 2, 2))
        for i, left in enumerate(word):
            for j, right in enumerate(powered):
                new[i+j] += left @ right
        word = new
    return word[threshold:, 0, :].sum()


class ExtendedControlTests(unittest.TestCase):
    def test_new_degree12_tail_against_positive_noncommuting_enumeration(self):
        columns = np.zeros((2, 13, 2, 2))
        columns[0, 0] = [[1, 2], [0, 1]]
        columns[0, 7] = [[0, 1], [3, 1]]
        columns[0, 12] = [[2, 0], [1, 1]]
        columns[1, 0] = [[1, 0], [2, 1]]
        columns[1, 5] = [[1, 3], [0, 2]]
        columns[1, 12] = [[0, 2], [1, 1]]
        for groups, threshold in ((1, 0), (1, 19), (2, 27), (2, 48)):
            got = float(control.positive_pair_deficit_tail(logs(columns), groups, threshold))
            expected = direct_tail(columns, groups, threshold)
            self.assertAlmostEqual(got, log2(expected), places=12)

    def test_overlap_with_frozen_degree11_dispatcher(self):
        columns = logs(np.arange(1, 1+3*5*4).reshape(3, 5, 2, 2)/20)
        for threshold in (0, 9, 24):
            actual = control.positive_pair_deficit_tail(columns, 2, threshold)
            expected = control.generic.positive_deficit_tail(columns, 2, threshold)
            np.testing.assert_array_equal(actual, expected)

    def test_tail_guards_and_callback(self):
        valid = np.zeros((8, 13, 2, 2))
        for groups, floor in ((0, 0), (19, 0), (1, 97), (1, -1)):
            with self.assertRaises(ValueError):
                control.positive_pair_deficit_tail(valid, groups, floor)
        with self.assertRaises(ValueError):
            control.positive_pair_deficit_tail(np.zeros((8, 14, 2, 2)), 1, 0)
        with self.assertRaises(ValueError):
            control.positive_pair_deficit_tail(np.zeros((9, 2, 2, 2)), 1, 0)
        with self.assertRaises(TimeoutError):
            control.positive_pair_deficit_tail(valid, 1, 0,
                check_time=lambda: (_ for _ in ()).throw(TimeoutError('test')))

    def test_controller_normalization_with_mocked_no_job_interfaces(self):
        bank = dict(groups=18, source_class='Fmax<=1', mark='9927/10000', source_count='123')
        sources = (tuple(range(5)), tuple(range(6)))
        builder = lambda a, **kw: SimpleNamespace(arity=a, q=256, positive_mass_kernels=3, counts=(1, 2))
        with patch.object(control, 'dependencies', return_value={'mock': 'pin'}), \
             patch.object(control.frozen, 'load_measure', return_value=(bank, {}, 1000.)), \
             patch.object(control.ratio, 'prepare_measure', return_value=SimpleNamespace(groups=18, rows=16)) as prepared, \
             patch.object(control.extended, 'build_catalog', side_effect=builder) as built, \
             patch.object(control, 'prepare_source', return_value=(np.zeros((8, 13, 2, 2)), [])), \
             patch.object(control, 'positive_pair_deficit_tail', side_effect=(-100., -200.)):
            result = control.run(sources=sources)
        self.assertTrue(result['complete'])
        self.assertEqual(prepared.call_count, 1)
        self.assertEqual(built.call_count, 2)
        expected_roots = np.array([-100/18-40*log2(255), -200/18-48*log2(255)])
        np.testing.assert_allclose([row['normalized_group_root_log2'] for row in result['rows']], expected_roots)
        self.assertAlmostEqual(result['selected_rowset_sum_bound_log2'],
            1000+18*float(np.logaddexp2.reduce(expected_roots)))
        self.assertFalse(result['all_labelled_class_covered'])
        self.assertFalse(result['whole_code_certificate'])

    def test_partial_budget_receipt_and_fixed_scope(self):
        bank = dict(groups=18, source_class='Fmax<=1', mark='9927/10000', source_count='123')
        catalog = SimpleNamespace(arity=5, q=256, positive_mass_kernels=3, counts=(1, 2))
        with patch.object(control, 'dependencies', return_value={}), \
             patch.object(control.frozen, 'load_measure', return_value=(bank, {}, 1000.)), \
             patch.object(control.ratio, 'prepare_measure', return_value=object()), \
             patch.object(control.extended, 'build_catalog', return_value=catalog), \
             patch.object(control, 'prepare_source', side_effect=control.extended.BudgetExceeded('over limit')):
            result = control.run(sources=(tuple(range(5)),))
        self.assertFalse(result['complete'])
        self.assertIsNone(result['selected_rowset_sum_bound_log2'])
        self.assertEqual(result['stop_reason'], 'over limit')
        self.assertEqual(result['rows'], [])
        for source in ((0, 1, 2, 3), (0, 1, 2, 3, 4, 5, 6), (0, 1, 2, 4, 3)):
            with self.assertRaises(ValueError):
                control.validate_sources((source,))

    def test_source_and_bank_pins_read_only(self):
        pins = control.dependencies()
        self.assertEqual(pins['hd_lambda_shear_shared_ratio_extended.py'],
            'b51f0869028327137f58c95d260a86653cc2b4acfb50275b6ab346cd8db0c01f')
        self.assertEqual(pins[control.frozen.INPUT], control.frozen.INPUT_SHA)


if __name__ == '__main__':
    unittest.main()
