"""Only tiny/mock controller checks; geometry owns helper's exact toy proofs."""
from math import log2
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_shared_control as target


class SharedRatioControlTests(unittest.TestCase):
    def test_preserves_returned_root_polynomials_and_pair_chronology(self):
        prepared = SimpleNamespace(groups=18, rows=16)
        seen = []
        def pair(source, column, same_prepared, **kwargs):
            self.assertEqual(source, (0, 1, 4, 8))
            self.assertIs(same_prepared, prepared)
            self.assertEqual(kwargs['q'], 256)
            seen.append(column)
            return np.full((9, 2, 2), float(36+column))
        with patch.object(target.ratio, 'pair_polynomial', side_effect=pair):
            data = target.prepare_source((0, 1, 4, 8), prepared)
        self.assertEqual(seen, list(range(0, 16, 2)))
        self.assertTrue(np.array_equal(data[:, 0, 0, 0], np.arange(36., 52., 2.)))
        self.assertEqual(data.shape, (8, 9, 2, 2))

    def test_mock_normalization_single_measure_and_prefix(self):
        bank = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        prepared = SimpleNamespace(groups=18, rows=16)
        with patch.object(target.frozen, 'load_measure', return_value=(bank, {}, 17.)), \
             patch.object(target.ratio, 'prepare_measure', return_value=prepared) as once, \
             patch.object(target, 'prepare_source', return_value=np.zeros((8, 9, 2, 2))), \
             patch.object(target.tail_backend, 'positive_deficit_tail', return_value=18*(32*log2(255)+2)) as tail:
            data = target.run(sources=target.SOURCES[1:])
        self.assertEqual(once.call_count, 1)
        self.assertEqual(tail.call_args.args[1:], (18, 648))
        self.assertTrue(data['complete'])
        self.assertAlmostEqual(data['selected_rowset_sum_bound_log2'], 17+18*(2+1))
        self.assertEqual(data['rows'][0]['RS_factor_whole_word'], '255^(-576)')
        self.assertFalse(data['helper_result_rooted_again'])
        self.assertFalse(data['full_a4_class_covered'])

    def test_optional_reference_may_be_above_or_below_new_bound(self):
        bank = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        for old_bound in (-999., 999.):
            with patch.object(target.frozen, 'load_measure', return_value=(bank, {}, 17.)), \
                 patch.object(target, 'load_reference', return_value={(0, 1, 2): dict(fixed_A_repeated18_bound_log2=old_bound)}), \
                 patch.object(target.ratio, 'prepare_measure', return_value=SimpleNamespace(groups=18, rows=16)), \
                 patch.object(target, 'prepare_source', return_value=np.zeros((8, 7, 2, 2))), \
                 patch.object(target.tail_backend, 'positive_deficit_tail', return_value=18*(24*log2(255)+2)):
                data = target.run(sources=((0, 1, 2),))
            self.assertAlmostEqual(data['rows'][0]['new_minus_prior_bound_log2'], 53-old_bound)
            self.assertFalse(data['rows'][0]['ordering_between_two_valid_majorants_asserted'])

    def test_timeout_keeps_completed_fixedA_but_no_requested_set_sum(self):
        bank = dict(groups=18, source_class='Fmax<=1', source_count='7', mark='9927/10000')
        calls = 0
        def tail(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls > 1: raise TimeoutError('mock budget')
            return 0.
        with patch.object(target.frozen, 'load_measure', return_value=(bank, {}, 17.)), \
             patch.object(target.ratio, 'prepare_measure', return_value=SimpleNamespace(groups=18, rows=16)), \
             patch.object(target, 'prepare_source', return_value=np.zeros((8, 9, 2, 2))), \
             patch.object(target.tail_backend, 'positive_deficit_tail', side_effect=tail):
            data = target.run(sources=target.SOURCES[1:])
        self.assertFalse(data['complete'])
        self.assertEqual(data['completed_sources'], 1)
        self.assertIsNone(data['selected_rowset_sum_bound_log2'])

    def test_scope_guards(self):
        for sources in ((), ((0, 1, 2),)*2, ((0, 1, 2, 3, 4),), ((0, 2, 1),),
                        tuple((0, 1, i) for i in range(2, 11))):
            with self.assertRaises(ValueError): target.validate_sources(sources)
        with self.assertRaises(ValueError): target.run(max_seconds=121)


if __name__ == '__main__':
    unittest.main()
