"""Frozen witness A/B orchestration tests; no production provider calls."""
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_paired_low_control as target


class PairedLowControlTests(unittest.TestCase):
    def test_reference_pins_and_complete_forty_five_scope(self):
        data = target.load_reference()
        self.assertEqual(len(data['witnesses']), 45)
        self.assertEqual([row['active_groups'] for row in data['selected_anchor_bounds']], [18, 24, 32])
        self.assertFalse(data['continuous_weight_band_covered'])

    def test_three_banks_keep_every_original_numeric_witness(self):
        reference = target.load_reference()
        curves = [np.array([[(-np.inf if value is None else value-.125) for value in row]
                           for row in bank['curves']]) for bank in reference['banks']]
        seen, prepared = [], []
        def prepare(c, w, floors, slopes):
            prepared.append((list(map(str, floors)), list(slopes)))
            return floors, slopes
        def evaluate(pair, r, a0, mark, log_t):
            original = reference['witnesses'][len(seen)]
            self.assertEqual(pair[0], list(map(target.F, original['log2_floor_markers'])))
            self.assertEqual(list(pair[1]), original['log2_slopes'])
            self.assertEqual((r, a0, str(mark), log_t),
                             (original['active_groups'], original['high_arity_cutoff'], original['mark'],
                              original['log2_total_arity_marker']))
            seen.append(original)
            return dict(original, failure_log2_screen=original['failure_log2_screen']-1)
        with patch.object(target, 'load_reference', return_value=reference), \
                patch.object(target.low.fractional, 'fractional_dense_weights', return_value={}), \
                patch.object(target.paired, 'paired_nz_regions') as provider, \
                patch.object(target.low.state, 'regional_curves', side_effect=curves), \
                patch.object(target.low, 'prepared_envelope', side_effect=prepare), \
                patch.object(target.low, 'evaluate_prepared', side_effect=evaluate):
            result = target.run()
        self.assertEqual(provider.call_count, 3)
        self.assertEqual(len(seen), 45)
        self.assertEqual(len(prepared), 9)
        self.assertTrue(result['complete_requested_anchors'])
        self.assertTrue(result['identical_numeric_slopes_to_reference'])
        self.assertFalse(result['directed'])
        self.assertEqual(result['fitted_parameters'], 0)
        self.assertTrue(all(row['change_from_packetwise'] == -1 for row in result['witnesses']))

    def test_provider_regression_fails_before_any_control(self):
        reference = target.load_reference()
        curves = np.array([[(-np.inf if value is None else value+1) for value in row]
                           for row in reference['banks'][0]['curves']])
        with patch.object(target, 'load_reference', return_value=reference), \
                patch.object(target.low.fractional, 'fractional_dense_weights', return_value={}), \
                patch.object(target.paired, 'paired_nz_regions'), \
                patch.object(target.low.state, 'regional_curves', return_value=curves), \
                patch.object(target.low, 'evaluate_prepared') as evaluate, self.assertRaises(ArithmeticError):
            target.run()
        evaluate.assert_not_called()

    def test_guards_and_timeout_do_not_build_banks(self):
        for budget in (31, True, float('nan')):
            with patch.object(target.paired, 'paired_nz_regions') as provider, self.assertRaises(ValueError):
                target.run(max_seconds=budget)
            provider.assert_not_called()
        with patch.object(target, 'load_reference') as load:
            result = target.run(max_seconds=1e-12)
        load.assert_not_called()
        self.assertFalse(result['complete_requested_anchors'])


if __name__ == '__main__':
    unittest.main()
