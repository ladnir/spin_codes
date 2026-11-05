"""Mocked driver/scope tests: never run a bank or production contraction."""
from copy import deepcopy
from unittest.mock import patch
import unittest

import numpy as np

import hd_two_sum_gf16_mean_zero as gate


def fixture():
    moments = [-2-m for m in range(19)]
    moments[0] = 0
    curves = [moments, [None]*19, [None]*19, [None]*19]
    return dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', alpha='1', regional_alpha=1, all_sources=False,
        all_sources_in_declared_stratum=True, all_fixed_arity_rowset_vectors=True,
        rowsets_resampled=False, source_count='17', initial_state=0,
        terminal_vector=[1, 1], regional_matrix_count=256,
        matrix_entry_order=['00', '01', '10', '11'], order_zero_moment=1,
        aggregate_input_weight_floor=1296, raw_regional_banks=1,
        new_moment_candidates=0, parameter_fits=0, directed=False,
        whole_code_certificate=False, complete_requested_control=True,
        regional_full_pair_restriction='f<=min(1,floor(m/2))',
        source_sha256=gate.SOURCE_PINS['hd_two_sum_late_alpha_live.py'],
        regional_log2_curves=curves,
        two_atom_envelope=dict(secant_endpoints=[5, 18], moment_order_zero=1,
                              log2_probability=-2., log2_nonzero_atom=-1.))


class MeanZeroTests(unittest.TestCase):
    def test_rechecks_raw_envelope_and_zero_path(self):
        matrix, log_p, log_b = gate.checked_inputs(fixture(), 17)
        np.testing.assert_array_equal(matrix, [[0., -np.inf], [-np.inf, -np.inf]])
        self.assertEqual((log_p, log_b), (-2., -1.))
        changed = fixture()
        changed['regional_log2_curves'][0][7] = 100
        with self.assertRaises(ArithmeticError):
            gate.checked_inputs(changed, 17)

    def test_scope_guards(self):
        for key, value in (('source_full_pair_max', 2), ('alpha', '.2'),
                           ('source_count', '18'), ('whole_code_certificate', True)):
            changed = deepcopy(fixture())
            changed[key] = value
            with self.assertRaises(ValueError):
                gate.checked_inputs(changed, 17)

    def test_exactly_two_controls_and_no_complete_class_claim(self):
        data = fixture()
        inputs = data, 17, gate.checked_inputs(data, 17)
        with patch.object(gate, 'read_inputs', return_value=inputs), \
             patch.object(gate, 'raw_holder_weights', return_value=np.zeros(17)) as weights, \
             patch.object(gate, 'whole_word_screen', side_effect=[
                 {'log2_screen': -400.}, {'log2_screen': -200.}]) as screen:
            result = gate.run()
        self.assertTrue(result['complete_requested_controls'])
        self.assertFalse(result['complete_homogeneous_class'])
        self.assertFalse(result['all_regional_state_paths_included'])
        self.assertEqual(result['raw_regional_banks_built'], 0)
        self.assertEqual([call.kwargs['quiet'] for call in weights.call_args_list], [True, False])
        self.assertEqual(screen.call_count, 2)

    def test_partial_result_preserved_and_budget_guards(self):
        data = fixture()
        with patch.object(gate, 'read_inputs', return_value=(data, 17, gate.checked_inputs(data, 17))), \
             patch.object(gate, 'raw_holder_weights', return_value=np.zeros(17)), \
             patch.object(gate, 'whole_word_screen', side_effect=[
                 {'log2_screen': -400.}, TimeoutError('toy')]):
            result = gate.run()
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(len(result['controls']), 1)
        for invalid in (True, 0, 31, float('nan')):
            with self.assertRaises(ValueError):
                gate.run(max_seconds=invalid)


if __name__ == '__main__':
    unittest.main()
