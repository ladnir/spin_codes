from copy import deepcopy
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_post_refresh_moment_gate as gate


def fixture():
    metadata = dict(construction=gate.CONSTRUCTION, state_order=list(gate.STATE_ORDER),
        regional_alpha=1, total_packets=128, total_pairs=4096, global_normalization_count=1,
        coefficient_maximum=18, coefficient_full_max=1, orientation_averaged_in_packet=True,
        source_geometry_resampled=False, states_persist_across_regions=True,
        cone_P_arbitrary_diagonal_probability=True, cone_G_arbitrary_first_nonzero_probability=True,
        cone_F_conditional_second_flat_density=True, probability_rows_assumed=False,
        directed=False, whole_code_certificate=False)
    # Raw identity caps test strict serialization only, not a new operator.
    rows = [dict(active_weight=m, full_pairs=f,
        log2_matrix=[[0. if i == j else None for j in range(8)] for i in range(8)])
        for m in range(19) for f in range(min(1, m//2)+1)]
    return dict(construction=gate.CONSTRUCTION, mark='691/700', maximum_regional_occupancy=18,
        full_pair_maximum=1, regional_alpha=1, state_order=list(gate.STATE_ORDER), initial_state='Z',
        terminal_vector=[1]*8, regional_histogram_count=36, zero_matrix_entries_encoded_as_none=True,
        raw_packet_banks_built=1, parameter_fits=0, complete_requested_bank=True,
        directed=False, whole_code_certificate=False, source_sha256=gate.PROVIDER_SHA,
        reviewed_source_sha256=gate.PROVIDER_DEPENDENCIES, provider_metadata=metadata,
        regional_log2_matrices=rows)


class PostRefreshMomentGateTests(TestCase):
    def test_complete_saved_bank_and_exact_zero_round_trip(self):
        bank = gate.decode_saved_bank(fixture())
        self.assertEqual(len(bank), 36)
        np.testing.assert_array_equal(np.exp2(bank[18, 1]), np.eye(8))
        self.assertTrue(np.isneginf(bank[0, 0][0, 1]))
        raw = json.dumps(fixture(), allow_nan=False).encode()
        with patch.object(Path, 'read_bytes', return_value=raw):
            _, again = gate.load_saved_bank('not-created.json', sha256(raw).hexdigest())
            self.assertEqual(set(again), set(bank))
            with self.assertRaises(ValueError):
                gate.load_saved_bank('not-created.json', '0'*64)

    def test_scope_provenance_and_histograms_fail_closed(self):
        for key, value in (('regional_alpha', .2), ('full_pair_maximum', 2),
            ('complete_requested_bank', False), ('construction', 'old-construction'),
            ('state_order', list(reversed(gate.STATE_ORDER))), ('source_sha256', '0'*64)):
            data = fixture()
            data[key] = value
            with self.assertRaises(ValueError):
                gate.decode_saved_bank(data)
        data = fixture()
        data['provider_metadata']['source_geometry_resampled'] = True
        with self.assertRaises(ValueError):
            gate.decode_saved_bank(data)
        for mutation in ('missing', 'duplicate', 'wrongshape'):
            data = fixture()
            if mutation == 'missing':
                data['regional_log2_matrices'].pop()
            elif mutation == 'duplicate':
                data['regional_log2_matrices'].append(deepcopy(data['regional_log2_matrices'][0]))
            else:
                data['regional_log2_matrices'][0]['log2_matrix'] = [[0.]]
            with self.assertRaises(ValueError):
                gate.decode_saved_bank(data)

    def test_mocked_gate_fixed_weights_prefix_and_no_provider(self):
        data = fixture()
        bank = gate.decode_saved_bank(data)
        curve = np.array([[0.]+[-20.-5*m for m in range(1, 19)],
            [-np.inf]+[-4.-m for m in range(1, 19)],
            [-np.inf]+[-130.-2*m for m in range(1, 19)], [-120.]*19])
        def control(weights, measure, prefix, **kwargs):
            self.assertEqual(measure['zero_secant_endpoints'], [5, 18])
            self.assertEqual(measure['cross_secant_endpoints'], [2, 18])
            self.assertEqual(measure['log2_atom01'], -1.)
            self.assertEqual(measure['log2_atom10'], -2.)
            return dict(log2_screen=prefix-5000., positive_tail_log2=-5000.)
        with patch.object(gate, 'check_source_pins'), \
             patch.object(gate, 'load_saved_bank', return_value=(data, bank)), \
             patch.object(gate, 'load_source_count', return_value=({'source_count': str(gate.SOURCE_COUNT)}, None)), \
             patch.object(gate, 'entrywise_curves', return_value=curve), \
             patch.object(gate, 'raw_holder_weights', return_value=np.zeros(17)), \
             patch.object(gate, 'polynomial_control', side_effect=control) as evaluate:
            result = gate.run('not-created.json', '0'*64)
        evaluate.assert_called_once()
        expected = log2(gate.SOURCE_COUNT)-gate.LIMIT*log2(float(gate.MARK))-1153*log2(255)
        self.assertAlmostEqual(result['prefix_log2'], expected, places=10)
        self.assertEqual(result['fixed_live_log_weights'], [0., -128., 0., 0., -128., -128., -256.])
        self.assertEqual(result['raw_regional_banks_built'], 0)
        self.assertEqual(result['extra_terminal_factor_log2'], 0.)
        self.assertEqual(result['rs_normalizer_exponent'], 1152)
        self.assertEqual(result['aggregate_surrogate_floor'], 1296)
        self.assertFalse(result['extra_rowset_factor'])
        self.assertTrue(result['complete_requested_control'])
        self.assertFalse(result['whole_code_certificate'])

    def test_failed_moment_family_is_incomplete_not_a_bound(self):
        data = fixture()
        with patch.object(gate, 'check_source_pins'), \
             patch.object(gate, 'load_saved_bank', return_value=(data, gate.decode_saved_bank(data))), \
             patch.object(gate, 'load_source_count', return_value=({'source_count': str(gate.SOURCE_COUNT)}, None)), \
             patch.object(gate, 'entry_measure', side_effect=gate.MomentFitError('prescribed p>1')), \
             patch.object(gate, 'polynomial_control') as evaluate:
            result = gate.run('not-created.json', '0'*64)
        evaluate.assert_not_called()
        self.assertFalse(result['complete_requested_control'])
        self.assertEqual(result['incomplete_kind'], 'MomentFitError')
        self.assertNotIn('log2_screen', result)

    def test_budget_hash_and_shared_engine_constants(self):
        with self.assertRaises(ValueError):
            gate.run('unused', '0'*64, max_seconds=31)
        with self.assertRaises(ValueError):
            gate.load_saved_bank('unused', 'not-a-sha')
        import hd_coupled_state_moment_gate as engine
        self.assertEqual((gate.GROUPS, gate.ARITY, gate.COLUMNS, gate.FLOOR),
                         (engine.GROUPS, engine.ARITY, engine.COLUMNS, engine.FLOOR))
        gate.check_source_pins()  # Read/hash only; no bank or gate evaluation.


if __name__ == '__main__':
    main()
