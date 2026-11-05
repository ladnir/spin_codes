"""Tiny inverse/scalar and mocked cutoff tests; no q128 terminal gate."""

from fractions import Fraction as F
import unittest
from unittest.mock import patch

import hd_field128_full_rs_state_terminal_certificate as wrapper
from hd_blockwise_maps import matvec, mul


def mock_terminal(*, q, distance, first, failure_bits, precision):
    first = 7000 if first is None else first
    return dict(construction='old identity must be replaced', field_order=q, field_bits=q.bit_length()-1,
                weight_range=[first, 8192], target_distance=str(distance), witness={'mark': '1/256'},
                log2_upper_rational='-80', log2_upper_enclosure='-80',
                covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                no_support_subset_union=True, meets_requested_range_target=True,
                requested_failure_bits=failure_bits, directed=True, whole_code_certificate=False)


class StateTerminalTests(unittest.TestCase):
    def test_toy_identity_feedthrough_recovers_stream_and_state(self):
        q, state = 16, 3
        a = wrapper.local.lagrange(range(8), range(8, 11), q)
        c = tuple(zip(*a))
        h = wrapper.local.lagrange(range(state, 2*state), range(state), q)
        first, second = (2, 3, 5), (7, 9, 11)
        def advance(incoming, x):
            syndrome = matvec(c, x, q)
            refreshed = matvec(h, tuple(mul(d, s^v, q) for d, s, v in zip(first, incoming, syndrome)), q)
            return tuple(mul(d, v, q) for d, v in zip(second, refreshed))
        inputs = [tuple((5*i+3*t)%q for i in range(8)) for t in range(4)]
        incoming, outputs = (0,)*state, []
        for x in inputs:
            outputs.append(tuple(u^v for u, v in zip(x, matvec(a, incoming, q))))
            incoming = advance(incoming, x)
        recovered_state, recovered = (0,)*state, []
        for y in outputs:
            x = tuple(u^v for u, v in zip(y, matvec(a, recovered_state, q)))
            recovered.append(x)
            recovered_state = advance(recovered_state, x)
        self.assertEqual(inputs, recovered)
        self.assertEqual(incoming, recovered_state)

    def test_fresh_small_field_scalar_formula_matches_frozen_helper(self):
        # This explicit q256 tail need not pass; it checks formula/metadata only.
        actual = wrapper.terminal_receipt(q=256, state=24, distance=F(2, 5), first=8192)
        frozen = wrapper.terminal.terminal_receipt(q=256, distance=F(2, 5), first=8192)
        self.assertEqual(actual['log2_upper_rational'], frozen['log2_upper_rational'])
        self.assertEqual(actual['witness'], frozen['witness'])
        self.assertEqual(actual['construction'], wrapper.local.construction(24))
        self.assertEqual(actual['target_distance'], '2/5')
        self.assertTrue(actual['covers_all_weights_in_range_simultaneously'])
        self.assertTrue(actual['no_union_over_tail_weights'])
        self.assertFalse(actual['old_numerical_bounds_reused'])
        self.assertFalse(actual['whole_code_certificate'])

    def test_dynamic_cutoff_uses_requested_distance_and_directed_gate(self):
        with patch.object(wrapper.terminal, 'terminal_receipt', side_effect=mock_terminal) as call:
            result = wrapper.terminal_receipt(state=24, distance=F(2, 5), failure_bits=60)
        call.assert_called_once_with(q=1 << 128, distance=F(2, 5), first=None, failure_bits=60, precision=256)
        self.assertEqual(result['weight_range'], [7000, 8192])
        self.assertTrue(result['automatic_terminal_cutoff'])
        self.assertTrue(result['directed'])
        self.assertEqual(result['state_dimension'], 24)

    def test_threshold_rounding_can_move_auto_cutoff_one_weight(self):
        def rounded(**kwargs):
            result = mock_terminal(**kwargs)
            if kwargs['first'] is None:
                result.update(meets_requested_range_target=False, log2_upper_rational='-39')
            return result
        with patch.object(wrapper.terminal, 'terminal_receipt', side_effect=rounded) as call:
            result = wrapper.terminal_receipt()
        self.assertEqual(call.call_count, 2)
        self.assertEqual(result['weight_range'], [7001, 8192])
        self.assertTrue(result['meets_requested_range_target'])

    def test_source_pins_and_state_validation(self):
        for name in ('hd_field128_full_rs_state_terminal_certificate.py', 'hd_field128_full_rs_state_certificate.py',
                     'hd_field128_full_rs_state_screen.py', 'hd_field128_lift_tail_certificate.py'):
            self.assertIn(name, wrapper.source_pins())
        with self.assertRaises(ValueError):
            wrapper.terminal_receipt(state=48)
        with self.assertRaises(ValueError):
            wrapper.terminal_receipt(distance=.4)


if __name__ == '__main__':
    unittest.main()
