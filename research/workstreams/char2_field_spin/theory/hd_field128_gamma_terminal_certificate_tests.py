"""Toy causal inverse and mocked terminal adapter checks; no full-size gates."""

from fractions import Fraction as F
import unittest
from unittest.mock import patch

import hd_field128_gamma_terminal_certificate as wrapper
from hd_blockwise_maps import matvec, mul


def mock_terminal(**kwargs):
    first = 7000 if kwargs['first'] is None else kwargs['first']
    return dict(construction='prior construction', weight_range=[first, 8192],
                field_order=kwargs['q'], target_distance=str(kwargs['distance']),
                directed=True, whole_code_certificate=False,
                log2_upper_rational='-80', meets_requested_range_target=True,
                covers_all_supports_in_range=True, no_support_subset_union=True)


class GammaTerminalTests(unittest.TestCase):
    def test_causal_inverse_with_nontrivial_gamma(self):
        q, state, gamma = 16, 3, 7
        a = wrapper.local.frozen.lagrange(range(8), range(8, 11), q)
        c = tuple(zip(*a))
        h = wrapper.local.frozen.lagrange(range(state, 2*state), range(state), q)
        first, second = (2, 3, 5), (7, 9, 11)
        def advance(incoming, x):
            syndrome = matvec(c, x, q)
            v = tuple(mul(d, s^mul(gamma, r, q), q)
                      for d, s, r in zip(first, incoming, syndrome))
            return tuple(mul(d, r, q) for d, r in zip(second, matvec(h, v, q)))
        inputs = [tuple((5*i+3*t)%q for i in range(8)) for t in range(6)]
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

    def test_fresh_formula_call_and_metadata(self):
        with patch.object(wrapper.terminal, 'terminal_receipt', side_effect=mock_terminal) as call:
            result = wrapper.terminal_receipt(state=30, distance=F(2, 5), failure_bits=54, precision=384)
        call.assert_called_once_with(q=2**128, distance=F(2, 5), first=None,
                                     failure_bits=54, precision=384)
        self.assertEqual(result['construction'], wrapper.local.construction(30))
        self.assertEqual(result['feedback'], 'gamma A transpose')
        self.assertEqual(result['gamma_extension_degree'], 16)
        self.assertTrue(result['covers_all_weights_in_range_simultaneously'])
        self.assertTrue(result['no_union_over_tail_weights'])
        self.assertTrue(result['outer_independent_of_route_and_inner'])
        self.assertFalse(result['old_numerical_bounds_reused'])
        self.assertFalse(result['whole_code_certificate'])

    def test_automatic_cutoff_guard_and_fixed_cutoff_no_relabel(self):
        def rounded(**kwargs):
            result = mock_terminal(**kwargs)
            if kwargs['first'] is None:
                result['meets_requested_range_target'] = False
            return result
        with patch.object(wrapper.terminal, 'terminal_receipt', side_effect=rounded) as call:
            result = wrapper.terminal_receipt()
        self.assertEqual(call.call_count, 2)
        self.assertEqual(result['weight_range'], [7001, 8192])
        with patch.object(wrapper.terminal, 'terminal_receipt', return_value=dict(
                mock_terminal(first=7000, q=2**128, distance=F(2, 5)),
                meets_requested_range_target=False)) as call:
            result = wrapper.terminal_receipt(first=7000)
        self.assertEqual(call.call_count, 1)
        self.assertFalse(result['meets_requested_range_target'])
        self.assertFalse(result['automatic_terminal_cutoff'])

    def test_scope_and_source_guards(self):
        for state in (32, 48, True):
            with self.assertRaises(ValueError):
                wrapper.terminal_receipt(state=state)
        with self.assertRaises(ValueError):
            wrapper.terminal_receipt(distance=.4)
        pins = wrapper.source_pins()
        for name in ('hd_field128_gamma_terminal_certificate.py',
                     'hd_field128_full_rs_gamma_screen.py', 'hd_field128_lift_tail_certificate.py'):
            self.assertIn(name, pins)
        with patch.object(wrapper.terminal, 'terminal_receipt', side_effect=mock_terminal), \
             patch.object(wrapper, 'source_pins', side_effect=[pins, {}]):
            with self.assertRaises(RuntimeError):
                wrapper.terminal_receipt()


if __name__ == '__main__':
    unittest.main()
