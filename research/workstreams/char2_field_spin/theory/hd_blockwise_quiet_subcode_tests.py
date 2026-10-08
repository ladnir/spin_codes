"""Exact parameter arithmetic and actual finite-field quiet-state witnesses."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_blockwise_quiet_subcode as proof
from hd_blockwise_inner import dot


class QuietSubcodeTests(unittest.TestCase):
    def test_exact_k20_primary_caps(self):
        for b, s in product((128, 256), (4, 8, 16)):
            receipt = proof.quiet_cap(outer_length=b, outer_dimension=b//2, state=s)
            self.assertEqual(receipt['output_symbols'], 1 << 21)
            self.assertEqual(receipt['relative_distance_cap'], str(F(s, 32)))
            self.assertEqual(receipt['quiet_dimension_lower'], s)
            self.assertEqual(receipt['feedback_equations'], s*(32768-1))
            self.assertEqual(receipt['rejects_40_percent'], s < 13)
            self.assertTrue(receipt['deterministic_every_setup'])

    def test_exact_general_rounding(self):
        for b, k, rows, width, state in product((4, 6), (1, 2), (4, 8), (2, 4), range(5)):
            receipt = proof.quiet_cap(rows*k, b, k, width, state)
            equations = state*(b*rows//width-1)
            self.assertEqual(receipt['selected_rows'], equations//k+1)
            if receipt['guaranteed_nonzero']:
                self.assertGreater(receipt['selected_rows']*k, equations)
                self.assertLessEqual((receipt['selected_rows']-1)*k, equations)

    def test_strict_distance_target_and_state_thresholds(self):
        for b in (128, 256):
            for distance, expected in ((F(2, 5), 13), (F(9, 20), 15), (F(1, 2), 16)):
                self.assertEqual(proof.minimum_unexcluded_state(b, b//2, distance=distance), expected)
        self.assertEqual(proof.quiet_cap(state=16)['relative_distance_cap'], '1/2')
        self.assertFalse(proof.quiet_cap(state=16)['rejects_45_percent'])

    def test_terminal_constraint_is_separate(self):
        free = proof.quiet_cap(state=8)
        closed = proof.quiet_cap(state=8, terminal_free=False)
        self.assertEqual(closed['feedback_equations']-free['feedback_equations'], 8)
        self.assertEqual(closed['selected_rows'], free['selected_rows']+1)

    def test_smaller_packet_alternatives(self):
        for b, width, state, minimum in ((128, 32, 8, 7), (256, 16, 4, 4)):
            receipt = proof.quiet_cap(outer_length=b, outer_dimension=b//2,
                                      width=width, state=state)
            self.assertEqual(receipt['relative_distance_cap'], '1/2')
            self.assertEqual(receipt['quiet_dimension_lower'], state)
            self.assertEqual(proof.minimum_unexcluded_state(b, b//2, width), minimum)

    def test_alternative_paired_window_gates(self):
        for b, width, state, minimum, mean in (
                (128, 32, 8, 6, 344.6784748603134),
                (256, 16, 4, 1, 10.49452346472164)):
            gate = proof.alternative_window_gate(b, width, state)
            self.assertEqual(gate['components'], b//2)
            self.assertEqual(gate['exit_constraints'], 512)
            self.assertEqual(gate['total_weight'], 838848)
            self.assertEqual(gate['total_weight'] % width, 0)
            self.assertEqual(gate['failure_probability_lower'], '0')
            self.assertFalse(gate['rejects_40_percent_40_bits'])
            self.assertEqual(gate['minimum_paired_gate_unexcluded_state'], minimum)
            self.assertAlmostEqual(gate['expected_shortened_dimension'], mean, places=10)

    def test_one_packet_and_zero_state_edges(self):
        receipt = proof.quiet_cap(8, 4, 2, 16, 3)
        self.assertEqual(receipt['packets'], 1)
        self.assertEqual(receipt['feedback_equations'], 0)
        self.assertEqual(receipt['selected_rows'], 1)
        self.assertEqual(proof.quiet_cap(state=0)['selected_rows'], 1)

    def test_too_many_constraints_do_not_claim_a_witness(self):
        receipt = proof.quiet_cap(8, 4, 2, 2, 8)
        self.assertFalse(receipt['guaranteed_nonzero'])
        self.assertIsNone(receipt['relative_distance_cap'])

    def test_actual_GF4_witness(self):
        receipt = proof.finite_field_witness()
        self.assertEqual(receipt['variables'], 8)
        self.assertEqual(receipt['constraints'], 7)
        self.assertGreaterEqual(receipt['nullity'], 1)
        self.assertTrue(any(receipt['CA_values']))
        self.assertTrue(any(witness['terminal_state'] for witness in receipt['witnesses']))
        for witness in receipt['witnesses']:
            self.assertTrue(any(witness['message']))
            self.assertFalse(any(witness['emitted_states']))
            self.assertGreater(witness['weight'], 0)
            self.assertLessEqual(witness['weight'], 16)

    def test_small_RREF_basis(self):
        matrix = ((1, 2, 3), (2, 3, 1))
        basis = proof.nullspace(matrix, 3)
        self.assertTrue(basis)
        for vector in basis:
            for row in matrix:
                self.assertEqual(dot(row, vector, 4), 0)

    def test_validation(self):
        for options in (dict(message_symbols=3), dict(width=3), dict(state=-1), dict(outer_dimension=129)):
            with self.assertRaises(ValueError):
                proof.quiet_cap(**options)


if __name__ == '__main__':
    unittest.main()
