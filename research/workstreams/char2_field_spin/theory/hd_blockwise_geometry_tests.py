"""Bounded exact regressions for the restarted direct-SPIN geometry gate."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_blockwise_geometry as gate
from high_distance_window_obstruction import lower_bound


class BlockwiseGeometryTests(unittest.TestCase):
    def test_no_padding_and_strict_packet_windows(self):
        for b, w, distance in product(gate.OUTERS, gate.WIDTHS, gate.DISTANCES):
            profile = gate.window_profile(b, w, distance)
            self.assertEqual(profile['k'], 1 << 20)
            self.assertEqual(profile['n'], 1 << 21)
            self.assertEqual(profile['rows'] % w, 0)
            self.assertEqual(profile['window'] % w, 0)
            self.assertLess(F(profile['window'], profile['rows']), distance)
            self.assertGreaterEqual(F(profile['window']+w, profile['rows']), distance)

    def test_original_gate_agreement(self):
        for b, w, distance in product(gate.OUTERS, gate.WIDTHS, gate.DISTANCES):
            profile = gate.window_profile(b, w, distance)
            for s in gate.STATES:
                mean, bound = lower_bound(b, b//2, profile['rows'], profile['window'], s)
                self.assertEqual(mean, profile['mean'])
                self.assertEqual(bound, gate.failure_lower(mean, profile['maximum'], b*s))

    def test_terminal_suffix_improves_old_gate(self):
        for b, w, s, distance in product(gate.OUTERS, gate.WIDTHS, gate.STATES, gate.DISTANCES):
            profile = gate.window_profile(b, w, distance)
            old = gate.failure_lower(profile['mean'], profile['maximum'], b*s)
            new = gate.failure_lower(profile['mean'], profile['maximum'], (b-1)*s)
            self.assertGreaterEqual(new, old)

    def test_exhaustive_small_failure_probabilities(self):
        for b, k, rows, window in ((4, 2, 4, 2), (4, 2, 5, 3), (6, 3, 3, 2)):
            law = gate.exhaustive_dimension_law(b, k, rows, window)
            mean = sum(size*mass for size, mass in law.items())
            expected = rows*gate.mean_row_dimension(b, k, F(window, rows))
            self.assertEqual(mean, expected)
            self.assertEqual(sum(law.values()), 1)
            for cost in range(k*window+2):
                actual = sum(mass for size, mass in law.items() if size > cost)
                self.assertLessEqual(gate.failure_lower(mean, k*window, cost), actual)

    def test_prefix_suffix_causal_witness_with_nonzero_final_state(self):
        # Four GF4 repetition-code rows: [m_i,m_i], packet width2, state1.
        # The first two source rows land in a first-region prefix and a
        # last-region suffix. Their shortening dimension2 exceeds one exit.
        # GF4 uses z^2=z+1; multiplication by z maps1->2,2->3,3->1.
        twice = (0, 2, 3, 1)
        message = (2, 1, 0, 0)
        routed = message + (0, 0, message[0], message[1])
        state, output, exits = 0, [], []
        for packet in range(4):
            x0, x1 = routed[2*packet:2*packet+2]
            output.extend((x0 ^ state, x1 ^ state))
            # Each local input-to-output block is identity and hence invertible.
            state ^= x0 ^ (twice[x1] if packet == 0 else x1)
            exits.append(state)
        self.assertEqual(exits[0], 0)
        self.assertNotEqual(exits[-1], 0)
        self.assertTrue(any(output))
        self.assertEqual(output[2:6], [0, 0, 0, 0])

    def test_exact_minimum_state_threshold(self):
        epsilon = F(1, 1 << 40)
        for b, w, distance, suffix in product(gate.OUTERS, gate.WIDTHS, gate.DISTANCES, (False, True)):
            profile = gate.window_profile(b, w, distance)
            minimum = gate.minimum_unexcluded_state(profile, terminal_suffix=suffix)
            cost = b-int(suffix)
            self.assertLessEqual(gate.failure_lower(profile['mean'], profile['maximum'], cost*minimum), epsilon)
            if minimum:
                self.assertGreater(gate.failure_lower(profile['mean'], profile['maximum'], cost*(minimum-1)), epsilon)

    def test_frozen_k20_thresholds_and_full_menu(self):
        # Each row lists widths16,32,64,128; entries are40%/45% minimum state.
        expected = {
            64: ((46, 245), (45, 241), (43, 241), (43, 241)),
            128: ((3, 45), (3, 45), (3, 45), (3, 38)),
            256: ((1, 6), (1, 6), (1, 5), (1, 3)),
            512: ((1, 1), (1, 1), (1, 1), (1, 1)),
        }
        checked = 0
        for b, w in product(gate.OUTERS, gate.WIDTHS):
            pair = tuple(gate.minimum_unexcluded_state(gate.window_profile(b, w, d)) for d in gate.DISTANCES)
            self.assertEqual(pair, expected[b][gate.WIDTHS.index(w)])
            for s, d in product(gate.STATES, gate.DISTANCES):
                receipt = gate.screen(b, w, s, d, balanced=False)
                self.assertEqual(receipt['excludes_requested_distance_and_margin'], s < pair[gate.DISTANCES.index(d)])
                self.assertTrue(receipt['zero_lower_bound_is_not_a_pass'])
                checked += 1
        self.assertEqual(checked, 160)

    def test_balanced_window_exact_poisson_binomial(self):
        for numerators in ((1, 2, 1, 2), (0, 1, 2), (2, 2, 2)):
            denominator = 3
            counts = [0]*(len(numerators)+1)
            for positions in product(range(denominator), repeat=len(numerators)):
                occupancy = sum(p < u for p, u in zip(positions, numerators))
                counts[occupancy] += 1
            self.assertEqual(gate.poisson_binomial_coefficients(numerators, denominator), counts)

    def test_balanced_budget_and_threshold(self):
        epsilon = F(1, 1 << 40)
        expected = {64: ((46, 247), (46, 247), (46, 247), (46, 246)),
                    128: ((3, 47),)*4, 256: ((1, 6),)*4, 512: ((1, 1),)*4}
        for b, w, distance in product(gate.OUTERS, gate.WIDTHS, gate.DISTANCES):
            profile = gate.balanced_profile(b, w, distance)
            self.assertLess(profile['total_window'], distance*profile['n'])
            self.assertGreaterEqual(profile['total_window']+w, distance*profile['n'])
            self.assertEqual(sum(profile['regional_packet_counts'])*w, profile['total_window'])
            self.assertEqual(profile['maximum'], F(profile['total_window'], 2))
            minimum = gate.minimum_unexcluded_state(profile)
            self.assertEqual(minimum, expected[b][gate.WIDTHS.index(w)][gate.DISTANCES.index(distance)])
            self.assertLessEqual(gate.failure_lower(profile['mean'], profile['maximum'], (b-1)*minimum), epsilon)
            self.assertGreater(gate.failure_lower(profile['mean'], profile['maximum'], (b-1)*(minimum-1)), epsilon)
            for state in gate.STATES:
                self.assertEqual(gate.screen(b, w, state, distance)['excludes_requested_distance_and_margin'],
                                 state < minimum)

    def test_validation(self):
        for args in ((63, 16, F(2, 5)), (64, 3, F(2, 5)), (64, 16, F(1, 2))):
            with self.assertRaises(ValueError):
                gate.window_profile(*args)
        with self.assertRaises(ValueError):
            gate.failure_lower(F(2), 1, 0)


if __name__ == '__main__':
    unittest.main()
