"""Exact and bounded checks for the partial-row connected-window obstruction."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_blockwise_partial as proof
from hd_blockwise_geometry import balanced_profile, exhaustive_dimension_law, failure_lower


class PartialWindowTests(unittest.TestCase):
    def test_paired_windows_have_half_as_many_components(self):
        lengths = proof.WITNESS_LENGTHS
        regions, blocks = len(lengths), 256
        intervals = [(c*blocks+(blocks-length if c % 2 == 0 else 0),
                      c*blocks+(blocks if c % 2 == 0 else length))
                     for c, length in enumerate(lengths)]
        components = 1+sum(left[1] < right[0] for left, right in zip(intervals, intervals[1:]))
        self.assertEqual(components, 64)
        self.assertLess(intervals[-1][1], regions*blocks)
        self.assertEqual(sum(end-start for start, end in intervals)*64, 838848)

    def test_state4_exact_failure_probability(self):
        receipt = proof.paired_window_lower(4)
        bound = F(receipt['probability_lower'])
        self.assertGreater(bound, F(2115171265784, 10**16))
        self.assertGreater(bound, F(1, 1 << 40))
        self.assertTrue(receipt['rejected_at_40_bits'])
        self.assertEqual(receipt['constraints'], 256)

    def test_state8_and16_zero_mean_gate_not_certificate(self):
        for state in (8, 16):
            self.assertEqual(F(proof.paired_window_lower(state)['probability_lower']), 0)

    def test_outer256_comparison(self):
        for state in (4, 8):
            forty = proof.paired_window_lower(state, 256)
            fortyfive = proof.paired_window_lower(state, 256, distance=F(9, 20))
            self.assertEqual(F(forty['probability_lower']), 0)
            self.assertTrue(fortyfive['rejected_at_40_bits'])

    def test_exact_small_pair_shortening_lower(self):
        # Changing window starts does not change the uniform subset law.
        for b, rows, window in ((4, 4, 2), (4, 5, 3), (6, 3, 2)):
            law = exhaustive_dimension_law(b, b//2, rows, window)
            mean = sum(value*mass for value, mass in law.items())
            cost = b//2
            actual = sum(mass for value, mass in law.items() if value > cost)
            self.assertLessEqual(failure_lower(mean, b*window//2, cost), actual)

    def test_frozen_partial_row_witness(self):
        receipt = proof.validate_frozen_witness()
        self.assertEqual(receipt['dimension'], 720)
        self.assertEqual(receipt['state8_bad_subspace_dimension_lower'], 208)
        self.assertEqual(receipt['fully_contained_rows'], 0)
        self.assertLess(receipt['low_weight'], F(2, 5)*(1 << 21))
        self.assertFalse(receipt['setup_probability_claim'])

    def test_one_bounded_search_recounts_exactly(self):
        receipt = proof.route_witness(0, steps=4, jumps=(1, 2))
        self.assertEqual(receipt['steps'], 4)
        self.assertGreater(receipt['final_dimension'], receipt['initial_dimension'])
        self.assertEqual(receipt['selection_order_sha256'], proof.WITNESS_HASH)
        self.assertTrue(receipt['setup_failure_probability_not_certified'])


if __name__ == '__main__':
    unittest.main()
