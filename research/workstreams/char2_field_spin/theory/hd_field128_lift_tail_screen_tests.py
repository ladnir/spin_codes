"""Small reference and coverage checks; no extension-field numerical gates."""

from fractions import Fraction as F
from math import log2
import unittest

from flint import arb, ctx
import numpy as np

import hd_field128_lift_tail_screen as tail
from hd_field128_lift_screen import Budget
from hd_randomized_outer_proof import Parameters, ball_log2
from hd_shifted_run_potential import run_log_screen
from hd_shifted_swap_potential import marginal_log_screen
from hd_shifted_swap_screen import positive_pair_matrix, regional, SwapPositiveWitness


class TailScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z, cls.alpha = .9, .125
        cls.odds = ((-4., -3.), (-3., -2.))
        cls.packet = tail.PositivePacket(cls.z, cls.alpha, q=256)
        cls.screen = tail.PositiveBankScreen(cls.z, cls.alpha, cls.odds, -8., q=256,
                                            outer_measure='inherited-diagonal', _packet=cls.packet)

    def test_positive_pair_q256_replay(self):
        for odds in self.odds:
            expected = np.log2(positive_pair_matrix(self.z, self.alpha, odds))
            np.testing.assert_allclose(self.packet.pair_logs(odds), expected, rtol=0, atol=2e-12)

    def test_regional_normalize_before_bank_min(self):
        weight = 3
        actual = self.screen.regional(weight)
        for index, full in enumerate(tail.full_counts(4096, weight)):
            histogram = (4096-weight+full, weight-2*full, full)
            expected = np.minimum.reduce([regional(self.z, self.alpha, histogram, odds, method='global')[0]
                                           for odds in self.odds])
            np.testing.assert_allclose(actual[index], expected, rtol=0, atol=2e-10)

    def test_small_support_sums_match_q_agnostic_helpers(self):
        for weight in (1, 2, 3, 6, 8):
            full = tail.full_counts(4, weight)
            regions = np.array([np.log2([[.7+f/100, .13], [.001, .2+f/200]]) for f in full])
            mapping = dict(zip(map(int, full), regions))
            for ratio in (-7., 0., 3.):
                for method, old in (('marginal', marginal_log_screen), ('runs', run_log_screen)):
                    actual = tail.support_sum_logs(4, weight, regions, ratio, method=method, repetitions=4)
                    expected = old(4, weight, mapping, ratio, repetitions=4)
                    self.assertAlmostEqual(actual, expected, places=11)

    def test_box_q256_replay_covers_clipped_vertices(self):
        ranges = ((1, 3), (0, 1), (0, 1))
        old = SwapPositiveWitness.create(self.z, self.alpha, self.odds).evaluate(ranges)[0]
        actual = self.screen.box(ranges, marker_indices=(0, 1))
        self.assertAlmostEqual(actual['log2_upper_screen'], old, places=7)
        self.assertTrue(actual['covers_entire_clipped_box'])
        self.assertFalse(actual['directed'])

    def test_band_checks_every_integer_weight_and_records_unfinished_work(self):
        result = self.screen.band(1, 3, max_seconds=10)
        self.assertEqual(result['evaluated_weights'], [1, 2, 3])
        self.assertTrue(result['every_requested_weight_evaluated'])
        self.assertIsNone(result['pending_weight_range'])
        self.assertFalse(result['whole_code_certificate'])
        incomplete = self.screen.band(1, 3, max_seconds=1e-12)
        self.assertEqual(incomplete['evaluated_weights'], [])
        self.assertFalse(incomplete['every_requested_weight_evaluated'])
        self.assertFalse(incomplete['requested_band_passes_float_threshold'])
        self.assertEqual(incomplete['pending_weight_range'], [1, 3])

    def test_terminal_q256_exact_mark_and_point_cap(self):
        result = tail.terminal_tail_screen(q=256, distance=F(7, 20), first=7649)
        with ctx.workprec(192):
            parameters = Parameters(field=256, state=16, limit=734003)
            ball = ball_log2(parameters, mark=F(result['mark']))
            expected = ball-7649*(96*8+32*arb(255).log()/arb(2).log())-arb(255).log()/arb(2).log()
            self.assertAlmostEqual(result['log2_upper_screen'], float(expected), places=10)
        self.assertEqual(result['weight_range'], [7649, 8192])
        self.assertTrue(result['no_support_subset_union'])
        self.assertFalse(result['directed'])

    def test_log_domain_extension_packet_handles_large_markers(self):
        packet = tail.PositivePacket(.5, F(1, 128), q=1 << 128)
        value = packet.pair_logs((100., 200.))
        self.assertTrue(np.all(np.isfinite(value)))
        self.assertGreater(value.max(), 1024)

    def test_scope_validation(self):
        for weight in (0, 8193):
            with self.assertRaises(ValueError):
                self.screen.regional(weight)
        with self.assertRaises(ValueError):
            self.screen.band(1, 4, max_weights=3)
        with self.assertRaises(ValueError):
            self.screen.box(((2, 2), (1, 1), (1, 1)))
        with self.assertRaises(TimeoutError):
            self.packet.pair_logs((0., 0.), _budget=Budget(1e-12))


if __name__ == '__main__':
    unittest.main()
