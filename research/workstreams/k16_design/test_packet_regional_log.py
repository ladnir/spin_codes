"""Small exact comparisons of the proposal-only log-semiring helper."""
from fractions import Fraction as Q
from math import log
import unittest

import numpy as np
from flint import arb, arb_mat, fmpq_mat, ctx
import packet_q1 as q1
import packet_regional_log as logarithmic
from packet_rs_length_proposal import scaled_placement, regional_uniform
from packet_rs_s20_proposal import log_power


class RegionalLogTests(unittest.TestCase):
    def test_noncommuting_log_products_match_direct_products(self):
        a = np.array([[.3, .2], [0, .4]])
        b = np.array([[.2, 0], [.5, .7]])
        logs = logarithmic.array_logs([a, b])
        ab = logarithmic.log_matmul(logs[0], logs[1])
        ba = logarithmic.log_matmul(logs[1], logs[0])
        np.testing.assert_allclose(np.exp(ab), a@b, rtol=1e-14, atol=1e-14)
        self.assertFalse(np.allclose(ab, ba))

    def test_ordered_placement_matches_exact_rationals(self):
        ops = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
               fmpq_mat([[2, 1], [0, 1]])]
        linear = np.array([[[float(m[i, j]) for j in range(2)] for i in range(2)] for m in ops])
        for epochs in (1, 2, 3, 7):
            degree = min(5, 2*epochs)
            exact = q1.placement(ops, epochs, 2, fmpq_mat, lambda x: x, maximum_groups=degree)
            expected = np.array([[[float(m[i, j]) for j in range(2)] for i in range(2)] for m in exact])
            actual = logarithmic.log_placement(logarithmic.array_logs(linear), degree,
                                               epochs=epochs, windows=2)
            np.testing.assert_allclose(np.exp(actual), expected, rtol=2e-13, atol=2e-13)

    def test_matches_scaled_method_in_its_stable_range(self):
        linear = np.array([[[.9, .01], [0, .7]], [[.2, .4], [.01, .5]],
                           [[.03, .02], [.003, .05]]])
        regional = logarithmic.log_placement(logarithmic.array_logs(linear), 7, epochs=5, windows=2)
        scaled = scaled_placement(linear, 7, epochs=5, windows=2)
        for q in (0, 3, 7):
            mixed = logarithmic.regional_uniform_log(regional, q)
            old, scale = regional_uniform(scaled, q)
            for regions in (1, 4, 64, 256):
                actual = logarithmic.log_power_matrix(mixed, regions)
                expected = log_power(old, regions)+regions*scale
                self.assertAlmostEqual(actual, expected, places=9)

    def test_tiny_arb_entries_are_logged_before_float_conversion(self):
        previous = ctx.prec
        try:
            ctx.prec = 192
            local = logarithmic.upper_logs([arb_mat([[arb('1e-1000'), 0], [0, 1]])])
            self.assertAlmostEqual(local[0, 0, 0], -1000*log(10), places=10)
            self.assertEqual(local[0, 0, 1], -np.inf)
            self.assertAlmostEqual(logarithmic.log_power_matrix(local[0], 64),
                                   -64000*log(10), places=8)
        finally:
            ctx.prec = previous

    def test_tiny_intermediate_edge_can_return_with_large_weight(self):
        # A linear scaling would lose M[0,1]=exp(-1000) relative to exp(1000).
        matrix = np.array([[0., -1000.], [1000., -np.inf]])
        squared = logarithmic.log_matmul(matrix, matrix)
        self.assertAlmostEqual(squared[0, 0], log(2), places=13)
        self.assertAlmostEqual(logarithmic.log_power_matrix(matrix, 2), log(2), places=13)

    def test_long_regional_paths_and_rare_packet_counts_survive(self):
        local = np.full((3, 1, 1), -10000.)
        regional = logarithmic.log_placement(local, 5, epochs=8, windows=2)
        np.testing.assert_allclose(regional[:, 0, 0], -80000., rtol=0, atol=1e-9)
        mixed = logarithmic.regional_uniform_log(regional, 5)
        self.assertAlmostEqual(logarithmic.log_power_matrix(mixed, 256), -80000.*256, places=6)
        q = 2000
        rare = -100*np.arange(q+1, dtype=float)[:, None, None]
        mixed = logarithmic.regional_uniform_log(rare, q)
        self.assertAlmostEqual(mixed[0, 0], -q*log(16), places=8)

    def test_continuous_state_is_not_reset_between_blocks(self):
        matrix = logarithmic.array_logs([[[.1, .9], [0, 2.]]])[0]
        actual = logarithmic.log_power_matrix(matrix, 8)
        expected = log(np.linalg.matrix_power(np.exp(matrix), 8)[0].sum())
        reset = 8*logarithmic.log_power_matrix(matrix, 1)
        self.assertAlmostEqual(actual, expected, places=12)
        self.assertGreater(actual, reset)

    def test_invalid_logs_and_structurally_zero_moment_reject(self):
        for bad in (np.array([[np.nan]]), np.array([[np.inf]]), np.zeros((2, 3))):
            with self.assertRaises(ValueError):
                logarithmic.log_power_matrix(bad, 2)
        with self.assertRaises(FloatingPointError):
            logarithmic.log_power_matrix(np.full((2, 2), -np.inf), 2)
        self.assertEqual(logarithmic.log_power_matrix(np.full((2, 2), -np.inf), 0), 0)
        with self.assertRaises(FloatingPointError):
            logarithmic.log_power_matrix(np.array([[1e308]]), 4)
        with self.assertRaises(ValueError):
            logarithmic.upper_logs([arb_mat([[-1]])])
        with self.assertRaises(ValueError):
            logarithmic.log_placement(np.zeros((3, 1, 1)), 5, epochs=2, windows=2)
        with self.assertRaises(ValueError):
            logarithmic.regional_uniform_log(np.zeros((3, 1, 1)), True)


if __name__ == '__main__':
    unittest.main()
