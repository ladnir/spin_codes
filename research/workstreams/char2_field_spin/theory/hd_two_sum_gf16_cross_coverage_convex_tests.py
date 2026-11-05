"""Small convex-marker controls; no regional/proof gate."""
from unittest.mock import patch
import unittest

import numpy as np

import hd_two_sum_gf16_cross_coverage_convex as gate


class ConvexCoverageTests(unittest.TestCase):
    def test_gapped_polynomial_defeats_adjacent_ratio(self):
        #16+X^4+2^-100 X^5+X^8: stationary marker is essentially logX=.5,
        #whereas the adjacent4/5 ratio suggests logX=100.
        word = np.full(9, -np.inf)
        word[[0, 4, 5, 8]] = [4., 0., -100., 0.]
        marker, diagnostics = gate.optimal_marker(word, 4)
        self.assertAlmostEqual(marker, .5, places=13)
        self.assertAlmostEqual(diagnostics['mean_degree'], 4., places=12)
        bad_marker = word[4]-word[5]
        bad_objective = float(np.logaddexp2.reduce(word+bad_marker*(np.arange(9)-4)))
        self.assertGreater(bad_objective-diagnostics['objective_log2'], 390.)
        self.assertEqual(diagnostics['bisections'], 48)

    def test_exact_gap_is_accepted(self):
        word = np.array([6., -np.inf, -np.inf, -np.inf, 0.])
        marker, diagnostics = gate.optimal_marker(word, 2)
        self.assertAlmostEqual(marker, 1.5, places=13)
        self.assertAlmostEqual(diagnostics['objective_log2'], 4., places=13)

    def test_endpoint_minimum_and_constant_objective(self):
        marker, diagnostics = gate.optimal_marker(np.array([0., 0., 3.]), 1)
        self.assertEqual(marker, 0.)
        self.assertTrue(diagnostics['endpoint_minimum'])
        marker, diagnostics = gate.optimal_marker(np.array([-np.inf, 7.]), 1)
        self.assertEqual(marker, 0.)
        self.assertEqual(diagnostics['objective_log2'], 7.)

    def test_absent_finite_bracket_and_invalid_inputs_fail(self):
        for word, floor in (([0., 0.], 1), ([-np.inf, -np.inf], 1),
                            ([0., np.nan], 1), ([0., np.inf], 1)):
            with self.assertRaises(ValueError):
                gate.optimal_marker(word, floor)
        with self.assertRaises(ValueError):
            gate.optimal_marker([100., -np.inf, 0.], 1, max_log_marker=1.)

    def test_timeout_and_resource_guards(self):
        with self.assertRaises(TimeoutError):
            gate.optimal_marker([0., 0., 0.], 1,
                check_time=lambda: (_ for _ in ()).throw(TimeoutError('toy')))
        with patch.object(gate, 'read_inputs', side_effect=AssertionError('must not read')):
            for invalid in (True, 0, 31, float('nan')):
                with self.assertRaises(ValueError):
                    gate.run(max_seconds=invalid)


if __name__ == '__main__':
    unittest.main()
