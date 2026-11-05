"""Tiny path/fold/selection tests; no bank, production or file writes."""
import unittest
from math import log2
from unittest.mock import patch
import numpy as np
import hd_field128_return_grid_paths_v38 as pilot


class Checks(unittest.TestCase):
    def test_paths(self):
        for a, b, d in ((.2, .7, .8), (.5, .3, .5), (2., 3., 4.)):
            matrix = np.array([[log2(a), log2(b)], [-np.inf, log2(d)]])
            for n in (1, 2, 7, 256):
                result = pilot.paths(matrix, n)
                quiet, birth = a**n, b*sum(a**k*d**(n-1-k) for k in range(n))
                self.assertAlmostEqual(result['log2_Rn00'], log2(quiet), places=10)
                self.assertAlmostEqual(result['log2_Rn01'], log2(birth), places=10)
                self.assertAlmostEqual(result['log2_Rn_row0'], log2(quiet+birth), places=10)
        with self.assertRaises(ValueError):
            pilot.paths(np.zeros((2, 2)))

    def test_winners_and_control(self):
        best = pilot.Best()
        matrices = [np.array([[2., 3.], [-np.inf, 4.]]), np.array([[1., 4.], [-np.inf, 5.]]),
                    np.array([[3., 2.], [-np.inf, 3.]])]
        for index, matrix in enumerate(matrices):
            best.add(index, matrix, [True, False, True, True])
        best.add(7, matrices[-1], [True]*4)  # Exact tie keeps first winner.
        self.assertEqual(best.winning_subset(), [1, 2])
        control = pilot.Best()
        for index in best.winning_subset():
            control.add(index, matrices[index], [True]*4)
        np.testing.assert_array_equal(control.logs, best.logs)
        result = pilot.snapshot(best, 19.)
        self.assertIsNone(result['winning_marker_indices'][1][0])
        self.assertIsNone(result['regional_log2'][1][0])
        self.assertAlmostEqual(result['total_score'], result['log2_Rn_row0']+19.)

    def test_nested_grid_fold(self):
        coefficients = np.array([[[.2], [.1]], [[.3], [.4]], [[0.], [0.]], [[.5], [.2]]])
        with np.errstate(divide='ignore'):
            logs = np.log2(coefficients)
        baseline = np.log2(np.array([[.5, 1.], [.4, .8]]))
        small = pilot.old.triangular_atlas(logs, (0., 0.), baseline, grid=4, blocks=8)
        large = pilot.old.triangular_atlas(logs, (0., 0.), baseline, grid=8, blocks=8)
        folded = np.zeros((4, 4, 4))
        for x in range(8):
            for y in range(8):
                folded[:, x % 4, y % 4] += np.exp2(large.logs[:, x, y])
        np.testing.assert_allclose(folded, np.exp2(small.logs), atol=2e-9, rtol=1e-9)
        self.assertTrue(np.all(folded[2] == 0))

    def test_platform_guard(self):
        with patch.object(pilot.sys, 'platform', 'win32'):
            with self.assertRaises(RuntimeError):
                with pilot.remote_limits():
                    pass


if __name__ == '__main__':
    unittest.main()
