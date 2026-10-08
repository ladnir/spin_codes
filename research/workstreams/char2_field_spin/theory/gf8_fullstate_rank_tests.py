"""Exact transfer checks, not a whole-code certificate."""
from fractions import Fraction as F
import math
import unittest

import numpy as np

import gf8_fullstate_rank as refined
from gf8_fullstate_screen import local_matrix as coarse
from gf8_rank_history_screen import cone_preserved


class FullStateRankTests(unittest.TestCase):
    def test_cone_and_coarse_domination(self):
        for state in (7, 8):
            for y, tau in ((F(3, 4), F(8)), (F(999, 1000), F(255)), (F(1), F(1))):
                for weight in range(9):
                    matrix = refined.local_matrix(weight, y, tau, state)
                    control = coarse(weight, y, tau, state)
                    self.assertTrue(cone_preserved(np.asarray(matrix, dtype=object)))
                    self.assertEqual(matrix[0][0], control[0][0])
                    self.assertEqual(sum(matrix[0][1:]), control[0][1])
                    for cap in range(1, state+1):
                        self.assertLessEqual(matrix[cap][0], control[1][0])
                        self.assertLessEqual(sum(matrix[cap][1:]), control[1][1])

    def test_aggregate(self):
        x, y, tau = .4, .7, 11.
        for state in (7, 8):
            direct = sum(math.comb(8, w)*x**w*np.asarray(refined.local_matrix(w, y, tau, state))
                         for w in range(9))
            np.testing.assert_allclose(refined.packet(math.log(x), math.log(y), math.log(tau), state),
                                       direct, rtol=2e-13, atol=1e-15)

    def test_live_return_bound_and_cone(self):
        y, tau = F(97, 100), F(6)
        for state in (7, 8):
            for weight in range(9):
                old = refined.local_matrix(weight, y, tau, state)
                sharp = refined.local_matrix(weight, y, tau, state, sharp_return=True)
                self.assertTrue(cone_preserved(np.asarray(sharp, dtype=object)))
                for row in range(state+1):
                    for col in range(state+1):
                        self.assertLessEqual(sharp[row][col], old[row][col])
        self.assertLess(refined.live_kernel_bound(tau, 8, 2, 9), tau**(-8)/10**8)

    def test_sharp_aggregate(self):
        x, y, tau, state = .4, .7, 11., 8
        direct = sum(math.comb(8, w)*x**w*np.asarray(refined.local_matrix(w, y, tau, state, True))
                     for w in range(9))
        np.testing.assert_allclose(refined.packet(math.log(x), math.log(y), math.log(tau), state, True),
                                   direct, rtol=2e-13, atol=1e-15)


if __name__ == '__main__':
    unittest.main()
