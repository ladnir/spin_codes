"""Exact local checks for width8/state7 and state8 candidate screens."""
from fractions import Fraction as F
from itertools import combinations
import math
import unittest

import numpy as np

import gf8_fullstate_screen as target
from gf8_packet_maps import determinant


class FullStateTests(unittest.TestCase):
    def test_invertibility_and_exact_return_minima(self):
        for state in (7, 8, 12, 16):
            result = target.validate_maps(state)
            self.assertEqual(result['return_rank'], 8)
            self.assertEqual(result['state_schur_rank'], state)
            self.assertFalse(result['ca_is_zero'])

    def test_all_cauchy_minors(self):
        for state in (7, 8):
            a, c = target.maps(state)
            count = 0
            for size in range(1, state+1):
                for rows in combinations(range(8), size):
                    for columns in combinations(range(state), size):
                        self.assertNotEqual(determinant([[a[i][j] for j in columns] for i in rows]), 0)
                        self.assertNotEqual(determinant([[c[j][i] for i in rows] for j in columns]), 0)
                        count += 2
            self.assertEqual(count, 2*(math.comb(8+state, state)-1))

    def test_all_zero_outputs_not_accidentally_excluded(self):
        y, tau = F(1, 100), F(11)
        for state in (7, 8):
            for weight in range(1, 9):
                z = 8-weight
                if z >= state:
                    continue
                all_zero_term = (1-y)**8*target.gamma(tau, state)[z]/tau**weight
                self.assertGreater(all_zero_term, 0)
                self.assertGreaterEqual(target.local_matrix(weight, y, tau, state)[1][1], all_zero_term)
        self.assertEqual(target.local_matrix(1, y, tau, 8)[0][0], 0)
        for state in (12, 16):
            self.assertEqual(len(target.gamma(tau, state)), 9)
            all_zero_quiet = (1-y)**8*target.gamma(tau, state)[8]
            self.assertGreater(all_zero_quiet, 0)
            self.assertGreaterEqual(target.local_matrix(0, y, tau, state)[1][1], all_zero_quiet)

    def test_aggregate_exact_sum(self):
        x, y, tau = F(3, 5), F(7, 9), F(17, 2)
        for state in (7, 8, 12, 16):
            expected = [[sum(math.comb(8, w)*x**w*target.local_matrix(w, y, tau, state)[i][j]
                             for w in range(9)) for j in range(2)] for i in range(2)]
            np.testing.assert_allclose(target.packet(math.log(x), math.log(y), math.log(tau), state),
                                       np.asarray(expected, dtype=float), rtol=2e-13, atol=1e-15)


if __name__ == '__main__':
    unittest.main()
