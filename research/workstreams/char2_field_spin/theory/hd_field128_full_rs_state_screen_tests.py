"""Exact fixed-map/local tests for variable state; no distance gate."""

from fractions import Fraction as F
from itertools import combinations
from math import comb
import unittest

import numpy as np

import hd_field128_full_rs_state_screen as variable
import hd_field128_full_rs_screen as fixed
from hd_field128_full_rs_screen_tests import toy_moments
from hd_blockwise_maps import maps, matmul, mul, power


def rank(rows, q=256):
    matrix = [list(row) for row in rows]
    row = 0
    for column in range(len(matrix[0]) if matrix else 0):
        pivot = next((i for i in range(row, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inverse = power(matrix[row][column], q-2, q)
        matrix[row] = [mul(value, inverse, q) for value in matrix[row]]
        for i in range(row+1, len(matrix)):
            factor = matrix[i][column]
            matrix[i] = [a^mul(factor, b, q) for a, b in zip(matrix[i], matrix[row])]
        row += 1
        if row == len(matrix):
            break
    return row


class VariableStateTests(unittest.TestCase):
    def test_actual_subfield_maps_ca_zero_and_nonsingular_refresh(self):
        for state in variable.STATES:
            a, c, h = variable.fixed_maps(state)
            self.assertEqual(len(a), 64)
            self.assertEqual(len(a[0]), state)
            self.assertFalse(any(x for row in matmul(c, a, 256) for x in row))
            self.assertEqual(rank(h), state)
            self.assertEqual(rank(a), state)
            for size in (1, 2, 3, state):
                selected = tuple(tuple(h[i][j] for j in range(state-size, state)) for i in range(size))
                self.assertEqual(rank(selected), size)
        self.assertEqual(variable.fixed_maps(16)[:2], maps(width=64, state=16, base=128, q=256))

    def test_non_power_of_two_state_toy_and_all_refresh_minors(self):
        for state in (2, 3, 4):
            a = variable.lagrange(range(8), range(8, 8+state), 16)
            self.assertFalse(any(x for row in matmul(tuple(zip(*a)), a, 16) for x in row))
            for positions in ((), (0,), (0, 1), (1, 6)):
                quiet, actual_return = toy_moments(a, positions, q=16)
                k, birth, cap, live = variable.exact_moments(len(positions), F(1, 2), 16, 8, state)
                self.assertEqual(quiet, k)
                self.assertGreaterEqual(cap, actual_return)
                self.assertGreaterEqual(birth, actual_return)
                self.assertGreaterEqual((16**state-1)*live, actual_return)
            h = variable.lagrange(range(state, 2*state), range(state), 16)
            for size in range(1, state+1):
                for rows in combinations(range(state), size):
                    for columns in combinations(range(state), size):
                        self.assertEqual(rank(tuple(tuple(h[i][j] for j in columns) for i in rows), 16), size)

    def test_positive_formulas_and_logs_match_exact_moments(self):
        for state in (2, 3, 4):
            z, q, n = F(7, 8), 16, 8
            logs = variable.local_logs(float(z), q, n, state).reshape((n+1, 4))
            t = z+(1-z)/q
            for j in range(n+1):
                k, birth, ret, live = variable.exact_moments(j, z, q, n, state)
                self.assertEqual(birth, t**j-k)
                m = n-j
                flat = sum((F(comb(m, e), q**min(state, e))*z**(m-e)*(1-z)**e for e in range(m+1)), F(0))
                self.assertEqual(live, t**j*(q**state*flat-1)/(q**state-1))
                exact = variable.exact_packet(j, z, q, n, state)
                np.testing.assert_allclose(logs[j], [variable.old.log_fraction(x) for x in exact], rtol=0, atol=1e-11)

    def test_state16_replays_frozen_local_and_regional_prefix(self):
        for joint in (False, True):
            np.testing.assert_array_equal(variable.local_logs(.5, 256, state=16, joint=joint),
                                          fixed.local_logs(.5, 256, joint=joint))
            for j in (0, 1, 16, 33, 64):
                self.assertEqual(variable.exact_packet(j, F(1, 2), 256, joint=joint),
                                 fixed.exact_packet(j, F(1, 2), 256, joint=joint))
        a = variable.regional_screen(.5, .1, 3, q=256, state=16)
        b = fixed.regional_screen(.5, .1, 3, q=256)
        self.assertEqual(a, b)

    def test_large_state_field128_tiny_marks_keep_positive_logs(self):
        z, q = F(1, 1 << 200), 1 << 128
        for state in (24, 32):
            logs = variable.local_logs(float(z), q, state=state)[1].ravel()
            exact = variable.exact_packet(1, z, q, state=state)
            np.testing.assert_allclose(logs, [variable.old.log_fraction(x) for x in exact], rtol=0, atol=1e-9)
            self.assertTrue(np.all(np.isfinite(logs)))
            self.assertLess(logs[2], -state*128)

    def test_alpha1_pairs_equal_single_region_all_public_states(self):
        for state in variable.STATES:
            a = variable.regional_screen(.5, 1., 3, q=256, state=state, paired=True)
            b = variable.regional_screen(.5, 1., 3, q=256, state=state, paired=False)
            self.assertEqual(a, b)
        with self.assertRaises(ValueError):
            variable.regional_screen(.5, .1, 2, state=20)
        with self.assertRaises(ValueError):
            variable.local_logs(.5, 256, state=33)
        self.assertEqual(len({variable.construction(s) for s in variable.STATES}), 3)
        pins = variable.source_pins()
        self.assertIn('hd_field128_full_rs_state_screen.py', pins)
        self.assertIn('hd_field128_full_rs_screen.py', pins)
        self.assertIn('hd_blockwise_maps.py', pins)


if __name__ == '__main__':
    unittest.main()
