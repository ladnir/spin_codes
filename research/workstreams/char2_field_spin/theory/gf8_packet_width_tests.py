"""Exact transfer regressions and numerical rank tests for packet widths."""
from fractions import Fraction as F
from functools import lru_cache
import math
from random import Random
import unittest

import numpy as np

from gf8_packet_width_screen import (bound, gamma, joint_terms, local_matrix,
                                     minimum_return, packet, return_bound, validate)
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from cauchy_firstmoment_screen import geometry
from rate_design_coset import mul, power, xor_sum
from state_capacity_reference import rank
from state_rank_cap_transfer import return_minimum_weight
from symmetric_density_screen import local_matrix as old_local
from symmetric_density_screen import return_bound as old_return


@lru_cache(maxsize=None)
def maps(width, state, phase=0):
    inputs = tuple((phase*width)^u for u in range(width))
    nodes = tuple(128+v for v in range(state))
    def quotient(u, j):
        value = 1
        for k, v in enumerate(nodes):
            if k != j:
                value = mul(value, u^v)
        return value
    a = tuple(tuple(quotient(u, j) for j in range(state)) for u in inputs)
    return a, tuple(zip(*a))


class PacketWidthTests(unittest.TestCase):
    def test_exact_known_width16_state8_local_agreement(self):
        minima = tuple(return_minimum_weight(w) for w in range(1, 17))
        y, tau = F(3, 4), F(7)
        for w in range(17):
            expected = old_local(w, y, tau, 8, 16)
            if w:
                expected[1][0] = old_return(w, y, tau, 8)
            self.assertEqual(local_matrix(w, y, tau, 16, 8, minima), expected)

    def test_generic_return_bound_is_complete_minimum(self):
        for width, state in ((8, 4), (16, 4), (32, 8), (32, 16)):
            y, tau = F(2, 3), F(11, 2)
            for w in (1, 2, state, width):
                pure = (sum(math.comb(state, k)/tau**k for k in range(min(w-1, state-1)+1))
                        *y**max(1, width-state+1-w)/tau**state)
                actual = return_bound(w, y, tau, width, state)
                self.assertGreaterEqual(actual, 0)
                self.assertLessEqual(actual, pure)
                self.assertTrue(all(charge >= state and charge >= i+z
                                    for i, z, _, charge in joint_terms(w, width, state)))

    def test_batched_packet_matches_exact_local_sum(self):
        x, y, tau = F(4, 5), F(2, 3), F(17, 2)
        for width, state, minima in ((8, 4, ()), (8, 4, (5, 5, 4, 2, 1, 1, 1, 1)),
                                     (16, 8, ()), (32, 8, ()), (32, 16, ())):
            expected = sum((np.asarray(local_matrix(w, y, tau, width, state, minima), dtype=object)
                            *math.comb(width, w)*x**w for w in range(width+1)),
                           np.zeros((2, 2), dtype=object))
            actual = np.exp(packet(math.log(x), math.log(y), math.log(tau), width, state, minima))
            np.testing.assert_allclose(actual, np.asarray(expected, dtype=float), rtol=2e-12, atol=1e-300)

    def test_polynomial_maps_self_orthogonal_and_joint_rank(self):
        rng = Random(832168)
        checks = 0
        for width, state in ((8, 4), (16, 8), (32, 8), (32, 16)):
            a, c = maps(width, state)
            for i in range(state):
                for j in range(state):
                    self.assertEqual(xor_sum(mul(c[i][k], a[k][j]) for k in range(width)), 0)
            for _ in range(128):
                h = rng.randrange(1, state+1)
                w = rng.randrange(state-h+1, width+1)
                hs = sorted(rng.sample(range(state), h))
                js = sorted(rng.sample(range(width), w))
                inactive = sorted(set(range(width))-set(js))
                i = rng.randrange(w+1)
                z = rng.randrange(min(h-1, len(inactive))+1)
                selected = sorted(rng.sample(js, i)+rng.sample(inactive, z))
                numerical = [[int(v == k) for k in hs]+[c[v][j] for j in js]
                             for v in range(state)]
                numerical += [[mul(2, a[v][k]) for k in hs]+[int(v == j) for j in js]
                              for v in selected]
                actual = rank(numerical)
                expected = max(state, i+z, h+w-(width-i-z))
                self.assertGreaterEqual(actual, expected)
                self.assertGreaterEqual(expected, max(state, i+z, w-width+i+2*z+1))
                checks += 1
        self.assertEqual(checks, 512)

    def test_density_accounting_and_invalid_shapes(self):
        point = (-.2, -.7, math.log(12), -.1)
        g = geometry(20, 128, 84)
        weight = round(.1*g['n'])
        value, terms = bound(point, weight, dimension=84, components=True)
        self.assertAlmostEqual(value, sum(terms.values()), places=8)
        self.assertAlmostEqual(terms['route'], 128*modal_log_penalty(g['rows'])/math.log(2), places=9)
        self.assertAlmostEqual(terms['source'], g['rows']*source_row_log(weight/g['n'], point[2], point[3],
                                                                 128, 84, True)/math.log(2), places=8)
        for shape in ((7, 4), (8, 5), (32, 17), (8, 0)):
            with self.assertRaises(ValueError):
                validate(*shape)
        with self.assertRaises(ValueError):
            minimum_return(1, 8, 4, (3,))


if __name__ == '__main__':
    unittest.main()
