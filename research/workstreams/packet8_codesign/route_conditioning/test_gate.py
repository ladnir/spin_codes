"""Tiny exact chronology checks for the route-conditioning comparison."""
from itertools import combinations, product
from math import comb
import unittest

import numpy as np

import gate


class ConditioningTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, .2], [0., .7]], [[.4, .3], [.1, .6]], [[.2, .5], [.3, .4]]])

    def test_thinning_equivalence_at_zero_route_tilt(self):
        original, _ = gate.prior.placement(self.local, 6, epochs=3, windows=2)
        mixed = gate.potential_operators(self.local, 0., packet_bits=2)
        potential, _ = gate.prior.placement(mixed, 6, epochs=3, windows=2)
        for q in range(7):
            expected = gate.prior.uniform_mixture(original, q, packet_bits=2)
            np.testing.assert_allclose(np.exp(potential[q]), np.exp(expected), rtol=2e-13, atol=2e-13)

    def test_marked_ordered_support_average(self):
        nu, q = .7, 3
        mixed = gate.potential_operators(self.local, nu, packet_bits=2)
        regional, _ = gate.prior.placement(mixed, q, epochs=3, windows=2)
        bare = gate.potential_operators(self.local, 0., packet_bits=2)
        direct = np.zeros((2, 2))
        for support in combinations(range(6), q):
            counts = [sum(slot//2 == step for slot in support) for step in range(3)]
            value = np.eye(2)
            for j in counts:
                value = value @ bare[j]
            direct += np.exp(nu*sum(j > 0 for j in counts))*value/comb(6, q)
        np.testing.assert_allclose(np.exp(regional[q]), direct, rtol=2e-13, atol=2e-13)

    def test_good_route_indicator_majorant(self):
        # Two regions, ordered noncommuting matrices, continuous state.
        q, threshold, nu = 2, 4, .8
        bare = gate.potential_operators(self.local, 0., packet_bits=2)
        paths = []
        for support in combinations(range(6), q):
            occupancies = [sum(slot//2 == step for slot in support) for step in range(3)]
            value = np.eye(2)
            for j in occupancies:
                value = value @ bare[j]
            paths.append((sum(j > 0 for j in occupancies), value))
        restricted, tilted = np.zeros((2, 2)), np.zeros((2, 2))
        for (h1, m1), (h2, m2) in product(paths, repeat=2):
            value = (m1 @ m2)/(len(paths)**2)
            if h1+h2 >= threshold:
                restricted += value
            tilted += np.exp(nu*(h1+h2-threshold))*value
        self.assertTrue(np.all(restricted <= tilted+1e-14))


if __name__ == '__main__':
    unittest.main()
