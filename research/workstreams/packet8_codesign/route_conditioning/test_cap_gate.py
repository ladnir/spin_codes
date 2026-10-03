from itertools import combinations
from math import comb
import unittest

import numpy as np

import cap_gate


class CapGateTests(unittest.TestCase):
    def test_cap_one_is_first_conditioning(self):
        local = np.arange(1, 17, dtype=float).reshape(4, 2, 2)/20
        for nu in (0., .5, 2.):
            np.testing.assert_allclose(cap_gate.marked(local, nu, cap=1, packet_bits=2),
                cap_gate.gate.potential_operators(local, nu, packet_bits=2))

    def test_cap_two_ordered_supports(self):
        local = np.array([[[.8, .2], [0., .7]], [[.4, .3], [.1, .6]],
                          [[.2, .5], [.3, .4]], [[.1, .4], [.4, .3]]])
        q, nu = 4, .7
        bare = cap_gate.marked(local, 0., packet_bits=2)
        marked = cap_gate.marked(local, nu, packet_bits=2)
        regional, _ = cap_gate.gate.prior.placement(marked, q, epochs=3, windows=3)
        direct = np.zeros((2, 2))
        for support in combinations(range(9), q):
            js = [sum(slot//3 == step for slot in support) for step in range(3)]
            matrix = np.eye(2)
            for j in js:
                matrix = matrix @ bare[j]
            direct += np.exp(nu*sum(min(j, 2) for j in js))*matrix/comb(9, q)
        np.testing.assert_allclose(np.exp(regional[q]), direct, atol=2e-13, rtol=2e-13)


if __name__ == '__main__':
    unittest.main()
