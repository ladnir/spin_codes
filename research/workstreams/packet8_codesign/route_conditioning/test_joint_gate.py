from fractions import Fraction
from itertools import combinations, product
from math import comb
import unittest

import numpy as np

import joint_gate


class JointTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, .2], [0., .7]], [[.4, .3], [.1, .6]], [[.2, .5], [.3, .4]]])

    def test_single_constraint_regressions(self):
        for nu in (0., .6, 2.):
            np.testing.assert_allclose(joint_gate.marked(self.local, nu, 0., packet_bits=2),
                joint_gate.gate.potential_operators(self.local, nu, packet_bits=2))
            np.testing.assert_allclose(joint_gate.marked(self.local, 0., nu, packet_bits=2),
                joint_gate.cap_gate.marked(self.local, nu, packet_bits=2))

    def test_joint_ordered_average_and_indicator(self):
        q, nu1, nu2, h1, h2 = 2, .4, .7, 4, 4
        bare = joint_gate.gate.potential_operators(self.local, 0., packet_bits=2)
        paths = []
        for support in combinations(range(6), q):
            occupancies = [sum(slot//2 == step for slot in support) for step in range(3)]
            value = np.eye(2)
            for j in occupancies:
                value = value @ bare[j]
            paths.append((sum(j > 0 for j in occupancies), sum(min(j, 2) for j in occupancies), value))
        regional, _ = joint_gate.gate.prior.placement(
            joint_gate.marked(self.local, nu1, nu2, packet_bits=2), q, epochs=3, windows=2)
        expected = sum(np.exp(nu1*a+nu2*b)*m for a, b, m in paths)/comb(6, q)
        np.testing.assert_allclose(np.exp(regional[q]), expected, rtol=3e-13, atol=3e-13)
        restricted, tilted = np.zeros((2, 2)), np.zeros((2, 2))
        for (a1, b1, m1), (a2, b2, m2) in product(paths, repeat=2):
            value = m1 @ m2/len(paths)**2
            if a1+a2 >= h1 and b1+b2 >= h2:
                restricted += value
            tilted += np.exp(nu1*(a1+a2-h1)+nu2*(b1+b2-h2))*value
        self.assertTrue(np.all(restricted <= tilted+1e-14))

    def test_exact_bad_route_sum(self):
        total, witness = joint_gate.exact_bad_union()
        one = joint_gate.route_counts.rational_bad_route_bound(119, 1483, Fraction(3, 20))
        two = joint_gate.cap_counts.rational_bad_route_bound(119, 2559, Fraction(9, 35), cap=2)
        self.assertEqual(total, one+two)
        self.assertLessEqual(total, Fraction(1, 1 << 59))
        self.assertTrue(witness['exact_59bit_sum_check'])
        self.assertFalse(witness['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
