from fractions import Fraction as F
import math
import unittest

import numpy as np

from symmetric_density_screen import bound, packet, return_bound
from superregular_return_screen import local_matrix
from state8_joint_return import joint_return_bound
from state8_density_screen import packet as state8_packet


class SymmetricDensityTests(unittest.TestCase):
    def test_packet_matches_direct_positive_sum(self):
        x, y, tau = .4, .7, 5.
        expected = sum(math.comb(16,w)*x**w*np.asarray(local_matrix(w,y,tau,state=6))
                       for w in range(17))
        actual = np.exp(packet(math.log(x), math.log(y), math.log(tau), 6))
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_threshold_monotonicity(self):
        point = (-2., -.2, math.log(16), 0.)
        self.assertGreater(bound(point, 10000, k_log=16, distance=F(1,5)),
                           bound(point, 10000, k_log=16, distance=F(1,10)))

    def test_restricted_state_domain(self):
        for state in (0, 9):
            with self.assertRaises(ValueError):
                packet(-2., -.2, math.log(16), state)

    def test_exact_eight_state_return_reproduction(self):
        for y,tau in ((F(1,3),F(3)),(F(7,8),F(19)),(F(1),F(255))):
            for w in range(1,17):
                self.assertEqual(return_bound(w,y,tau,8),joint_return_bound(w,y,tau))

    def test_joint_never_exceeds_pure(self):
        for state in (2,4,6):
            for w in range(1,17):
                pure = local_matrix(w,F(2,5),F(7),state=state)[1][0]
                self.assertLessEqual(return_bound(w,F(2,5),F(7),state),pure)

    def test_joint_packet_reproduces_eight_state(self):
        point = (-.6,-.3,math.log(11))
        np.testing.assert_allclose(packet(*point,8,True),state8_packet(*point,joint=True),rtol=1e-13)


if __name__ == '__main__':
    unittest.main()
