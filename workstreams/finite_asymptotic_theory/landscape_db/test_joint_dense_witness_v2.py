import math
import unittest

import numpy as np

import joint_dense_witness_v2 as joint
import test_joint_dense_witness_v1 as previous


class SmoothedWitnessTests(previous.JointWitnessTests):
    @classmethod
    def setUpClass(cls):
        cls.model = joint.JointCoverRefiner({2: 2, 3: 1, 4: 1}, 4, 4, 2, {2: 2, 4: 1},
                                           [1, 0, 2, 0, 1], 8, [[2, 3], [4]])

    def test_gradient_at_shell_crossing(self):
        corners = np.array([[2, 3, 3], [3, 3, 2]], dtype=float)
        objective, _, _, _ = self.model.objective(corners, 1.3, 0, np.array([0., 3/7, .81]), .01)
        coordinates = np.array([-2.113, -.431, -.893, math.log(.75), 1.317])
        _, gradient = objective(coordinates)
        for j in range(len(coordinates)):
            step = np.zeros(len(coordinates)); step[j] = 1e-6
            expected = (objective(coordinates+step)[0]-objective(coordinates-step)[0])/2e-6
            self.assertAlmostEqual(gradient[j], expected, places=7)


if __name__ == '__main__':
    unittest.main()
