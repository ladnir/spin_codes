import math
import unittest

import numpy as np

import joint_dense_witness_v1 as joint


class JointWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = joint.JointCoverRefiner({2: 2, 4: 1}, 4, 4, 2, {2: 2, 4: 1},
                                           [1, 0, 2, 0, 1], 8, [[2], [4]])

    def test_joint_gradient(self):
        corners = np.array([[2, 3, 3], [3, 3, 2]], dtype=float)
        objective, _, _, _ = self.model.objective(corners, 1.3, 0, np.array([0., .45, .81]))
        coordinates = np.array([-2.113, -.431, -.893, -.211, 1.317])
        _, gradient = objective(coordinates)
        for j in range(len(coordinates)):
            step = np.zeros(len(coordinates)); step[j] = 1e-6
            expected = (objective(coordinates+step)[0]-objective(coordinates-step)[0])/2e-6
            self.assertAlmostEqual(gradient[j], expected, places=7)

    def test_fixed_witness_does_not_regress(self):
        p = np.array([0., .5, 1.]); pi = np.array([.5, .3, .2])
        witness = dict(log_surprisal=-1., proposal=pi.tolist(), probabilities=p.tolist(),
                       log_density_costs=self.model.costs(p).tolist())
        lo, hi = [0, 0, 0], [5, 8, 8]
        value = self.model.direct(lo, hi, witness)
        box = dict(lower=lo, upper=hi, witness=witness, own_log_bound=value)
        result = self.model.refine_box(box, -20*math.log(2))
        actual = self.model.direct(lo, hi, result['witness'])
        self.assertAlmostEqual(actual, result['own_log_bound'], places=11)
        self.assertLessEqual(actual, value+1e-11)


if __name__ == '__main__':
    unittest.main()
