from itertools import combinations
from math import comb
import unittest

import numpy as np

import capped_gate as cg


class CappedGateTests(unittest.TestCase):
    def test_joint_marker_matches_ordered_supports(self):
        local = np.array([[[.8,.2],[0.,.7]], [[.4,.3],[.1,.6]],
                          [[.2,.5],[.3,.4]], [[.1,.4],[.4,.3]]])
        caps, nus, q = [2,3], [.4,.7], 4
        bare = cg.marked(local,caps,[0.,0.],packet_bits=2)
        marked = cg.marked(local,caps,nus,packet_bits=2)
        regional,_ = cg.gate.prior.placement(marked,q,epochs=3,windows=3)
        direct = np.zeros((2,2))
        for support in combinations(range(9),q):
            js = [sum(slot//3==step for slot in support) for step in range(3)]
            matrix = np.eye(2)
            for j in js:
                matrix = matrix@bare[j]
            penalty = sum(nu*sum(min(j,cap) for j in js) for nu,cap in zip(nus,caps))
            direct += np.exp(penalty)*matrix/comb(9,q)
        np.testing.assert_allclose(np.exp(regional[q]),direct,atol=3e-13,rtol=3e-13)

    def test_disabled_second_marker(self):
        local = np.arange(1,37,dtype=float).reshape(9,2,2)/40
        np.testing.assert_array_equal(cg.marked(local,[3],[2.]),
                                     cg.marked(local,[2,3],[0.,2.]))

    def test_invalid_marker(self):
        local = np.ones((9,2,2))
        for caps,nus in (([2,2],[1.,2.]),([3],[-1.]),([3],[float('nan')]),([9],[1.])):
            with self.assertRaises(ValueError):
                cg.marked(local,caps,nus)


if __name__ == '__main__':
    unittest.main()
