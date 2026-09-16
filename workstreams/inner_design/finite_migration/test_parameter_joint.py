"""Check generalized activation and fixed witnesses after joint tuning."""
from collections import Counter
import math
import unittest

import numpy as np
import parameter_joint_dense as joint
import parameter_hybrid_dense as hybrid
import parameter_hybrid_safe as safe


class JointTests(unittest.TestCase):
    def test_hybrid_preserves_its_selected_moment_family(self):
        with hybrid.maps.use():
            checker = safe.Dense(64,10,128,16,[-2.,-1.,0.,1.])
            self.assertEqual(checker.moment(1.,0.),math.inf)
            self.assertEqual(checker.mixed.moment(0.,0.),math.inf)
            for point in (np.array([0,0,1024]),np.array([600,423,1])):
                result = checker.evaluate(point,point)
                w = result['witness']
                proposal = np.array(w['proposal'])
                moment = checker.replay_moment(float(proposal@checker.ps),w)
                value = float(hybrid.base.typed.point_logs(point[None],checker.length,128,checker.log_gammas,
                    proposal,moment,checker.cutoff,math.exp(w['tilt']))[0])
                self.assertAlmostEqual(value,result['own_log_bound'],places=8)

    def test_activation_density_against_all_toy_syndromes(self):
        t,s = 6,3
        columns = [1,2,4,3,5,6]
        spectrum = Counter(sum((q&c).bit_count()%2 for c in columns) for q in range(1,1 << s))
        levels = sorted(spectrum)
        for theta,lam in ((.1,.3),(.5,1.),(.9,2.)):
            weights = np.zeros(1 << s)
            for x in range(1 << t):
                syndrome = 0
                for j,c in enumerate(columns):
                    if (x >> j) & 1:
                        syndrome ^= c
                w = x.bit_count()
                weights[syndrome] += theta**w*(1-theta)**(t-w)*math.exp(-lam*w)
            matrix = np.full((len(levels)+2,)*2,-np.inf)
            matrix[0,0] = math.log(weights[0])
            result = joint.activated.activate(matrix,spectrum,spectrum,t,theta,lam)
            self.assertEqual(result[0,0],matrix[0,0])
            self.assertEqual(result[0,1],-np.inf)
            for q in range(1,1 << s):
                w = sum((q&c).bit_count()%2 for c in columns)
                cap = math.exp(result[0,levels.index(w)+2])/spectrum[w]
                self.assertLessEqual(weights[q],cap+1e-14)

    def test_joint_witness_recomputes_on_every_vertex(self):
        with joint.maps.use():
            checker = joint.Dense(64,10,128,16,[-2.,-1.,0.,1.])
            lower,upper = np.array([598,420,0]),np.array([600,424,4])
            result = checker.evaluate(lower,upper)
            w = result['witness']
            p = np.array(w['proposal'])
            moment = checker.moment(float(p@checker.ps),w['tilt'])
            values = joint.base.typed.point_logs(result['corners'],checker.length,128,checker.log_gammas,
                p,moment,checker.cutoff,math.exp(w['tilt']))
            expected = float(max(values))+joint.base.typed.lattice_log_count(lower,upper)
            self.assertAlmostEqual(expected,result['own_log_bound'],places=8)


if __name__ == '__main__':
    unittest.main()
