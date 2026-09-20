"""Check whole-moment selection and the extended reference bank."""
import math
import unittest

import numpy as np
import parameter_mixed_extended as extended


class MixedTests(unittest.TestCase):
    def test_whole_moment_minimum(self):
        base = extended.prior.base
        checker = extended.prior.Dense(64,10,128,16,[])
        _,_,_,caps,low = base.prepare(64,10)
        for theta,tilt in ((.01,-4.),(.4,0.),(.999,1.)):
            epoch = base.g.epochs(checker.record['spectrum'],checker.kernel,
                checker.record['feedback_columns'],caps,low,math.exp(tilt),maximum=64)
            matrix = np.full(epoch.shape[1:],-np.inf)
            for j in range(65):
                term = math.log(math.comb(64,j))+j*math.log(theta)+(64-j)*math.log1p(-theta)
                np.logaddexp(matrix,epoch[j]+term,out=matrix)
            occupation = base.g.terminal(matrix,checker.block*checker.length//64)
            fourier = base.Dense.moment(checker,theta,tilt)
            self.assertAlmostEqual(checker.moment(theta,tilt),min(occupation,fourier),places=10)

    def test_uniform_bank_dominates_every_shell(self):
        checker = extended.Dense(64,10,128,16,[])
        counts,_ = extended.prior.base.grid.outer(128)
        scale,ps,costs = checker.probability_banks[-1]
        self.assertEqual(scale,0.)
        for band,p,cost in zip(checker.bands,ps[1:],costs[1:]):
            for w in band:
                if p == 1.:
                    self.assertEqual(w,128)
                    value = math.log(counts[w])
                else:
                    value = math.log(counts[w])-math.log(math.comb(128,w))-w*math.log(p)-(128-w)*math.log1p(-p)
                self.assertLessEqual(value,cost+1e-12)


if __name__ == '__main__':
    unittest.main()
