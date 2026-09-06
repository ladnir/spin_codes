"""Numerical checks against direct long finite sums; these are not proofs."""
import math
import unittest
import numpy as np
from scipy.special import gammaln, logsumexp
from convex_partition_screen import Partition, poisson_tail_series


class PartitionTests(unittest.TestCase):
    def test_tail_and_gradients(self):
        for start in [3, 20, 898]:
            for lam in [.01, start/3, start*.95, start*2]:
                k = np.arange(start, int(max(start+100, 4*lam+100)))
                logs = k*math.log(lam)-gammaln(k+1)
                total = logsumexp(logs)
                expected = np.exp(logs-total)@k
                got, mean = poisson_tail_series(start, math.log(lam))
                self.assertAlmostEqual(got, total, delta=1e-9)
                self.assertAlmostEqual(mean, expected, delta=1e-8)
        part = Partition(10)
        psi = .01*(np.arange(10)-4)**2
        slope, theta = .2, 1.
        value, mean, grad, ds = part.evaluate(psi, slope, theta)
        for i in [0, 5, 9]:
            direction = np.zeros(10)
            direction[i] = 1e-5
            numerical = (part.evaluate(psi+direction, slope, theta)[0]-part.evaluate(psi-direction, slope, theta)[0])/2e-5
            self.assertAlmostEqual(numerical, grad[i], delta=1e-8)
        numerical = (part.evaluate(psi, slope+1e-5, theta)[0]-part.evaluate(psi, slope-1e-5, theta)[0])/2e-5
        self.assertAlmostEqual(numerical, ds, delta=1e-8)


if __name__ == '__main__':
    unittest.main()
