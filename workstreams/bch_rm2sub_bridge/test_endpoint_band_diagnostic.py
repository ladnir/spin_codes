import math
import unittest
from fractions import Fraction as F
import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom
import endpoint_band_diagnostic as diag
import tightened_occupancy as tight


class EndpointTests(unittest.TestCase):
    def test_log_power(self):
        a=np.array([[.7,.3,0],[.1,.2,.7],[.1,0,.9]])
        with np.errstate(divide='ignore'):log=np.log(a)
        for n in (0,1,2,7,256):
            np.testing.assert_allclose(np.exp(diag.logpower(log,n)),np.linalg.matrix_power(a,n),rtol=1e-12,atol=1e-14)

    def test_full_occupancy_mixture_identity(self):
        t,s,length=8,4,32;spectrum={4:14,8:1};kernel={0:1,4:14,8:1};z=F(3,4);p=F(2,5)
        epoch=tight.epoch_matrices(t,s,spectrum,kernel,z,F,t)
        region=tight.regions(t,s,spectrum,kernel,z,F,length,length=length)
        def average(table,n):
            return np.array([float(sum((math.comb(n,j)*p**j*(1-p)**(n-j)*r[k]
                for j,r in enumerate(table)),F(0))) for k in range(9)]).reshape(3,3)
        np.testing.assert_allclose(average(region,length),np.linalg.matrix_power(average(epoch,t),length//t),rtol=1e-13)

    def test_zero_only_is_not_zero_terminal(self):
        matrix=np.array([[.5,.5,0],[.25,0,.75],[.1,0,.9]])
        with np.errstate(divide='ignore'):a=np.log(matrix)
        final=diag.logpower(a,10)[0]
        self.assertGreater(final[0],10*a[0,0])
        self.assertAlmostEqual(float(logsumexp(final)),0.,places=12)


if __name__=='__main__':unittest.main()
