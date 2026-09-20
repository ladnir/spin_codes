import unittest
from fractions import Fraction as F
import math
import numpy as np
from flint import arb,ctx
from scipy.special import logsumexp
import endpoint_log_replay_diagnostic as matrix
import endpoint_zero_log_diagnostic as scalar
import endpoint_zero_exact_diagnostic as exact


class LogTests(unittest.TestCase):
    def test_entrywise_recurrence(self):
        values=np.arange(1,46,dtype=float).reshape(5,3,3)/50
        left=np.array([.4,.6]);right=np.array([.7,.3]);current=values
        while len(current)>1:
            current=np.maximum(left[0]*current[:-1]+right[0]*current[1:],left[1]*current[:-1]+right[1]*current[1:])
        np.testing.assert_allclose(np.exp(matrix.adaptive_matrix(np.log(values),left,right)),current[0],rtol=1e-13)

    def test_sparse_coefficients_at_long_length(self):
        ctx.prec=256;t=8;length=2048;kernel={0:1,8:1};z=arb(1)/8;a=.77;b=1.
        values=exact.zero_region(t,kernel,z,length)
        expected=length//t*float(logsumexp([t*math.log(a),t*math.log(b/8)]))
        actual=scalar.adaptive_log(values,np.array([a]),np.array([b]))
        self.assertAlmostEqual(actual,expected,places=8)


if __name__=='__main__':unittest.main()
