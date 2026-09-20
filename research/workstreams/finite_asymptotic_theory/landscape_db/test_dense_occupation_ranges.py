"""Exact region-mixture and integer-interval checks for dense coverage."""
import math
import unittest

import numpy as np

import activation_occupation as general
import dense_occupation_ranges as dense


class DenseOccupationTest(unittest.TestCase):
    def test_eligible_coefficient_bound_dominates_each_region_mixture(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.4,8)
        regions=general.region_logs(epoch,8,16,16)
        for q in (0,1,3,8,15,16):
            for p in (.2,.5,.8):
                exact=np.zeros((3,3))
                for j in range(q+1):
                    exact+=math.comb(q,j)*p**j*(1-p)**(q-j)*np.exp(regions[j])
                for eta in (-2.,0.,2.):
                    bound=np.exp(dense.region_mixture_upper(epoch,8,16,q,p,eta))
                    self.assertTrue(np.all(bound+1e-12>=exact))

    def test_endpoint_maximum_covers_every_integer_and_their_sum(self):
        for block in (2,8,32):
            for eta in (-3.,0.,3.):
                values=dense.point_logs(np.arange(1,65),block,64,11.,eta,-25.,12,.2)
                for lo,hi in ((1,64),(2,9),(17,48),(63,64)):
                    actual=float(np.logaddexp.reduce(values[lo-1:hi]))
                    bound=dense.interval_log(lo,hi,block,64,11.,eta,-25.,12,.2)
                    self.assertLessEqual(actual,bound+1e-10)
                self.assertTrue(np.all(np.diff(values,n=2)>=-1e-10))

    def test_serialized_moment_matches_region_power(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.3,8)
        ps=np.array([.2,.5,.8]);etas=np.array([-1.,0.,1.])
        actual=dense.serialized_moments(epoch,8,128,ps,etas)
        matrices=dense.epoch_mixture_logs(epoch,8,ps*np.exp(-np.logaddexp(0.,-etas)))
        expected=[general.terminal_log(m,16) for m in matrices]
        np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-13)


if __name__=='__main__':unittest.main()
