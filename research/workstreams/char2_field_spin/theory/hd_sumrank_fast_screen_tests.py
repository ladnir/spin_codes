from fractions import Fraction as F
import unittest

import numpy as np

from hd_sumrank_fast_screen import moments,profile_data,fractional,packet_logs,bound
from hd_sumrank_transfer import profiles,fractional_packets,screen


class FastScreenTests(unittest.TestCase):
    def test_exact_moment_reference(self):
        for state,tiles,z in ((4,1,F(1,4)),(4,2,F(7,8)),(8,2,F(1,4)),(4,8,F(1,2))):
            k,t,v=moments(float(z),state,tiles)
            for index,(profile,j,multiplicity,ek,et,ev) in enumerate(profiles(z,state=state,tiles=tiles)):
                for actual,expected in ((k[index],ek),(t[index],et),(v[index],ev)):
                    self.assertAlmostEqual(actual/float(expected),1,places=12)

    def test_uniform_pre_transfer_reference(self):
        z,alpha,state,tiles=F(1,4),F(1,8),4,2
        actual=fractional(packet_logs(float(z),state,'uniform_pre',tiles),float(alpha),state,tiles)
        expected=np.array(fractional_packets(z,alpha,mixer='uniform_pre',state=state,tiles=tiles))
        np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)
        selected=screen(512,z,alpha,mixer='uniform_pre',state=state,tiles=tiles)
        value=bound(512,float(z),float(alpha),selected['marker_log_odds'],actual,tiles)
        self.assertAlmostEqual(value,selected['log2_upper_screen'],places=7)

    def test_marked_moment_reference(self):
        z,v,state,tiles=F(5,6),F(2),8,2
        k,t,flat=moments(float(z),state,tiles,float(v))
        for index,(profile,j,multiplicity,ek,et,ev) in enumerate(profiles(z,state=state,tiles=tiles,input_mark=v)):
            for actual,expected in ((k[index],ek),(t[index],et),(flat[index],ev)):
                self.assertAlmostEqual(actual/float(expected),1,places=8)


if __name__=='__main__':
    unittest.main()
