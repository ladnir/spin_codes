import math
import unittest

import numpy as np
from scipy.integrate import quad

import activation_q1_refresh as reference
import bch_growth_model as onset
import continuous_state_model as continuous


class ContinuousStateTest(unittest.TestCase):
    def test_zero_cancellation_recovers_onset_moment(self):
        tilts=np.array([.3,1.,4.])
        for b in (8,32,128):
            got=reference.coefficient_logs(*continuous.region_logs(0,.5,tilts),1,b)
            for i,a in enumerate(tilts):
                for w in (1,b//4,b):
                    self.assertAlmostEqual(got[i,w],onset.onset_log_moment(b,w,a),places=10)

    def test_transitions_match_location_integrals(self):
        p=1/31;mu=16/31;a=2.3
        zero,one=continuous.region_logs(p,mu,np.array([a]))
        self.assertAlmostEqual(math.exp(one[0,0,2]),quad(lambda v:math.exp(-a*mu*(1-v)),0,1)[0],places=13)
        self.assertAlmostEqual(math.exp(one[0,2,0]),quad(lambda v:p*math.exp(-a*mu*v),0,1)[0],places=13)
        self.assertAlmostEqual(math.exp(one[0,2,2]),(1-p)*math.exp(-a*mu),places=13)
        self.assertAlmostEqual(math.exp(zero[0,2,2]),math.exp(-a*mu),places=13)


if __name__=='__main__': unittest.main()
