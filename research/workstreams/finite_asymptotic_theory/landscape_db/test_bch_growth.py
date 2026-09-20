from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np

import activation_q1 as old
import activation_q1_refresh as refresh
import bch_growth_model as growth
from test_activation_q1 import actual_epochs, exact_coefficients, positive_regions


def exact_refresh(z):
    m0=(2*z*z+z**4)/3;m1=z/3+2*z**3/3
    zero=[[F(0) for _ in range(4)] for _ in range(4)]
    one=[[F(0) for _ in range(4)] for _ in range(4)]
    zero[0][0]=1;zero[1][2]=z*z;zero[2][2]=m0;zero[3][2]=min(z*z,F(3,2)*m0)
    one[0][1]=z
    for state,moment in ((1,z),(2,m1),(3,min(z,F(3,2)*m1))):
        one[state][0]=moment/3;one[state][3]=2*moment/3
    return zero,one


class RefreshTest(unittest.TestCase):
    def test_actual_shared_setup_is_dominated_and_old_bound_improves(self):
        for z in (F(1,3),F(1,2),F(3,4)):
            lam=np.array([-math.log(float(z))])
            for epochs in (1,2,3):
                actual=exact_coefficients(*positive_regions(*actual_epochs(z),epochs),4)
                exact=exact_coefficients(*positive_regions(*exact_refresh(z),epochs),4)
                self.assertTrue(all(a<=b for a,b in zip(actual,exact)))
                got=refresh.coefficient_logs(*refresh.epoch_logs(4,2,{2:2,4:1},lam),epochs,4)[0]
                prior=old.coefficient_logs(*old.region_logs(*old.epoch_logs(4,2,{2:2,4:1},lam),epochs),4)[0]
                np.testing.assert_allclose(got,[math.log(float(x)) for x in exact],atol=2e-13,rtol=2e-13)
                self.assertTrue(np.all(got<=prior+2e-13))

    def test_zero_epoch_refresh_has_no_repeated_density_penalty(self):
        z=.75;epochs=1000
        zero,one=refresh.epoch_logs(4,2,{2:2,4:1},np.array([-math.log(z)]))
        rz,_=refresh.region_logs(zero,one,epochs)
        self.assertAlmostEqual(rz[0,2,2],epochs*math.log((2*z*z+z**4)/3),places=10)
        self.assertAlmostEqual(rz[0,3,2],zero[0,3,2]+(epochs-1)*zero[0,2,2],places=10)


class GrowthTest(unittest.TestCase):
    def test_onset_probability_matches_support_enumeration(self):
        for b in (4,8,12):
            for w in range(1,b+1):
                threshold=b*.8
                actual=sum(min(1.,max(0.,min(support)+1-threshold))
                           for support in itertools.combinations(range(b),w))/math.comb(b,w)
                self.assertAlmostEqual(actual,growth.onset_probability(b,w),places=14)

    def test_moment_matches_independent_quadrature(self):
        from scipy.integrate import quad
        for b,w,a in ((8,2,.7),(64,12,4.4),(128,22,3.6)):
            exact=sum(math.comb(h-1,w-1)/math.comb(b,w)*
                      quad(lambda u:math.exp(-a*(h-u)/2),0,1,epsabs=1e-100)[0]
                      for h in range(w,b+1))
            self.assertAlmostEqual(math.log(exact),growth.onset_log_moment(b,w,a),places=11)

    def test_chernoff_dominates_onset_and_detects_distinct_regimes(self):
        for b,w in ((8,4),(32,8),(64,12),(128,22),(128,24),(128,26)):
            p=growth.onset_probability(b,w);bound=growth.onset_chernoff(b,w)['log_upper']
            if p:
                self.assertGreaterEqual(bound+1e-12,math.log(p))
            else:
                self.assertEqual(bound,-math.inf)


if __name__=='__main__':unittest.main()
