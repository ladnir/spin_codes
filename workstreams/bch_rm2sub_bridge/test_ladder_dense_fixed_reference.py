from fractions import Fraction as F
import unittest

from flint import arb,ctx

import ladder_dense_fixed_reference as fixed


class FixedTests(unittest.TestCase):
    def test_scalar_dominates_interior(self):
        ctx.prec = 256
        rows,lo,hi = 128,30,80
        nlo,nhi = F(1,4),F(7,8)
        r,eta,alpha = F(3,7),F(-50),arb(123)
        for beta in (F(0),F(71,4),F(128),F(256)):
            upper = fixed.scalar_secant(rows,lo,hi,nlo,nhi,r,eta,alpha,beta)
            for q in range(lo,hi+1):
                for j in range(11):
                    nu = nlo+(nhi-nlo)*F(j,10)
                    theta = fixed.aa(F(q,rows)*nu)
                    x = theta/(1-theta)/(fixed.aa(r)/(1-fixed.aa(r)))
                    actual = q*(alpha+fixed.aa(beta)*x.log()+fixed.aa(eta*nu))+256*rows*((1-theta)/(1-fixed.aa(r))).log()
                    self.assertGreaterEqual(upper,actual.lower())

    def test_fixed_point_matches_two_tilt_with_composition_cost(self):
        ctx.prec = 256
        checker = fixed.Checker(fixed.labels.identity.instance(20),fixed.labels.row_probabilities(.5))
        q,nu = 512,F(1,3)
        witness = dict(tilt=-80,input_tilt=0,eta=fixed.base.encode(F(10)))
        theta = F(q,checker.rows)*checker.density(nu,nu)[0]
        witness['fixed_r'] = fixed.base.encode(theta)
        actual = checker.bound(q,q,nu,nu,witness)
        expected = fixed.tilted.Checker.bound(checker,q,q,nu,nu,witness)
        import math
        expected += arb(math.comb(q+len(checker.bands),len(checker.bands))).log()
        self.assertLess(abs(float(actual-expected)),1e-7)


if __name__ == '__main__':
    unittest.main()
