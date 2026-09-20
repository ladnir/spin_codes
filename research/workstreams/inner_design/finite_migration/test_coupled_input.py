"""Test coupled input normalizers and the logit curvature bound."""
from fractions import Fraction as F
import itertools
import math
import unittest

from flint import arb,ctx
import coupled_input_dense as coupled


class CoupledTests(unittest.TestCase):
    def test_subset_sums_with_fixed_input_tilt(self):
        ctx.prec = 256
        ps = [F(1,4),F(1,2),F(7,8),F(1)]
        gammas = [F(1,10),F(3,100),F(1,10000),F(1,1000000)]
        n,rows,x,z = 32,8,F(3,2),F(1,3)
        nlo,nhi,xi = F(3,10),F(9,10),-3
        for eta in (-2,2):
            weighted = [coupled.model.number(g)*(-coupled.model.number(eta*p)).exp()
                        for p,g in zip(ps,gammas)]
            low = sum(weighted[:2],arb(0))
            high = sum(weighted[2:],arb(0))
            for q in (1,2,3):
                lower,upper = F(q,rows)*nlo,F(q,rows)*nhi
                def value(theta):
                    return n*coupled.model.number(1-theta+theta*z/x).log()+coupled.model.number(eta*rows*theta)
                core = coupled.coupled_majorant([value(lower),value(upper)],lower,upper,n,eta*rows)
                upper_low = (core+q*low.log()).exp()
                upper_high = (core+q*(low+high*arb(-xi).exp()).log()+xi).exp()
                actual = [F(0),F(0)]
                for labels in itertools.product(range(4),repeat=q):
                    mean = sum((ps[j] for j in labels),F(0))/q
                    if not nlo <= mean <= nhi:
                        continue
                    theta = F(q,rows)*mean
                    term = math.prod(gammas[j] for j in labels)*(1-theta+theta*z/x)**n
                    actual[int(any(j >= 2 for j in labels))] += term
                self.assertTrue(coupled.model.number(actual[0]) < upper_low)
                self.assertTrue(coupled.model.number(actual[1]) < upper_high)

    def test_transformed_moment_identity_and_majorant(self):
        ctx.prec = 256
        n = 8
        values = [F(j*j+1,83) for j in range(n+1)]
        def moment(p):
            return sum((math.comb(n,j)*v*p**j*(1-p)**(n-j) for j,v in enumerate(values)),F(0))
        for x in (F(1,5),F(1),F(7)):
            for slope in (-100,F(-1,3),0,F(1,3),100):
                for lo,hi in ((F(1,20),F(1,10)),(F(1,4),F(3,4)),(F(9,10),F(19,20))):
                    def target(theta):
                        r = theta/(x+(1-x)*theta)
                        transformed = (1-theta+theta/x)**n*moment(r)
                        exact = sum((math.comb(n,j)*v/x**j*theta**j*(1-theta)**(n-j)
                                     for j,v in enumerate(values)),F(0))
                        self.assertEqual(transformed,exact)
                        return coupled.model.number(transformed).log()+coupled.model.number(slope*theta)
                    upper = coupled.coupled_majorant([target(lo),target(hi)],lo,hi,n,slope)
                    for j in range(1,20):
                        self.assertTrue(target(lo+(hi-lo)*F(j,20)) < upper)


if __name__ == '__main__':
    unittest.main()
