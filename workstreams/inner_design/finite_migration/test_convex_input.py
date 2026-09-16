"""Check the positive-polynomial logit interpolation majorant."""
from fractions import Fraction as F
import unittest

from flint import arb,ctx
import convex_input_dense as convex


class InterpolationTests(unittest.TestCase):
    def test_positive_polynomials(self):
        ctx.prec = 256
        for degree in (4,8,16):
            coefficients = [F((j*13+7)%19+1,23) for j in range(degree+1)]
            def moment(p):
                return sum((c*p**j*(1-p)**(degree-j) for j,c in enumerate(coefficients)),F(0))
            for a,b in ((F(1,100),F(1,20)),(F(1,4),F(3,4)),(F(9,10),F(99,100))):
                bound = convex.interpolate(convex.model.number(moment(a)).log(),
                    convex.model.number(moment(b)).log(),a,b,degree)
                for j in range(1,20):
                    p = a+(b-a)*F(j,20)
                    self.assertTrue(convex.model.number(moment(p)).log() < bound)

    def test_point_has_no_curvature_penalty(self):
        ctx.prec = 256
        value = arb(3).log()
        bound = convex.interpolate(value,value,F(1,3),F(1,3),131072)
        self.assertTrue(abs(bound-value) < arb('1e-60'))


if __name__ == '__main__':
    unittest.main()
