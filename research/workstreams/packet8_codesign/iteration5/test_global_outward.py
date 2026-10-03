import unittest
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import outward_positive as op
import scaled_positive as sp
import global_outward as go


def exact_matrix(value):
    return np.array([[F(float(x))*F(2)**value.exponent for x in row] for row in value.value], dtype=object)


class GlobalTests(unittest.TestCase):
    def test_scaled_extremes(self):
        for e in (-10000,0,10000):
            x=sp.from_fraction(F(3,7)*F(2)**e)
            y=sp.from_fraction(F(1,9)*F(2)**(e-2000))
            self.assertGreaterEqual(sp.scalar_fraction(sp.add(x,y)),F(3,7)*F(2)**e+F(1,9)*F(2)**(e-2000))
            product=sp.multiply(x,y)
            self.assertGreaterEqual(sp.scalar_fraction(product),F(1,21)*F(2)**(2*e-2000))
            for a in (F(1,2),F(7,20)):
                bound=sp.scalar_fraction(sp.fractional_power(x,a))
                self.assertGreaterEqual(bound**a.denominator,sp.scalar_fraction(x)**a.numerator)

    def test_underflow_matrix_support(self):
        a=np.array([[1.,op.ETA],[0.,op.ETA]])
        b=np.array([[1.,0.],[op.ETA,op.ETA]])
        upper=exact_matrix(sp.matmul(sp.Scaled(a),sp.Scaled(b)))
        exact=np.array([[F(float(x)) for x in row] for row in a],dtype=object)@np.array([[F(float(x)) for x in row] for row in b],dtype=object)
        self.assertTrue(np.all(upper>=exact))
        self.assertGreater(upper[0,1],0)

    def test_placement_exact_small(self):
        local=np.array([[[1.,0.],[0.,.75]],[[.25,.5],[.125,.5]],[[.0625,.25],[.25,.25]]])
        matrices=[np.array([[F(float(x)) for x in row] for row in m],dtype=object) for m in local]
        computed=go.placement(local,4,epochs=2,windows=2)
        for q in range(5):
            exact=np.full((2,2),F(0),dtype=object)
            supports=list(combinations(range(4),q))
            for support in supports:
                k=sum(i<2 for i in support)
                exact+=matrices[k]@matrices[q-k]
            exact/=len(supports)
            bound=exact_matrix(computed[q])
            self.assertTrue(np.all(bound>=exact))
            self.assertLess(float(np.max(bound-exact)),1e-12)
            # All regions keep the same state. No scalar/reset replacement.
            moment=go.continuous_moment(computed[q],3)
            self.assertGreaterEqual(sp.scalar_fraction(moment),sum((exact@exact@exact)[0]))

    def test_prefactor_exact(self):
        regional=sp.Scaled(np.array([[.125,.25],[.0625,.5]]))
        for alpha in (F(1),F(2,5),F(1,2)):
            result=go.occupancy_upper(regional,2,F(3,4),alpha,groups=3,regions=2,cutoff=5,beta=F(5,4),a=F(2,3))
            matrix=exact_matrix(regional)
            moment=sum((matrix@matrix)[0])
            bound=sp.scalar_fraction(result)/(3*moment)
            exact=F(5,4)**2*F(4,3)**5*F(2,3)**4
            self.assertGreaterEqual(bound**alpha.denominator,exact**alpha.numerator)


if __name__=='__main__':unittest.main()
