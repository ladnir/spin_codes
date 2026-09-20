import unittest
from fractions import Fraction as F
import numpy as np
from flint import arb
import scaled_adaptive as scaled
from test_adaptive_range import exact_adaptive


class ScaledTests(unittest.TestCase):
    def test_extreme_scales_against_exact(self):
        for exps in ([0,-500,-1000,-1500,-2000],[-3000,-2600,-2200,-1800,-1400]):
            exact=[tuple(F(k+1,16)*F(2)**e for k in range(9)) for e in exps]
            regions=[tuple((arb(v.numerator)/v.denominator).upper() for v in row) for row in exact]
            m,e=scaled.initial(regions)
            ps=[F(1,4),F(3,4)];roots=[F(5,4),F(7,4)]
            left=np.nextafter(np.array([float(r*(1-p)) for r,p in zip(roots,ps)]),np.inf)
            right=np.nextafter(np.array([float(r*p) for r,p in zip(roots,ps)]),np.inf)
            for q,matrix,exponent in scaled.matrices(m,e,left,right):
                expected=exact_adaptive(exact[:q+1],ps,roots)
                for v,w in zip(matrix.flat,expected):
                    bound=F.from_float(float(v))*F(2)**exponent
                    self.assertGreaterEqual(bound,w)
                    self.assertLess(bound-w,w*F(1,1<<40))


if __name__=='__main__':unittest.main()
