import unittest
from fractions import Fraction as F
import numpy as np
import shared_range_certificate as shared
from test_adaptive_range import exact_adaptive


class SharedTests(unittest.TestCase):
    def test_every_depth_dominates_exact(self):
        region=[tuple(F((j+1)*(k+2),317) for k in range(9)) for j in range(7)]
        ps=[F(1,5),F(2,3),F(9,10)];roots=[F(3,2),F(7,5),F(6,5)]
        initial=shared.up(np.array([[float(v) for v in row] for row in region])).reshape(-1,3,3)
        left=shared.up(np.array([float(r*(1-p)) for r,p in zip(roots,ps)]))
        right=shared.up(np.array([float(r*p) for r,p in zip(roots,ps)]))
        for q,matrix,exponent in shared.adaptive_matrices(initial,left,right):
            exact=exact_adaptive(region[:q+1],ps,roots)
            for v,w in zip(matrix.flat,exact):
                upper=F.from_float(float(v))*F(2)**exponent
                self.assertGreaterEqual(upper,w)
                self.assertLess(upper-w,F(1,1<<40))

    def test_positive_underflow_rounds_up(self):
        smallest=np.nextafter(0.,1.)
        self.assertGreaterEqual(F.from_float(float(shared.up(smallest*smallest))),F.from_float(smallest)**2)


if __name__=='__main__':unittest.main()
