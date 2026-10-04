import unittest
from itertools import product
from math import comb
import numpy as np
import fractional_screen as fs


class FractionalTests(unittest.TestCase):
    def test_thinning_normalizes(self):
        family=np.repeat(np.eye(2)[None],17,axis=0)
        np.testing.assert_allclose(fs.thin(family),family,rtol=1e-14,atol=0)

    def test_ordered_grouping(self):
        physical=np.array([[[1.,.2],[.1,.7]],[[.4,.1],[.2,.5]],[[.1,.2],[.3,.2]]])
        # Tiny two-slot family retains the same four-bit labels.
        for steps in (1,2,4):
            for alpha in (1.,.7):
                actual,a=fs.grouped(physical,.085,alpha,steps)
                local=fs.thin(physical)/a**np.arange(3)[:,None,None]
                expected=np.zeros_like(actual)
                for positions in product(range(3),repeat=steps):
                    matrix=np.eye(2)
                    weight=1
                    for j in positions:
                        matrix=matrix@local[j]
                        weight*=comb(2,j)
                    total=sum(positions)
                    expected[total]+=weight/comb(2*steps,total)*matrix**alpha
                np.testing.assert_allclose(actual,expected,rtol=2e-14,atol=1e-15)

    def test_fractional_path_bound(self):
        left=np.array([[.2,.3],[.4,.1]])
        right=np.array([[.5,.2],[.1,.6]])
        for alpha in (.1,.5,.9,1.):
            exact=float((left@right)[0].sum())**alpha
            envelope=float(((left**alpha)@(right**alpha))[0].sum())
            self.assertGreaterEqual(envelope+1e-15,exact)


if __name__=='__main__':
    unittest.main()
