import unittest
import numpy as np
from fractions import Fraction as Q
import birth_refresh_probe as probe


class BirthRefreshProbeTests(unittest.TestCase):
    def test_newton_coefficients_against_direct_products(self):
        for values in ([.2,-.1,.3,-.05,.4],[0.,1.,-1.,.25,-.25]):
            records=np.array([[4,0,2,1,1],[0,1,0,4,3],[2,2,2,2,0]],dtype=float)
            computed=probe.elementary(records,np.array(values),8)
            for histogram,row in zip(records,computed):
                direct=np.array([1.])
                for n,value in zip(histogram,values):
                    for _ in range(int(n)):direct=np.convolve(direct,[1.,value])
                np.testing.assert_allclose(row,direct,rtol=1e-12,atol=1e-12)

    def test_zero_state_feedback_density_and_repeated_domination(self):
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            data=probe.prepare(probe.rank_return.prepare(images,columns,3,2));z=.75;L=7
            for p in (0.,.2,1.):
                actual=np.zeros((8,8))
                for a in range(8):
                    for x in range(256):
                        probability=(p/15 if x&15 else 1-p)*(p/15 if x>>4 else 1-p)
                        weighted=probability*z**(images[a]^x).bit_count();feedback=0
                        for b,c in enumerate(columns):
                            if x>>b&1:feedback^=c
                        if not a:actual[a,feedback]+=weighted
                        else:
                            actual[a,a^feedback]+=.25*weighted
                            for v in range(1,8):actual[a,v^feedback]+=.75/L*weighted
                original=probe.rank_return.floating(data,probe.sc.probabilities(Q(p)),-np.log(z))
                for matrix in probe.candidates(data,p,z,original).values():
                    exact=np.zeros(8);exact[0]=1;envelope=np.array([1.,0.,0.])
                    for _ in range(8):
                        exact=exact@actual;envelope=envelope@matrix
                        tolerance=1e-10*max(1.,envelope.sum())
                        self.assertLessEqual(exact[0],envelope[0]+tolerance)
                        self.assertLessEqual(np.maximum(0.,exact[1:]-envelope[2]/L).sum(),envelope[1]+tolerance)


if __name__=='__main__':unittest.main()
