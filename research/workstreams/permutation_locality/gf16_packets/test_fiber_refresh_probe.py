import unittest
import numpy as np
import rank_return
import fiber_refresh_probe as probe
from fractions import Fraction as Q


class FiberRefreshProbeTests(unittest.TestCase):
    def test_arbitrary_feedback_targets_and_domination(self):
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        columns=[1,2,4,3,5,7,6,1]
        data=probe.prepare(rank_return.prepare(images,columns,3,2));z=.75;L=7
        for p in (0.,.2,1.):
            actual=np.zeros((8,8));fibers=np.zeros((8,8))
            for a in range(8):
                for x in range(256):
                    probability=(p/15 if x&15 else 1-p)*(p/15 if x>>4 else 1-p)
                    weighted=probability*z**(images[a]^x).bit_count();feedback=0
                    for b,c in enumerate(columns):
                        if x>>b&1:feedback^=c
                    fibers[a,feedback]+=weighted
                    if not a:actual[a,feedback]+=weighted
                    else:
                        actual[a,a^feedback]+=.25*weighted
                        for v in range(1,8):actual[a,v^feedback]+=.75/L*weighted
            maximum,mean=probe.density(data,p,z)
            self.assertGreaterEqual(maximum+1e-12,fibers[1:].max())
            for target in range(8):
                shifted=sum(fibers[a,a^target] for a in range(1,8))/L
                self.assertGreaterEqual(mean+1e-12,shifted)
            original=rank_return.floating(data,probe.sc.probabilities(Q(p)),-np.log(z))
            for matrix in probe.candidates(data,p,z,original).values():
                for a in range(8):
                    exact=np.zeros(8);exact[a]=1
                    envelope=np.array([float(a==0),float(a!=0),0.])
                    for _ in range(6):
                        exact=exact@actual;envelope=envelope@matrix
                        tolerance=1e-10*max(1.,envelope.sum())
                        self.assertLessEqual(exact[0],envelope[0]+tolerance)
                        residual=np.maximum(0.,exact[1:]-envelope[2]/L).sum()
                        self.assertLessEqual(residual,envelope[1]+tolerance)


if __name__=='__main__':unittest.main()
