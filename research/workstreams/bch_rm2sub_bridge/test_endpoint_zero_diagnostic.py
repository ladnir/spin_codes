import math
import unittest
from fractions import Fraction as F
import numpy as np
from flint import arb,ctx
import endpoint_zero_diagnostic as diag


class ZeroTests(unittest.TestCase):
    def test_polynomial(self):
        ctx.prec=256;t=4;length=8;kernel={0:1,4:1};z=F(3,4)
        values=diag.zero_region(t,kernel,arb(3)/4,length)
        expected={0:F(1),4:2*z**4/math.comb(8,4),8:z**8}
        for j,v in enumerate(values):self.assertAlmostEqual(float(v),float(expected.get(j,0)),places=14)

    def test_scalar_recurrence(self):
        values=[F(2),F(1,2),F(1,3),F(1,7),F(1,11)]
        left=[F(1,2),F(1,4)];right=[F(1,2),F(3,4)]
        current=values
        while len(current)>1:current=[max(x*a+y*b for x,y in zip(left,right)) for a,b in zip(current,current[1:])]
        result=diag.adaptive_scalar([arb(v.numerator)/v.denominator for v in values],np.array(list(map(float,left))),np.array(list(map(float,right))))
        self.assertAlmostEqual(result,math.log(float(current[0])),places=13)


if __name__=='__main__':unittest.main()
