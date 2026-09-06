import unittest
from fractions import Fraction as F
from flint import arb,ctx
import tightened_occupancy as tight
import polynomial_regions as poly


class PolynomialTests(unittest.TestCase):
    def test_every_toy_coefficient(self):
        ctx.prec=256
        exact=tight.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},F(3,4),F,32,length=32)
        result=poly.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},arb(3)/4,32,length=32)
        from audit_bch_q1_full_arb import rational
        for row,want in zip(result,exact):
            for value,truth in zip(row,want):
                upper=rational(value)
                self.assertGreaterEqual(upper,truth)
                self.assertLessEqual(upper-truth,F(1,1<<180))
        truncated=poly.regions(8,4,{4:14,8:1},{0:1,4:14,8:1},arb(3)/4,9,length=32)
        for row,want in zip(truncated,exact):
            for value,truth in zip(row,want):self.assertGreaterEqual(rational(value),truth)


if __name__=='__main__':unittest.main()
