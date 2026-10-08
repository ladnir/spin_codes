from fractions import Fraction as F
import unittest

from flint import arb, arb_poly, ctx

from cauchy_firstmoment_screen import geometry
from field128_sparse_certificate import complete_cover, region, certify_zone
from joint_outer_screen import WITNESSES_128_87
from state8_sparse_certificate import regional


class ReplayTests(unittest.TestCase):
    def test_geometry_and_overlap(self):
        g = geometry(20,128,87)
        self.assertEqual((g['rows'],g['k'],g['n'],g['threshold']),
                         (12096,1052352,1548288,154828))
        self.assertGreaterEqual(42*1106,(3*g['n']+99)//100)
        self.assertLess(42*1105,(3*g['n']+99)//100)

    def test_complete_cover_rejects_gap_overlap(self):
        good = [dict(lower=1,upper=3),dict(lower=4,upper=8)]
        complete_cover(good,1,8)
        for second in (3,5):
            with self.assertRaises(ValueError):
                complete_cover([good[0],dict(lower=second,upper=8)],1,8)

    def test_binomial_recurrence(self):
        old = ctx.prec
        try:
            ctx.prec = 256
            polynomials = tuple(arb_poly([1,2,3,4,5]) for _ in range(4))
            a = region(polynomials,4,64,F(7,3))
            b = regional(polynomials,4,64,F(7,3))
            for x,y in zip(a,b):
                self.assertTrue(x.overlaps(y))
        finally:
            ctx.prec = old

    def test_direct_singleton_receipt(self):
        upper,receipt = certify_zone(1,1,1,WITNESSES_128_87[1],
                                     geometry(20,128,87),precision=512)
        self.assertLess(upper,arb(2)**-200)
        self.assertEqual(receipt['entries'][0]['lower'],1)
        self.assertEqual(receipt['certified_integer_bits'],208)


if __name__ == '__main__':
    unittest.main()
