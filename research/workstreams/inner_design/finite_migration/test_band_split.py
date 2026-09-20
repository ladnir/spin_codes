"""Exact subset-counting checks for the high-band refinement."""
from fractions import Fraction as F
from itertools import product
import unittest

from flint import arb,ctx
import band_split_point
import high_band_dense
import retune_fixed_input


class BandSplitTests(unittest.TestCase):
    def test_marked_subset_sum(self):
        weights = (F(2,7),F(3,11),F(1,19),F(1,31))
        for q in range(1,6):
            low,high = F(0),F(0)
            for indices in product(range(4),repeat=q):
                mass = F(1)
                for i in indices:
                    mass *= weights[i]
                if any(i >= 2 for i in indices):
                    high += mass
                else:
                    low += mass
            self.assertEqual(low,sum(weights[:2])**q)
            self.assertEqual(low+high,sum(weights)**q)
            # exp(-xi)=z>=1. Every high assignment receives a multiplier
            # z^h/z >=1 in the marked counting sum.
            for z in (F(1),F(2),F(7)):
                self.assertLessEqual(high,(sum(weights[:2])+z*sum(weights[2:]))**q/z)

    def test_point_and_rectangle_backend_agree(self):
        ctx.prec = 256
        c = high_band_dense.Checker(16)
        q,v = 205,F(201,1024)
        w = retune_fixed_input.proposal(c,q,v,-12,-0.3)
        point,*_ = band_split_point.bound(c,q,v,w,high_band_dense.CUTOFF,high_band_dense.XI)
        box = c.split_bound(q,q,v,v,w)
        self.assertTrue(abs(point-box) < arb('1e-50'))

    def test_combined_retains_each_complete_bound(self):
        import combined_dense
        ctx.prec = 192
        c = combined_dense.Checker(16)
        q,v = 205,F(201,1024)
        w = retune_fixed_input.proposal(c,q,v,-12,-0.3)
        actual = c.bound(q,q,v,v,w)
        high = high_band_dense.Checker.bound(c,q,q,v,v,w)
        zero = c.zero_checker.bound(q,q,v,v,w)
        self.assertTrue(actual <= high and actual <= zero)

    def test_pruned_density_agrees(self):
        import fast_density
        import constant_density
        ctx.prec = 192
        for rows,lo,hi,nlo,nhi,ps in (
                (16,4,6,F(1,3),F(2,3),(F(1,4),F(1,2),F(3,4))),
                (512,205,205,F(2,5),F(2,5),(F(1,4),F(1,2),F(3,4)))):
            expected = constant_density.factor(rows,lo,hi,nlo,nhi,ps,F(0))
            actual = fast_density.factor(rows,lo,hi,nlo,nhi,ps,F(0))
            self.assertTrue(abs(actual-expected) < arb('1e-40'))

    def test_pruned_unrestricted_density(self):
        import fast_density
        import density_comparison
        ctx.prec = 192
        ps = (F(1,4),F(1,2),F(3,4))
        for lo,hi,nlo,nhi in ((205,205,F(2,5),F(2,5)),(203,207,F(39,100),F(41,100))):
            maximum = (nhi-min(ps))/(1-min(ps))
            expected = density_comparison.factor(512,lo,hi,nlo,nhi,ps)
            actual = fast_density.factor(512,lo,hi,nlo,nhi,ps,maximum)
            self.assertTrue(abs(actual-expected) < arb('1e-40'))


if __name__ == '__main__':
    unittest.main()
