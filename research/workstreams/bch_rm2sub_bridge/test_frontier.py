"""Exact tests for the size parameterization and continuous-density enclosure."""
from fractions import Fraction as F
import itertools
import unittest
from flint import arb,ctx
import frontier_dense as dense
from frontier_sparse import parameters


class FrontierTest(unittest.TestCase):
    def test_size_parameters(self):
        for m in (20,24,30):
            rows,cutoff=parameters(m)
            self.assertEqual(128*rows,1<<m)
            self.assertEqual(cutoff,(1<<(m+1))//10)
            self.assertEqual((rows*256)//64,1<<(m-5))

    def test_exact_monotone_decomposition(self):
        for co in itertools.product(map(F,(0,1,3)),repeat=5):
            a,d=dense.monotone_parts(co)
            self.assertEqual([x-y for x,y in zip(a,d)],list(co))
            self.assertTrue(all(x<=y for x,y in zip(a,a[1:])))
            self.assertTrue(all(x<=y for x,y in zip(d,d[1:])))

    def test_exact_binomial_interval_enclosure(self):
        arrays=((0,1,4,2,0),(5,4,3,2,1),(0,0,0,0,7),(7,0,0,0,0))
        for raw in arrays:
            co=list(map(F,raw));a,d=dense.monotone_parts(co)
            for lo,hi in ((F(0),F(1)),(F(1,7),F(3,7)),(F(1,2),F(1,2)),(F(9,10),F(1))):
                upper=(sum((w*v for w,v in zip(dense.bernstein_weights(4,hi),a)),F(0))-
                    sum((w*v for w,v in zip(dense.bernstein_weights(4,lo),d)),F(0)))
                for j in range(9):
                    r=lo+(hi-lo)*F(j,8)
                    exact=sum((w*v for w,v in zip(dense.bernstein_weights(4,r),co)),F(0))
                    self.assertLessEqual(exact,upper)
                    if lo==hi:self.assertEqual(exact,upper)

    def test_arb_interval_encloses_exact_expression(self):
        ctx.prec=256
        co=list(map(F,(0,4,1,7,2)));a,d=dense.monotone_parts(co)
        lo,hi=F(3,11),F(5,11)
        exact=(sum((w*v for w,v in zip(dense.bernstein_weights(4,hi),a)),F(0))-
            sum((w*v for w,v in zip(dense.bernstein_weights(4,lo),d)),F(0)))
        got=(sum((w*dense.aa(v) for w,v in zip(dense.bernstein_weights(4,dense.aa(hi)),a)),arb(0))-
            sum((w*dense.aa(v) for w,v in zip(dense.bernstein_weights(4,dense.aa(lo)),d)),arb(0)))
        self.assertLessEqual(dense.rational(got.lower()),exact)
        self.assertGreaterEqual(dense.rational(got.upper()),exact)

    def test_exact_cost_hull_range_maximum(self):
        hull=[(F(1,4),F(1)),(F(1,2),F(3)),(F(1),F(0))]
        self.assertEqual(dense.hull_max(hull,F(3,8),F(3,4)),3)
        self.assertEqual(dense.hull_max(hull,F(3,4),F(1)),F(3,2))


if __name__=='__main__':unittest.main()
