"""Independent discrete verification of proposed quadratic witnesses."""
from fractions import Fraction as Q
import unittest
from flint import arb,ctx

from quadratic import witness


class QuadraticWitness(unittest.TestCase):
    def setUp(self):self.old=ctx.prec;ctx.prec=192
    def tearDown(self):ctx.prec=self.old

    def test_every_discrete_overlap_constraint(self):
        for W in (3,4,20,28,40):
            classes={48:Q(1,16),56:Q(5,16),64:Q(1,2)}
            value,c,intercepts=witness(classes,W,Q(W),'.072')
            ctx.prec=512
            slope=arb(c.numerator)/c.denominator
            for v,a in intercepts.items():
                for h in range(W+1):
                    d=2*h-W
                    actual=(-arb('.072')*(v-d)).exp()
                    self.assertGreaterEqual(a+slope*d*d,actual)
            expected=sum((arb(p.numerator)/p.denominator*intercepts[v] for v,p in classes.items()),arb(0))+slope*W
            self.assertGreaterEqual(value,expected)
            ctx.prec=192

    def test_direct_distributions_with_same_second_moment(self):
        W=4;classes={6:Q(1,4),8:Q(3,4)}
        cases=[[(6,0,Q(1,4)),(8,2,Q(1,2)),(8,4,Q(1,4))],
               [(6,3,Q(1,4)),(8,1,Q(1,2)),(8,2,Q(1,4))]]
        for law in cases:
            second=sum((p*(2*h-W)**2 for v,h,p in law),Q(0))
            value,_,_=witness(classes,W,second,'.17')
            ctx.prec=512
            actual=sum((arb(p.numerator)/p.denominator*(-arb('.17')*(v+W-2*h)).exp()
                        for v,h,p in law),arb(0))
            self.assertGreaterEqual(value,actual)
            ctx.prec=192

    def test_reject_nonexact_or_unsupported_inputs(self):
        for classes,W,second in (({2:Q(1)},4,Q(4)),({48:Q(2)},4,Q(4)),
                                 ({48:Q(1)},4,4.0),({48:Q(1)},4,Q(-1))):
            with self.assertRaises(ValueError):witness(classes,W,second,'.072')


if __name__=='__main__':unittest.main()
