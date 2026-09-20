import itertools
import unittest
from fractions import Fraction as F
from flint import arb,ctx
import bridge
import convex_region_majorant as convex


def expectation(sequence,ps):
    return sum((sequence[sum(bits)]*__import__('math').prod(p if bit else 1-p for bit,p in zip(bits,ps))
                for bits in itertools.product((0,1),repeat=len(ps))),F(0))


class ConvexTests(unittest.TestCase):
    def test_dominating_convex_sequence(self):
        ctx.prec=128
        original=[tuple(arb(v+i) for i in range(9)) for v in (0,5,2,8,1,0)]
        result=convex.convexify(original);convex.assert_convex(result)
        self.assertTrue(all(a>=b for x,y in zip(result,original) for a,b in zip(x,y)))

    def test_every_bernoulli_mixture(self):
        sequence=[F(17),F(8),F(4),F(2),F(1)]
        self.assertTrue(all(a-2*b+c>=0 for a,b,c in zip(sequence,sequence[1:],sequence[2:])))
        for ps in itertools.product((F(0),F(1,5),F(3,5),F(1)),repeat=4):
            mean=sum(ps)/4
            self.assertLessEqual(expectation(sequence,ps),expectation(sequence,[mean]*4))


if __name__=='__main__':unittest.main()
