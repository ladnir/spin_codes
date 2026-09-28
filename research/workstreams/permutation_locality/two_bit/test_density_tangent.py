from fractions import Fraction as Q
from itertools import product,permutations
from math import factorial
import unittest

from flint import arb,ctx

import density_tangent
from heterogeneous import tilted_reference,density_loss
from probe import aq


class DensityTangentTests(unittest.TestCase):
    def test_every_small_profile_and_anchor(self):
        for n in range(1,9):
            profiles=[(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]
            for anchors in [*profiles,(0,0,0),(n,n,n)]:
                base,tau=density_tangent.exact(n,anchors)
                for profile in profiles:
                    actual=Q(n**n,factorial(n));upper=base
                    for a,t in zip(profile,tau):
                        actual*=Q(factorial(a),a**a);upper*=t**a
                    self.assertLessEqual(actual,upper)
                    if profile==anchors:self.assertEqual(actual,upper)

    def test_tilted_shuffled_law_pointwise(self):
        laws=[(Q(1),Q(0),Q(0)),(Q(1,4),Q(1,2),Q(1,4)),(Q(3,5),Q(1,10),Q(3,10))]
        for tilt in ([1,1,1],[1,Q(2,7),Q(3,11)],[1,2,3]):
            reference,factor=tilted_reference(laws,[1,1,1],tilt)
            for anchors in ((1,1,1),(3,0,0),(2,1,0),(0,0,0)):
                base,tau=density_tangent.exact(3,anchors)
                for pattern in product(range(3),repeat=3):
                    actual=Q(0)
                    for order in permutations(range(3)):
                        term=Q(1)
                        for i,j in enumerate(order):term*=laws[j][pattern[i]]
                        actual+=term/6
                    upper=factor/density_loss(3)*base
                    for j in pattern:upper*=reference[j]*tau[j]
                    self.assertGreaterEqual(upper,actual)

    def test_arb_constants_enclose_exact_rationals(self):
        ctx.prec=192
        for n,anchors in ((3,(3,0,0)),(6,(4,1,1)),(9,(0,9,9)),(32,(27,4,1))):
            base,tau=density_tangent.exact(n,anchors)
            logbase,logtau=density_tangent.outward(n,anchors)
            self.assertTrue(logbase.exp().overlaps(aq(base)))
            for logt,t in zip(logtau,tau):self.assertTrue(logt.exp().overlaps(aq(t)))
            floating,_=density_tangent.floating(n,anchors)
            self.assertAlmostEqual(floating,float(logbase),places=11)
        for bad in ((3,()),(3,(4,0,0)),(3,(1.5,1,0)),(0,(0,0,0))):
            with self.assertRaises(ValueError):density_tangent.outward(*bad)


if __name__=='__main__':unittest.main()
