from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest

from flint import arb,ctx
import numpy as np

import model
from profiles import pair_profiles,conditioned_shells
from profile_dense import float_epoch,outward_epoch


class ProfileDenseTests(unittest.TestCase):
    def setUp(self):ctx.prec=192

    def test_scalar_iid_output_against_every_small_word(self):
        v,r=Q(3,5),Q(2,7)
        tilt=Q(17,100);za=(-arb(17)/100).exp()
        probabilities=[1-v,v*(1-r)/2,v*(1-r)/2,v*r]
        for y in range(64):
            expected=arb(0)
            for packets in product(range(4),repeat=3):
                x=sum(p<<(2*i) for i,p in enumerate(packets));mass=Q(1)
                for p in packets:mass*=probabilities[p]
                expected+=(arb(mass.numerator)/mass.denominator)*za**(x^y).bit_count()
            hist=[model.pattern_histogram(y,3)]
            actual=outward_epoch(hist,v,r,tilt)
            self.assertTrue(actual>=expected.lower())
            self.assertTrue(abs(actual-expected)<arb(2)**-160)
            self.assertAlmostEqual(float_epoch(np.array(hist),float(v),float(r),float(tilt)),float(expected.log()),places=13)

    def test_uniform_bits_state_independent(self):
        hist=[(64,0,0),(0,64,0),(0,0,64),(20,24,20)]
        z=(-arb(17)/100).exp();expected=((1+z)/2)**128
        actual=outward_epoch(hist,Q(3,4),Q(1,3),Q(17,100))
        self.assertTrue(abs(actual/expected-1)<arb(2)**-160)

    def test_profile_conditioning_exhaustive(self):
        n=4;p,r=Q(2,5),Q(3,7)
        profiles=pair_profiles([1,2,3,2,1]);caps=conditioned_shells(profiles,n,r)
        # Direct coefficients: the reference pair profile is multinomial.
        for u in range(n+1):
            actual=Q(0)
            for b in range(u+1):
                count=profiles.get((u,b),Q(0))
                mass=comb(n,u)*comb(u,b)*(1-p)**(n-u)*(p*(1-r))**(u-b)*(p*r)**b
                actual+=count/mass
            self.assertEqual(caps[u]/(comb(n,u)*p**u*(1-p)**(n-u)),actual)
        for bad in (0,1,-1):
            with self.assertRaises(ValueError):conditioned_shells(profiles,n,bad)


if __name__=='__main__':unittest.main()
