from fractions import Fraction
from math import comb
import unittest

from flint import arb, ctx

import two_epoch_return as candidate


class TwoEpochReturnTests(unittest.TestCase):
    def setUp(self):
        self.columns=[1,2,4,8]
        # Identity expansion; two disjoint two-bit windows.
        self.images=list(range(16))
        self.data=candidate.census(self.images,self.columns,width=2)

    def test_exact_zero_tilt_return_probability(self):
        for rounds in (1,2,3):
            values=candidate.moments(self.data,15,'0',rounds=rounds,width=2)
            alpha=Fraction(1,2**rounds)
            for (a,b),value in values.items():
                exact=(alpha/Fraction(2*comb(2,a)) if a==b else 0)+(1-alpha)/15
                self.assertTrue((value*exact.denominator-exact.numerator).contains(0))

    def test_direct_state_kernel_matches_histogram_moment(self):
        ctx.prec=192
        tilt=arb('.12');rho=arb('.8');alpha=arb(1)/4
        values=candidate.moments(self.data,15,'.12','.8',2,width=2)
        words={a:[x for x in range(1,16) if x.bit_count()==a and (x<4 or x%4==0)]
               for a in (1,2)}
        for (a,b),value in values.items():
            total=arb(0)
            for x in words[a]:
                for y in words[b]:
                    # Enumerate every refreshed state, independently of the
                    # equal-feedback histogram used by the candidate.
                    for refreshed in range(1,16):
                        probability=(alpha if refreshed==x else 0)+(1-alpha)/15
                        if refreshed^y==0:
                            total+=probability*(-tilt*(x.bit_count()+(x^y).bit_count())).exp()
            total*=rho**(int(a==2)+int(b==2))/(len(words[a])*len(words[b]))
            self.assertTrue(value.overlaps(total))

    def test_feedback_injectivity_is_checked(self):
        for columns in ([1,1,4,8],[0,2,4,8]):
            with self.assertRaisesRegex(ValueError,'globally injective'):
                candidate.census(self.images,columns,width=2)


if __name__=='__main__':unittest.main()
