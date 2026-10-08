from fractions import Fraction as F
from math import comb
import unittest

from gf8_punctured_middle_certificate import (WITNESSES,exact_transfer,source_upper,
                                             cell_upper,prefactor,certify)
from gf8_punctured_profile import metadata
from gf8_punctured_screen import packet
from return_cover_certificate import Witness,interval


class PuncturedMiddleTests(unittest.TestCase):
    def test_exact_transfer_matches_screen(self):
        import math
        for _,witness in WITNESSES[::3]:
            exact=exact_transfer(witness)
            floating=packet([math.log(getattr(witness,key)) for key in ('x','y','tau')],8,7,8)
            for i in range(2):
                for j in range(2):
                    self.assertAlmostEqual(float(exact[i][j]),floating[i,j],delta=1e-10*max(1.,floating[i,j]))

    def test_source_cell_bound_exact_samples(self):
        for lower,upper in ((F(1,50),F(3,100)),(F(1,3),F(3,8)),(F(2,5),F(9,20))):
            for _,witness in (WITNESSES[0],WITNESSES[-1]):
                bound=source_upper(lower,upper,witness)
                for alpha in (lower,(lower+upper)/2,upper):
                    beta=witness.tau*witness.tilt
                    exact=(1-alpha)**128+sum(F(comb(128,h))*beta**h*alpha**h*(1-alpha)**(128-h)
                                            for h in range(45,129))/witness.tau**44
                    self.assertLessEqual(exact,F(bound.hi))

    def test_cell_dominates_singleton_sum(self):
        g=metadata()
        witness=WITNESSES[0][1]
        factor=prefactor(witness,g)
        lo,hi=32113,32116
        envelope=cell_upper(lo,hi,g,witness,factor)
        total=sum((cell_upper(w,w,g,witness,factor) for w in range(lo,hi+1)),interval(F(0)))
        self.assertLessEqual(total.hi,envelope.hi)

    def test_reject_inexact_witness(self):
        with self.assertRaises(ValueError):
            exact_transfer(Witness(.5,F(1,2),F(5),F(1)))
        with self.assertRaises(ValueError):
            certify(.02,F(9,20))

    def test_directed_complete_replay(self):
        result=certify()
        self.assertEqual(result['weights'],[32113,722534])
        self.assertTrue(result['complete'] and result['rigorous'])
        self.assertFalse(result['whole_code'])
        self.assertGreaterEqual(result['margin_bits_lower'],40)


if __name__=='__main__':
    unittest.main()
