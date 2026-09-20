"""Check the causal XOR moment against exhaustive tiny encoders."""
import math
import unittest

import parameter_causal_dense as causal


class CausalTests(unittest.TestCase):
    def test_exhaustive_causal_maps(self):
        # Two-bit epochs, one-bit state; nonlinear feedback is also allowed
        # by this bound, since only output-before-update is needed.
        length = 6
        for expansion in (0,1,2,3):
            for rule in (0,3,5,11,15):
                for theta in (.05,.25,.5,.8,.99):
                    for lam in (.1,1.,3.):
                        moment = 0.
                        for x in range(1 << length):
                            state,weight = 0,0
                            for offset in range(0,length,2):
                                chunk = (x >> offset) & 3
                                weight += (chunk ^ (expansion if state else 0)).bit_count()
                                state ^= (rule >> chunk) & 1
                            probability = theta**x.bit_count()*(1-theta)**(length-x.bit_count())
                            moment += probability*math.exp(-lam*weight)
                        upper = math.exp(causal.causal_moment(theta,lam,length))
                        self.assertLessEqual(moment,upper+1e-14)
                        if theta == .5:
                            self.assertAlmostEqual(moment,upper,places=13)

    def test_reject_invalid_inputs(self):
        for args in ((-.1,1,1),(1.1,1,1),(.5,-1,1),(.5,1,-1)):
            with self.assertRaises(ValueError):
                causal.causal_moment(*args)


if __name__ == '__main__':
    unittest.main()
