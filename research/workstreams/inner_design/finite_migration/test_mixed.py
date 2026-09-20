from fractions import Fraction as F
import unittest
from flint import arb, ctx
import mixed_dense
import ladder


class MixedWitnessTests(unittest.TestCase):
    def test_occupation_proposal_rescues_low_density_point(self):
        ctx.prec = 128
        checker = mixed_dense.Checker(22)
        args = (512, 512, F(1, 4), F(1, 4))
        direct = ladder.Checker.witness(checker, *args)
        both = checker.witness(*args)
        self.assertGreater(checker.bound(*args, direct), arb(0))
        self.assertLess(checker.bound(*args, both), -1000 * arb(2).log())


if __name__ == '__main__':
    unittest.main()
