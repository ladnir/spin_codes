"""Check the length adapter against frozen finite producers."""
from fractions import Fraction as F
import unittest
from flint import ctx
import ladder


class LengthAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 128

    def test_existing_identities(self):
        for exponent in (16, 18, 20):
            self.assertEqual(ladder.Engine(exponent).identity(),
                             ladder.candidate.Engine(exponent).identity())

    def test_extended_dimensions(self):
        for exponent in (22, 24):
            engine = ladder.Engine(exponent)
            checker = ladder.Checker(exponent)
            self.assertEqual(engine.output_bits, 1 << (exponent + 1))
            self.assertEqual(checker.rows, engine.length)
            self.assertEqual(checker.cutoff, engine.cutoff)
            self.assertEqual(1 << checker.power, engine.output_bits // 128)
            self.assertEqual(engine.identity()['inner'], ladder.Engine(20).identity()['inner'])

    def test_q1_same_finite_arithmetic(self):
        old = ladder.candidate.Engine(16)
        new = ladder.Engine(16)
        expected, actual = old.q1(F(-7)), new.q1(F(-7))
        self.assertEqual(len(expected), len(actual))
        for a, b in zip(expected, actual):
            self.assertEqual(a.mid(), b.mid())
            self.assertEqual(a.rad(), b.rad())

    def test_dense_same_finite_arithmetic(self):
        old = ladder.candidate.Checker(16)
        new = ladder.Checker(16)
        self.assertEqual(old._fixed_moment(-100, F(1, 5)),
                         new._fixed_moment(-100, F(1, 5)))

    def test_reject_unknown_length(self):
        with self.assertRaises(ValueError):
            ladder.Engine(17)


if __name__ == '__main__':
    unittest.main()
