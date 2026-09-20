from collections import Counter
from fractions import Fraction as F
from types import SimpleNamespace
import unittest

from flint import arb, ctx
import bernoulli_activation as activation
import fixed_input


class NewTransferTests(unittest.TestCase):
    def test_bernoulli_activation_pointwise_on_embedded_toy(self):
        ctx.prec = 256
        columns = [1, 2, 4, 3, 5, 6]
        spectrum = {1: 3, 2: 3, 3: 1}
        b_spectrum = Counter(sum((a & c).bit_count() % 2 for c in columns)
                             for a in range(1, 8))
        engine = SimpleNamespace(m=7, levels=sorted(spectrum), spectrum=spectrum,
                                 b_spectrum=b_spectrum)
        # The other 122 coordinates have zero feedback. They contribute
        # independently to the emitted weight, keeping the 128-bit formula.
        for theta in (F(1, 10), F(1, 2), F(9, 10)):
            for z in (F(1, 4), F(3, 4)):
                actual = [F(0)] * 8
                for x in range(64):
                    syndrome = 0
                    for j, column in enumerate(columns):
                        if x >> j & 1:
                            syndrome ^= column
                    w = x.bit_count()
                    actual[syndrome] += (theta * z)**w * (1 - theta)**(6 - w)
                tail = (1 - theta + theta * z)**122
                actual = [v * tail for v in actual]
                matrix = [arb(0)] * 25
                matrix[0] = activation.model.number(actual[0])
                result = activation.activated(matrix, engine, activation.model.number(theta),
                                               -activation.model.number(z).log())
                self.assertIs(result[0], matrix[0])
                self.assertEqual(result[1], arb(0))
                for q in range(1, 8):
                    shell = q.bit_count()
                    cap = result[engine.levels.index(shell) + 2] / spectrum[shell]
                    self.assertLessEqual(activation.model.number(actual[q]), cap)

    def test_fixed_input_evaluation_and_uniform_moment(self):
        ctx.prec = 256
        c = fixed_input.Checker(16)
        q, a = 218, F(21, 128)
        # Fixed, deterministic witness: discovery is not part of the test.
        witness = dict(tilt=-4, input_tilt=0,
                       eta=activation.model.base.encode(F(0)),
                       fixed_r=activation.model.base.encode(F(1, 10)))
        bound = c.fixed_input_bound(q, q, a, a, witness)
        self.assertTrue(bound.is_finite())
        # Enlarging only the moment's input interval must dominate its
        # complete point moments; this checks the matrix interval path.
        p = arb('0.100 +/- 0.001')
        lam = (arb(-4) / 40).exp()
        matrix = c.engine.bernoulli(p, lam)
        n, epochs = c.engine.n, 1 << c.power
        upper = activation.model.independent.terminal(matrix, n, epochs)
        for r in (F(99, 1000), F(1, 10), F(101, 1000)):
            point = c.engine.bernoulli(activation.model.number(r), lam)
            moment = activation.model.independent.terminal(point, n, epochs)
            self.assertLess(moment, upper)


if __name__ == '__main__':
    unittest.main()
