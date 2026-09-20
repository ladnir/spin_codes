from fractions import Fraction as F
import unittest
import numpy as np
from flint import arb, arb_mat, ctx
import split_imt as split


class SplitTests(unittest.TestCase):
    def test_seven_state_terminal_encloses_arb(self):
        ctx.prec = 256
        rng = np.random.default_rng(19)
        matrices = rng.integers(0, 32, size=(8, 7, 7)).astype(float) / 128
        exponents = np.arange(-4, 4, dtype=np.int64)
        values, powers = split.terminal_256(matrices, exponents)
        for matrix, exponent, upper, power in zip(matrices, exponents, values, powers):
            exact = arb_mat([[arb(float(v)) for v in row] for row in matrix])**256
            value = sum((exact[0, j] for j in range(7)), arb(0)) * arb(2)**(256 * int(exponent))
            bound = arb(float(upper)) * arb(2)**int(power)
            self.assertLess(value, bound)
            self.assertLess(bound / value, arb('1.000000001'))

    def test_fold_all_offsets_against_scalar_arb(self):
        ctx.prec = 256
        region = [tuple(arb(i + j + 1) / 128 for j in range(49)) for i in range(5)]
        left, right = np.array([.25, .5]), np.array([.75, .25])
        current = region
        for d, matrices, powers in split.fold.folds(*split.sparse_ranges.initial(region, 7), left, right):
            for h, flat in enumerate(current):
                for k, exact in enumerate(flat):
                    bound = arb(float(matrices[h].flat[k])) * arb(2)**int(powers[h])
                    self.assertLessEqual(exact.upper(), bound)
            current = [tuple(max((arb(float(x)) * a[k] + arb(float(y)) * b[k]).upper()
                                 for x, y in zip(left, right)) for k in range(49))
                       for a, b in zip(current[:-1], current[1:])]

    def test_reject_four_state_terminal(self):
        with self.assertRaises(AssertionError):
            split.terminal_256(np.ones((1, 4, 4)), np.zeros(1, dtype=np.int64))


if __name__ == '__main__':
    unittest.main()
