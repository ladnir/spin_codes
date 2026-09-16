import itertools
import unittest
import numpy as np
from flint import arb, arb_mat, ctx
import holder_split as holder


class HolderTests(unittest.TestCase):
    def test_norm_encloses_high_precision(self):
        ctx.prec = 256
        rng = np.random.default_rng(256)
        values = rng.random((12, 10, 7, 7))
        upper = holder.norm256(values)
        for index in ((0, 0, 0), (3, 4, 5), (9, 6, 6)):
            exact = sum((arb(float(v))**256 for v in values[(slice(None),) + index]), arb(0)).root(256)
            self.assertLess(exact, arb(float(upper[index])))
            self.assertLess(arb(float(upper[index])) / exact, arb('1.000000000001'))

    def test_band_sequence_sum_enclosed(self):
        ctx.prec = 256
        rng = np.random.default_rng(23)
        region = [tuple(arb(float(v)) for v in (rng.integers(1, 17, 49) / 256)) for _ in range(4)]
        left, right = np.array([.25, .5, .75]), np.array([.5, .125, .25])
        for d, matrices, powers in holder.folds(*holder.base.sparse_ranges.initial(region, 7), left, right):
            h = 3 - d
            exact_sum = arb(0)
            for sequence in itertools.product(range(3), repeat=d):
                current = region
                for g in sequence:
                    current = [tuple(arb(float(left[g])) * a[k] + arb(float(right[g])) * b[k]
                                     for k in range(49)) for a, b in zip(current[:-1], current[1:])]
                matrix = arb_mat([list(current[h][7*i:7*i+7]) for i in range(7)])**256
                exact_sum += sum((matrix[0, j] for j in range(7)), arb(0))
            value, exponent = holder.base.terminal_256(matrices[h:h+1], powers[h:h+1])
            bound = arb(float(value[0])) * arb(2)**int(exponent[0])
            self.assertLess(exact_sum, bound)


if __name__ == '__main__':
    unittest.main()
