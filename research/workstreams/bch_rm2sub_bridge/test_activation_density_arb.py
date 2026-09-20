import math
import unittest

import numpy as np
from flint import arb, ctx

import activation_density_arb as engine


class DensityArbTests(unittest.TestCase):
    def test_epochs_agree_with_independent_log_engine(self):
        ctx.prec = 256
        for t, s, ac, kernel in ((4, 2, {2: 2, 4: 1}, {0: 1, 2: 2, 4: 1}),
                                (8, 4, {4: 14, 8: 1}, {0: 1, 4: 14, 8: 1})):
            model = engine.reference.Epochs(t, s, ac, [kernel.get(j, 0) for j in range(t+1)])
            for lam in (.15, .6, 2.):
                got = np.array([[float(v) for v in row] for row in
                    engine.epoch(t, s, ac, kernel, (-arb(lam)).exp())]).reshape(-1, 4, 4)
                np.testing.assert_allclose(got, np.exp(model.at(lam)), rtol=3e-14, atol=1e-16)

    def test_regions_match_direct_positive_product(self):
        ctx.prec = 256
        args = (4, 2, {2: 2, 4: 1}, {0: 1, 2: 2, 4: 1}, arb(3)/5)
        epoch = engine.epoch(*args)
        region = engine.regions(*args, 8, 8)
        for j in range(9):
            direct = [arb(0) for _ in range(16)]
            for a in range(max(0, j-4), min(4, j)+1):
                product = engine.mul(epoch[a], epoch[j-a])
                for k in range(16):
                    direct[k] += product[k]*math.comb(4, a)*math.comb(4, j-a)/math.comb(8, j)
            for actual, expected in zip(region[j], direct):
                self.assertLess(abs(float(actual-expected)), 1e-70)

    def test_bad_kernel_rejected(self):
        with self.assertRaises(ValueError):
            engine.epoch(4, 2, {2: 2, 4: 1}, {0: 1, 1: 2, 4: 1}, arb(1)/2)


if __name__ == '__main__':
    unittest.main()
