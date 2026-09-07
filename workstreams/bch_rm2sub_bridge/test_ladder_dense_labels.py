from fractions import Fraction as F
import itertools
import math
import unittest

from flint import arb, ctx

import ladder_dense_labels as dense


class LabelTests(unittest.TestCase):
    def test_label_sum_dominates_each_interval(self):
        ctx.prec = 256
        ps, costs, eta = [F(1,4), F(3,4)], [F(2), F(3)], F(-5,2)
        log_s = dense.log_label_sum(ps, costs, eta)
        probabilities, factors = ps+[F(1)], costs+[F(1)]
        for q in range(1, 6):
            for lo, hi in ((F(1,4),F(1,2)), (F(1,2),F(1)), (F(1),F(1))):
                actual, maximum = arb(0), arb(0)
                for labels in itertools.product(range(3), repeat=q):
                    nu = sum((probabilities[g] for g in labels), F(0))/q
                    if lo <= nu <= hi:
                        moment = dense.aa(1+nu*nu)
                        actual += dense.aa(math.prod(factors[g] for g in labels))*moment
                        maximum = max(maximum, ((dense.aa(eta*q*nu)).exp()*moment).upper())
                upper = ((q*log_s).exp()*maximum).upper()
                self.assertTrue(upper >= actual.lower())

    def test_constant_is_one_label(self):
        ctx.prec = 256
        value = dense.log_label_sum([F(1,2)], [F(7)], F(0)).exp()
        self.assertTrue(value.contains(8))

    def test_monotone_parts_bound_interior(self):
        ctx.prec = 256
        coefficients = [F(1), F(1,10), F(4,5), F(1,20), F(1,2)]
        a, d = dense.intervals.monotone_parts(coefficients)
        low = dense.intervals.bernstein_weights(4, dense.aa(F(1,5)))
        high = dense.intervals.bernstein_weights(4, dense.aa(F(4,5)))
        upper = (sum((w*dense.aa(c) for w,c in zip(high,a)), arb(0))-
                 sum((w*dense.aa(c) for w,c in zip(low,d)), arb(0))).upper()
        for r in (F(1,5), F(1,3), F(1,2), F(4,5)):
            weights = dense.intervals.bernstein_weights(4, dense.aa(r))
            actual = sum((w*dense.aa(c) for w,c in zip(weights,coefficients)), arb(0))
            self.assertTrue(upper >= actual.lower())


if __name__ == '__main__':
    unittest.main()
