from fractions import Fraction as F
import math
import unittest
import numpy as np
from flint import arb, ctx
import model
import ladder
import screen_dense


class Transfers(unittest.TestCase):
    def setUp(self):
        ctx.prec = 256

    def test_baseline_exact_regression(self):
        current, old = model.Engine(128,19), ladder.Engine(16)
        z = F(-5)
        self.assertEqual(current.epoch(z),old.epoch(z))
        self.assertEqual(current.region(z,2),old.region(z,2))
        # Q1 returns balls rather than point upper bounds. Arb's == asks for
        # equality of exact numbers, so even identical nonzero-radius balls
        # need not compare equal. Compare midpoints and radii separately.
        for x,y in zip(current.q1(z),old.q1(z)):
            self.assertEqual(x.mid(),y.mid())
            self.assertEqual(x.rad(),y.rad())
        self.assertEqual(current.bernoulli(arb(1)/7,arb(1)/100),
                         old.bernoulli(arb(1)/7,arb(1)/100))

    def test_generic_float_independent_formulas(self):
        for t,s in ((16,10),(32,15),(64,15)):
            e = model.Engine(t,s)
            z = F(-4)
            actual = np.array([[float(x) for x in a] for a in e.epoch(z)]).reshape(t+1,e.n,e.n)
            logs = model.base.independent.g.epochs(e.spectrum,e.kernel,e.columns,e.caps,e.low,
                                                   math.exp(float(z)),maximum=t)
            np.testing.assert_allclose(actual,np.exp(logs),rtol=3e-13,atol=1e-300)
            actual = np.array([float(x) for x in e.bernoulli(arb(1)/7,arb(1)/100)]).reshape(e.n,e.n)
            expected = np.exp(screen_dense.bernoulli(e.spectrum,e.b_spectrum,e.kernel,1/7,.01))
            np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=1e-300)
            # Uniformly mixed fixed-j inputs have the exact same zero row.
            epoch = e.epoch(z)
            theta, lam = arb(1)/7, model.number(z).exp()
            weights = [math.comb(t,j)*theta**j*(1-theta)**(t-j) for j in range(t+1)]
            direct = e.bernoulli(theta,lam)
            for k in (0,1):
                value = sum((weights[j]*epoch[j][k] for j in range(t+1)),arb(0))
                self.assertLess(abs(float(value-direct[k])),1e-60)

    def test_no_frozen_global_mutation(self):
        self.assertEqual(model.base.independent.T,128)
        self.assertEqual(model.base.independent.base.T,128)
        e = model.Engine(32,15)
        self.assertEqual(e.length//e.t,16)
        self.assertEqual(e.identity()['inner']['t'],32)
        self.assertEqual(model.base.independent.T,128)
        c = model.Checker(32,15)
        self.assertEqual(1 << c.power,e.output_bits//32)
        self.assertEqual(c.fixed_input_bound.__globals__['bernoulli_activation'].activated,
                         model.activated)
        self.assertIsNot(model.fixed_input.bernoulli_activation.activated,model.activated)

    def test_activation_regression(self):
        import refine
        import activation_density
        current,old=refine.Engine(128,19),activation_density.Engine(16)
        self.assertEqual(current.epoch(F(-3)),old.epoch(F(-3)))


if __name__ == '__main__':
    unittest.main()
