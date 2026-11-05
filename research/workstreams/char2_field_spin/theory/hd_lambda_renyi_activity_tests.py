"""Small exact checks for the activity-tail comparison; no provider runs."""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb
import unittest
from flint import ctx
import hd_lambda_renyi_activity as target


class ActivityTests(unittest.TestCase):
    def test_pins(self):
        self.assertEqual(target.dependencies()['hd_lambda_renyi_middle.py'], target.PINS['hd_lambda_renyi_middle.py'])

    def test_tail_matches_labelled_bernoulli_expansion(self):
        z, atoms = F(2, 3), (F(1, 7), F(3, 5), F(4, 3), F(7, 2))
        for floor in range(5):
            exact = F(0)
            for bits in product((0, 1), repeat=4):
                if sum(bits) < floor:
                    continue
                term = F(1)
                for bit, nonzero in zip(bits, atoms):
                    term *= nonzero if bit else z
                exact += term
            with ctx.workprec(256):
                self.assertTrue(target.activity_tail(z, atoms, floor).contains(target.A.uniform.rational(exact)))

    def test_exact_cancellation_preserves_small_nonzero_atom(self):
        z = F(2)**2000
        i = z+F(1, 255)
        self.assertEqual(i-z, F(1, 255))
        with ctx.workprec(256):
            self.assertTrue(target.activity_tail(z, (i-z,)*2, 2).contains(target.A.uniform.rational(F(1, 255**2))))

    def test_zero_floor_recovers_product_of_moments(self):
        z, nonzero = F(2, 5), (F(1, 3), F(5, 7), F(11, 13))
        exact = F(1)
        for n in nonzero:
            exact *= z+n
        with ctx.workprec(256):
            self.assertTrue(target.activity_tail(z, nonzero, 0).contains(target.A.uniform.rational(exact)))

    def test_fixed_original_rowsets_and_multiplicity(self):
        # Uniform law at full input dimension: activity truncation can lower its norm.
        moments = {c: {p: F(1) for p in target.ORDERS}
                   for b in range(3) for c in combinations(range(16), b)}
        with ctx.workprec(256):
            weights, diag = target.scheduled_weights(moments, 3, minimum=0, arities=(14, 15, 16))
            for a, w, d in zip((14, 15, 16), weights, diag):
                self.assertTrue(w.contains(comb(16, a)))
                self.assertEqual(d['labelled_rowsets'], comb(16, a))
        with ctx.workprec(2048):
            truncated, _ = target.scheduled_weights(moments, 3, arities=(16,))
            self.assertTrue(truncated[0] <= 1)

    def test_guards(self):
        with self.assertRaises(ValueError):
            target.checked_middle({})
        for z, atoms, minimum in [(F(-1), (F(1),), 1), (F(1), (F(-1),), 1), (F(1), (F(1),), 2)]:
            with self.assertRaises(ValueError):
                target.activity_tail(z, atoms, minimum)


if __name__ == '__main__':
    unittest.main()
