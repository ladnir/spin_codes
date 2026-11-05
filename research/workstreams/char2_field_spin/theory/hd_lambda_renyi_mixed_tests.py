"""Small structural checks for the mixed comparison; no provider runs."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from flint import arb, ctx

import hd_lambda_renyi_mixed as target


class MixedRenyiTests(unittest.TestCase):
    def test_reviewed_dependencies_merge(self):
        pins = target.dependencies()
        self.assertEqual(pins['hd_lambda_renyi_orbits.py'], target.PINS['hd_lambda_renyi_orbits.py'])

    def test_uniform_moments_count_every_original_rowset(self):
        moments = {c: {3: F(1), 4: F(1)} for b in range(5)
                   for c in combinations(range(16), b)}
        with ctx.workprec(256):
            for p in target.ORDERS:
                weights = target.scheduled_weights(moments, p)
                for a, weight in zip(range(12, 17), weights):
                    self.assertTrue(weight.contains(comb(16, a)))

    def test_fixed_source_column_products_not_resampled_sources(self):
        # I3(c)=t(c)^3 gives d_A=product_j t(P_j A^c), an integer.
        scales = {c: 1+sum(c) % 2 for b in range(5) for c in combinations(range(16), b)}
        moments = {c: {3: F(s**3), 4: F(1)} for c, s in scales.items()}
        with ctx.workprec(256):
            weights = target.scheduled_weights(moments, 3)
            for a, observed in zip(range(12, 17), weights):
                expected = 0
                for c in combinations(range(16), 16-a):
                    term = 1
                    for j in range(16):
                        before = target.binding.column_permutations('input-inverse-output-xor', j)[0]
                        term *= scales[tuple(sorted(before[i] for i in c))]
                    expected += term
                self.assertTrue(observed.contains(expected))

    def test_positive_polynomial_includes_ordered_cohorts_once(self):
        with ctx.workprec(256):
            raw = [2, 3, 5, 7, 11]
            coeff = target.total_coefficients(list(map(arb, raw)), groups=3)
            expected = [0]*13
            for choices in product(range(5), repeat=3):
                term = 1
                for c in choices:
                    term *= raw[c]
                expected[sum(choices)] += term
            self.assertEqual(len(coeff), len(expected))
            for i, value in enumerate(expected):
                self.assertTrue(coeff[i].contains(value))

    def test_prescribed_scope_guards(self):
        with self.assertRaises(ValueError):
            target.scheduled_weights({}, 2)
        for weights, groups in [([arb(1)]*4, 18), ([arb(1)]*5, 0),
                                ([arb(0)]*5, 18), ([arb(1)]*5, True)]:
            with self.assertRaises(ValueError):
                target.total_coefficients(weights, groups)
        with self.assertRaises(ValueError):
            target.checked_inventory({})


if __name__ == '__main__':
    unittest.main()
