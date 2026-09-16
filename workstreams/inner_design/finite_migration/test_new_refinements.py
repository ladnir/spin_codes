"""Exact toy checks for the joint-map and constrained-density refinements."""
from fractions import Fraction as F
from itertools import product
import math
import unittest

from flint import ctx
import constant_density as density
import joint_overlap as joint


class JointAuditTests(unittest.TestCase):
    def test_exact_audit(self):
        result = joint.audit()
        self.assertEqual(result['rank'], 38)
        self.assertEqual(result['minimum_distance_lower'], 7)
        self.assertEqual(result['checked_syndromes'], sum(math.comb(128, i) for i in range(4)))

    def test_reject_invalid_hypotheses(self):
        with self.assertRaisesRegex(AssertionError, 'intersect'):
            joint.distance_audit([3], [3], 2, 1)
        with self.assertRaisesRegex(AssertionError, 'All-one'):
            joint.distance_audit([3], [4], 3, 1)
        with self.assertRaisesRegex(AssertionError, 'Low-weight'):
            joint.distance_audit([1], [7], 3, 1)

    def test_overlap_inequalities(self):
        # Exhaustive supports, excluding sum weights below d or above n-d.
        length, distance = 7, 2
        for a, b in product(range(1 << length), repeat=2):
            if distance <= (a ^ b).bit_count() <= length-distance:
                lo, hi = joint.overlap_ends(a.bit_count(), b.bit_count(), length, distance)
                self.assertLessEqual(lo, (a & b).bit_count())
                self.assertLessEqual((a & b).bit_count(), hi)


class ConstantDensityTests(unittest.TestCase):
    def test_vertices_lower_bound_every_toy_composition(self):
        ps = (F(1, 4), F(1, 2), F(3, 4))
        for probabilities in product(ps+(F(1),), repeat=4):
            alpha = F(probabilities.count(F(1)), 4)
            if alpha > F(1, 2):
                continue
            mean = sum(probabilities)/4
            if mean == 1:
                continue
            candidates = density.variance_vertices(ps, mean, mean, F(1, 2))
            for x in (F(1, 4), F(1), F(4)):
                values = {p:density.base.variance_term(p, x) for p in ps}
                bound = min(a*values[p]+b*values[r] for p,r,a,b in candidates)
                actual = sum(density.base.variance_term(p, x) for p in probabilities)/4
                self.assertLessEqual(bound, actual)

    def test_exact_poisson_binomial_ratio(self):
        ctx.prec = 192
        ps = (F(1, 4), F(1, 2), F(3, 4))
        rows, q = 8, 4
        checked = set()
        for probabilities in product(ps+(F(1),), repeat=q):
            probabilities = tuple(sorted(probabilities))
            if probabilities in checked or probabilities.count(F(1)) > 1:
                continue
            checked.add(probabilities)
            mean = sum(probabilities)/q
            upper = density.factor(rows, q, q, mean, mean, ps, F(1, 4))
            distribution = [F(1)]
            for p in probabilities:
                updated = [F(0)]*(len(distribution)+1)
                for j, value in enumerate(distribution):
                    updated[j] += value*(1-p)
                    updated[j+1] += value*p
                distribution = updated
            theta = F(q, rows)*mean
            for j, value in enumerate(distribution):
                denominator = math.comb(rows, j)*theta**j*(1-theta)**(rows-j)
                self.assertTrue(density.base.number(value/denominator) <= upper)

    def test_interval_ratio(self):
        ps = (F(1, 4), F(1, 2), F(3, 4))
        upper = density.factor(8, 3, 4, F(1, 3), F(3, 4), ps, F(1, 3))
        for probabilities in ((F(1, 4), F(1, 2), F(1)),
                              (F(1, 4), F(1, 4), F(1, 2), F(3, 4))):
            distribution = [F(1)]
            for p in probabilities:
                updated = [F(0)]*(len(distribution)+1)
                for j, value in enumerate(distribution):
                    updated[j] += value*(1-p)
                    updated[j+1] += value*p
                distribution = updated
            theta = sum(probabilities)/8
            for j, value in enumerate(distribution):
                denominator = math.comb(8, j)*theta**j*(1-theta)**(8-j)
                self.assertTrue(density.base.number(value/denominator) <= upper)

    def test_dense_point_matches_point_producer(self):
        import constant_dense
        import constant_point
        import retune_fixed_input
        from flint import arb
        ctx.prec = 256
        c = constant_dense.Checker(16)
        q, coordinate = 206, F(49, 256)
        witness = retune_fixed_input.proposal(c, q, coordinate, -14, -0.31)
        old = c.fixed_input_bound(q, q, coordinate, coordinate, witness)
        nu = c.density(coordinate, coordinate)[0]
        ratio = density.factor(c.rows, q, q, nu, nu, c.ps, constant_dense.MAXIMUM)
        original = c.comparison(q, q, coordinate, coordinate, ctx.prec)
        ordinary = old+256*(ratio.log()-original.log())
        exceptional = constant_point.evaluate(c, q, coordinate, witness,
                                               constant_dense.MAXIMUM, constant_dense.XI)
        point = (ordinary.exp()+exceptional.exp()).log()
        box = c.split_bound(q, q, coordinate, coordinate, witness)
        self.assertTrue(abs(point-box) < arb('1e-50'))


if __name__ == '__main__':
    unittest.main()
