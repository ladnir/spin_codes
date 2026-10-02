from fractions import Fraction as Q
from unittest.mock import patch
import unittest

import numpy as np
from scipy.optimize import linprog
from flint import arb, ctx
import joint_count_dual as joint


class JointCountDualTests(unittest.TestCase):
    def setUp(self):
        old = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', old)
        ctx.prec = 256

    def assert_close_upper(self, upper, exact, tolerance=1e-10):
        self.assertGreaterEqual(upper, joint._arb(exact))
        self.assertLess(float(upper-joint._arb(exact)), tolerance)

    def test_atom_caps_and_normalization_improve_separate_sum(self):
        value, witness = joint.bound([1, 2, 3], [1, 1, 1])
        self.assert_close_upper(value, 3)
        self.assertEqual(value, joint.verify([1, 2, 3], [1, 1, 1], (), witness))

    def test_raw_mgf_and_mean_match_small_primal_lp(self):
        values = [Q(1), Q(2), Q(8), Q(3)]
        caps = [Q(4, 5), Q(3, 5), Q(2, 5), Q(1, 5)]
        raw = [(Q(1, 2), arb(2).log())]
        moments = [(list(range(4)), Q(6, 5)), ([-j for j in range(4)], Q(-4, 5))]
        primal = linprog(-np.array(list(map(float, values))),
            A_ub=[np.exp(np.arange(4)/2), np.arange(4), -np.arange(4)],
            b_ub=[2., 1.2, -.8], A_eq=[np.ones(4)], b_eq=[1.],
            bounds=[(0., float(cap)) for cap in caps], method='highs')
        self.assertTrue(primal.success)
        prepared = joint.prepare(caps, raw, moment_bounds=moments)
        upper, witness = prepared.bound(values)
        self.assertAlmostEqual(float(upper), -primal.fun, places=9)
        self.assertEqual(upper, prepared.verify(values, witness))

    def test_enumerated_feasible_distributions_are_bounded(self):
        values = [Q(1), Q(7), Q(2)]
        raw = [(Q(1), arb(3).log())]
        upper, _ = joint.bound(values, [Q(4, 5)]*3, raw)
        for a in range(9):
            for b in range(9):
                c = 10-a-b
                if not 0 <= c <= 8:
                    continue
                ps = [Q(a, 10), Q(b, 10), Q(c, 10)]
                if sum(joint._arb(p)*arb(j).exp() for j, p in enumerate(ps)) > 3:
                    continue
                self.assertGreaterEqual(upper, joint._arb(sum(p*f for p, f in zip(ps, values))))

    def test_tiny_atom_large_objective_never_divides_deficit_by_cap(self):
        tiny, huge = arb(-3000).exp(), arb(3000).exp()
        prepared = joint.prepare([1, tiny], [(Q(3000), arb(2).log())])
        upper, witness = prepared.bound([0, huge])
        self.assertTrue(upper.is_finite())
        self.assertLess(float(upper), 1.00000001)
        self.assertGreater(upper, .99999999)
        perturbed = dict(witness, alpha='1/1000000000')
        self.assertTrue(prepared.verify([0, huge], perturbed).is_finite())

    def test_zero_caps_and_zero_objective(self):
        upper, witness = joint.bound([10, 3], [0, 1])
        self.assert_close_upper(upper, 3)
        self.assertEqual(joint.bound([0, 0], [1, 1])[0], 0)
        self.assertEqual(joint.verify([10, 3], [0, 1], (), witness), upper)

    def test_rounded_infeasible_dual_is_repaired_coordinatewise(self):
        prepared = joint.prepare([1, 1, 1])
        witness = dict(schema=joint.SCHEMA, size=3, scale_exponent=3,
                       alpha='374999999999/1000000000000', beta=[])
        upper = prepared.verify([1, 2, 3], witness)
        self.assertGreaterEqual(upper, 3)
        self.assertLess(float(upper), 3.000000001)

    def test_bad_solver_or_precision_fails_safely(self):
        prepared = joint.prepare([1, 1])
        with patch.object(joint, 'linprog', return_value=type('Fit', (), {'success': False})()):
            upper, _ = prepared.bound([1, 2])
        self.assert_close_upper(upper, 3)
        ctx.prec = 384
        with self.assertRaises(ArithmeticError):
            prepared.bound([1, 2])

    def test_invalid_inputs_and_witnesses(self):
        for caps in ([], [-1, 2], [Q(1, 4)]*2, [arb('nan')], [True]):
            with self.assertRaises(ValueError):
                joint.prepare(caps)
        for mgfs in (({'tilted_atom': True},), ((True, 0),), ((Q(1), arb('nan')),)):
            with self.assertRaises(ValueError):
                joint.prepare([1, 1], mgfs)
        with self.assertRaises(ValueError):
            joint.prepare([1, 1], moment_bounds=[([0], 1)])
        prepared = joint.prepare([1, 1], [(Q(1), arb(3).log())])
        _, witness = prepared.bound([1, 2])
        for changes in ({'size': 3}, {'beta': ['-1']}, {'alpha': 1.0},
                        {'scale_exponent': True}, {'upper': 0}):
            with self.assertRaises(ValueError):
                prepared.verify([1, 2], dict(witness, **changes))
        with self.assertRaises(ValueError):
            prepared.bound([-1, 2])


if __name__ == '__main__':
    unittest.main()
