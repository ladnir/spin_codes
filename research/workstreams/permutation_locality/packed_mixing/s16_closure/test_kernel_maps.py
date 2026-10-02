"""Exhaustive tiny-map checks of uniform GL and explicit-map adapters."""
from fractions import Fraction as Q
from math import comb
from types import SimpleNamespace
import unittest

import numpy as np
from flint import arb, ctx
import kernel_maps as kernel
import scalar_cover as sc


def endpoint(value):
    mantissa, exponent = value.upper().man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


class KernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.images = tuple(sum(((s >> j)&1)*v for j, v in enumerate((1, 6, 120))) for s in range(8))
        cls.columns = [1, 2, 4, 3, 5, 7, 6, 1]
        cls.uniform = kernel.prepare_maps(cls.images, cls.columns, 3, distribution='uniform_gl', birth_density='classes')
        cls.transvections = kernel.prepare_maps(cls.images, cls.columns, 3, 4, 'transvections')

    def setUp(self):
        ctx.prec = 192

    def exact(self, activity, z):
        result = [[Q(0)]*8 for _ in range(8)]
        for state in range(8):
            for x in range(256):
                probability = Q(1)
                for offset in (0, 4):
                    probability *= activity/15 if (x >> offset)&15 else 1-activity
                feedback = 0
                for bit, column in enumerate(self.columns):
                    if x >> bit & 1:
                        feedback ^= column
                weight = probability*z**(self.images[state]^x).bit_count()
                if state:
                    for refreshed in range(1, 8):
                        result[state][refreshed^feedback] += weight/7
                else:
                    result[state][feedback] += weight
        return result

    def test_uniform_has_no_update_parameter(self):
        self.assertNotIn('updates', self.uniform)
        self.assertEqual(self.transvections['updates'], 4)
        changed = dict(self.uniform, columns=[2, *self.columns[1:]])
        with self.assertRaises(ArithmeticError):
            kernel.authenticate(changed)
        with self.assertRaises(ValueError):
            kernel.authenticate(dict(self.uniform, updates=32))

    def test_every_state_and_eight_steps(self):
        data = self.uniform
        for activity in (Q(0), Q(1, 5), Q(15, 16), Q(1)):
            matrix = kernel.outward_at_z(data, sc.probabilities(activity), arb(3)/4)
            exact = self.exact(activity, Q(3, 4))
            n = matrix.nrows()
            levels = list(map(int, data['birth_class_levels']))
            self.assertGreaterEqual(endpoint(matrix[0, 0]), exact[0][0])
            for target, level in enumerate(levels, 3):
                self.assertGreaterEqual(endpoint(matrix[0, target]),
                    sum(exact[0][s] for s in range(1, 8) if self.images[s].bit_count() == level))
            for state in range(1, 8):
                source = 3+levels.index(self.images[state].bit_count())
                self.assertGreaterEqual(endpoint(matrix[source, 0]), exact[state][0])
                for target in range(1, 8):
                    self.assertGreaterEqual(endpoint(matrix[source, 2])/7, exact[state][target])
            values = [Q(1)]*8
            for steps in range(1, 9):
                values = [sum(a*b for a, b in zip(row, values)) for row in exact]
                upper = matrix**steps
                sums = [endpoint(sum((upper[i, j] for j in range(n)), arb(0))) for i in range(n)]
                self.assertGreaterEqual(sums[0], values[0])
                self.assertGreaterEqual(sums[2], sum(values[1:])/7)
                for state in range(1, 8):
                    self.assertGreaterEqual(sums[1], values[state])
                    self.assertGreaterEqual(sums[3+levels.index(self.images[state].bit_count())], values[state])

    def test_float_matches_outward(self):
        for activity in (Q(1, 100), Q(1, 2), Q(15, 16), Q(1)):
            probabilities = sc.probabilities(activity)
            for tilt in (Q(1, 100), Q(1, 4), Q(1)):
                floating = kernel.floating(self.uniform, np.array(list(map(float, probabilities))), float(tilt))
                exact = kernel.outward(self.uniform, probabilities, tilt)
                np.testing.assert_allclose(floating, np.array([[float(exact[i, j]) for j in range(exact.ncols())]
                    for i in range(exact.nrows())]), rtol=2e-12, atol=1e-15)

    def test_uniform_iid_matches_conditional_exact_entries(self):
        import sparse_kernel
        import adapters_regional as adapters
        activity, z = Q(1, 5), arb(3)/4
        local = sparse_kernel.outward_at_z(self.uniform, z)
        n = local[0].nrows()
        mixture = sum((kernel.aq(comb(2, j)*activity**j*(1-activity)**(2-j))*m
                       for j, m in enumerate(local)), arb(0)*local[0])
        iid = kernel.outward_at_z(self.uniform, sc.probabilities(activity), z)
        for target in range(n):
            self.assertLess(abs(float(iid[0, target]-mixture[0, target])), 1e-45)
        self.assertLess(abs(float(iid[2, 2]-mixture[2, 2])), 1e-45)
        via_adapter = adapters.local_operators(SimpleNamespace(data=self.uniform),
            dict(parameters=['1/4'], regional_exact_zero=True))
        self.assertEqual(len(via_adapter), 3)

    def test_map_specific_refinements_do_not_call_actual(self):
        import adapters_regional as adapters
        import return_moment
        import lazy_density
        def forbidden(*args, **kwargs):
            raise AssertionError('legacy actual-map adapter was called')
        saved_joint, saved_lazy = return_moment.actual_census, lazy_density.actual_caps
        try:
            return_moment.actual_census = lazy_density.actual_caps = forbidden
            model = SimpleNamespace(data=self.transvections)
            witness = dict(parameters=['1/4'], regional_exact_zero=True,
                regional_joint_return_through=2, regional_lazy_density_through=2,
                regional_feedback_classes_from=1, regional_feedback_classes_through=2,
                regional_feedback_uniform_classes=True, regional_feedback_uniform_replace=True)
            local = adapters.local_operators(model, witness)
            self.assertEqual(len(local), 3)
            self.assertTrue(all(m.nrows() == 3+len(self.transvections['birth_class_levels']) for m in local))
        finally:
            return_moment.actual_census, lazy_density.actual_caps = saved_joint, saved_lazy


if __name__ == '__main__':
    unittest.main()
