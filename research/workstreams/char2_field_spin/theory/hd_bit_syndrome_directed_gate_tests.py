"""Small independent controls for final aggregation and ordered gap words."""
from fractions import Fraction as F
from unittest import TestCase, main

import numpy as np
from flint import arb, ctx

import hd_positive_upward as up
from hd_bit_syndrome_atoms_directed import export_scaled, rational, upper_fraction
from hd_bit_syndrome_directed_gate import (aggregate, covered_bounds, ordered_products,
    scaled_fraction, signed_float_upper, MARK, LIMIT)


def exact(value):
    return np.array([scaled_fraction(x, value.exponent) for x in value.mantissas.flat],
                    dtype=object).reshape(value.mantissas.shape)


class DirectedGateTests(TestCase):
    def test_root_sum_perfect_powers(self):
        inputs = up.Scaled(np.array([2.**18, float(3**18), 0.]), 0)
        answer = upper_fraction(aggregate(inputs, precision=256))
        target = F(5**18)
        self.assertGreaterEqual(answer, target)
        self.assertLess(answer-target, target/F(2**200))

    def test_group_one_and_huge_scale(self):
        inputs = up.Scaled(np.array([.5, .25, .125]), 5000)
        answer = upper_fraction(aggregate(inputs, groups=1, precision=384))
        target = F(7, 8)*2**5000
        self.assertGreaterEqual(answer, target)
        self.assertLess(answer-target, target/F(2**300))

    def test_ordered_gap_words(self):
        matrices = up.Scaled(np.array([[[1., .5], [.25, 1.]], [[.5, 1.], [1., .25]]]), 0)
        parents, labels = (-1, 0, 0, 1, 2), (-1, 0, 1, 1, 0)
        answer = exact(ordered_products(matrices, parents, labels))
        a, b = exact(matrices)
        expected = [np.eye(2, dtype=object), a, b, a@b, b@a]
        self.assertTrue(np.all(answer >= np.asarray(expected, dtype=object)))
        self.assertFalse(np.array_equal(a@b, b@a))

    def test_prefix_and_signed_log_export(self):
        counts = dict(F0=2, F1=3, F_le1=5)
        with ctx.workprec(256):
            result = covered_bounds(arb(2)**1100, counts)
        with ctx.workprec(512):
            for key, count in counts.items():
                independent = arb(count)*arb(2)**1100 * (arb(1000)/993)**LIMIT / arb(255)**1153
                upper = rational(F(result[key]['bound_upper_dyadic']))
                self.assertTrue(upper >= independent)
                log_target = independent.log()/arb(2).log()
                self.assertTrue(arb(result[key]['log2_upper']) >= log_target)
                self.assertTrue(result[key]['below_2_minus_40'])
            negative = -arb(1)/3
            self.assertTrue(arb(signed_float_upper(negative)) >= negative)

    def test_export_import_retains_endpoint(self):
        with ctx.workprec(384):
            values = np.array([arb(0), (arb(1)/3).upper(), (arb(2)**-5000).upper()], dtype=object)
            wanted = np.array([upper_fraction(x) for x in values], dtype=object)
        with ctx.workprec(32):
            received = up.from_export(export_scaled(values))
        self.assertTrue(np.all(exact(received) >= wanted))
        self.assertEqual(exact(received)[0], 0)
        self.assertGreater(exact(received)[2], 0)


if __name__ == '__main__':
    main()
