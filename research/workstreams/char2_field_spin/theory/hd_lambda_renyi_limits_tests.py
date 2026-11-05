from fractions import Fraction as F
from itertools import product
from math import comb
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_renyi_limits as limits
import hd_lambda_renyi_support as support


class RenyiLimitsTests(unittest.TestCase):
    def test_distance_table_and_exact_full_space(self):
        self.assertEqual([limits.minimum_distance_bound(a) for a in range(1, 17)],
                         [16, 6, 5, 4, 3]+[2]*10+[1])
        for q, p in product((4, 256), (2, 3, 4)):
            expected = F(1, q**16), 1-F(1, q**16)
            self.assertEqual(limits.column_atoms(16, p, limits.COARSE, q=q), expected)
            self.assertEqual(limits.column_atoms(16, p, limits.IDEAL, q=q), expected)

    def test_singleton_coarse_is_exact(self):
        for q, p in product((4, 256), (2, 3, 4)):
            z, n = limits.column_atoms(1, p, limits.COARSE, q=q)
            self.assertEqual(z, F(q)**(15*p-16))
            self.assertEqual(n, F(q)**(15*p-16)*F(q-1)**(16-15*p))

    def test_toy_exact_moments_between_ideal_and_coarse(self):
        for rows in ((0b1011,), (0b1011, 0b0101), (0b10101, 0b01110)):
            width = max(rows).bit_length()
            data = support.support_counts_from_binary(rows, width, q=4)
            a = width-data['codimension']
            distance = min(s.bit_count() for s, n in enumerate(data['counts']) if s and n)
            for p in (2, 3, 4):
                z = F(4)**(p*(width-a)-width)
                actual = support.density_moment_from_counts(data['counts'], data['codimension'], q=4, order=p)-z
                ideal = (1-F(4)**(-a))**p/(1-F(4)**(-width))**(p-1)
                coarse = limits.coarse_nonzero_moment(a, p, distance, q=4, width=width)
                self.assertLessEqual(ideal, actual)
                self.assertLessEqual(actual, coarse)

    def test_positive_tail_and_raw_normalization(self):
        with ctx.workprec(256):
            z, n = F(5, 7), F(3, 11)
            exact = sum((F(comb(4, k))*z**(4-k)*n**k for k in range(2, 5)), F(0))
            self.assertTrue(limits.activity.activity_tail(z, (n,)*4, 2).contains(limits.A.uniform.rational(exact)))
            for kind in limits.KINDS:
                weights, _ = limits.raw_group_weights(3, kind)
                for a in (1, 8, 16):
                    z, n = limits.column_atoms(a, 3, kind)
                    tail = sum((F(comb(16, k))*z**(16-k)*n**k for k in range(9, 17)), F(0))
                    exact_cube = F(comb(16, a))**3*F(256)**(48*a)/F(255)**(24*a)*tail
                    self.assertTrue((weights[a-1]**3).contains(limits.A.uniform.rational(exact_cube)))

    def test_raw_coefficient_receives_no_extra_beta(self):
        with ctx.workprec(256):
            for p in (2, 3, 4):
                got = limits.total_log_bound(arb(255)*2**15, p, arb(0), arb(0))
                expected = limits.A.uniform.rational(limits.ALPHA*F(p, p-1)*15)
                self.assertTrue(got.contains(expected))

    def test_mock_receipt_scope_one_provider(self):
        sentinel = object()
        terms = [(i, 0, 1, arb(1)) for i in range(99)]
        fake = dict(controls=[dict(kind=limits.COARSE), dict(kind=limits.IDEAL)], selected_totals=[])
        with patch.object(limits.A.uniform, 'regional_directed', return_value=sentinel) as provider, \
                patch.object(limits.A, 'source_profile_terms', return_value=terms) as source, \
                patch.object(limits, 'controls_from_uniform_log', return_value=fake):
            data = limits.receipt(256)
        provider.assert_called_once_with(limits.MARK, limits.ALPHA, 18, 60)
        source.assert_called_once_with(18, sentinel)
        self.assertEqual(data['total_arity_range'], [18, 245])
        self.assertTrue(data['coarse_is_failure_upper_bound'])
        self.assertFalse(data['ideal_is_failure_bound'])
        self.assertFalse(data['whole_code_certificate'])
        self.assertEqual(data['fresh_regional_banks'], 1)

    def test_guards_and_pins(self):
        self.assertEqual(limits.dependencies()['hd_lambda_renyi_activity.py'], limits.PINS['hd_lambda_renyi_activity.py'])
        for a, p, kind in ((0, 3, limits.COARSE), (17, 3, limits.IDEAL), (1, 1, limits.COARSE), (1, 3, 'ideal-upper')):
            with self.assertRaises(ValueError):
                limits.column_atoms(a, p, kind)
        with self.assertRaises(ValueError):
            limits.coarse_nonzero_moment(4, 3, 2, q=4, width=4)
        with self.assertRaises(ValueError):
            limits.receipt(128)


if __name__ == '__main__':
    unittest.main()
