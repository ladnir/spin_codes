"""Tiny local/coefficient and mocked receipt checks; no K20 distance gate."""

from fractions import Fraction as F
import unittest
from unittest.mock import patch

from flint import arb, ctx
import numpy as np

import hd_field128_full_rs_gamma_certificate as cert


base, local = cert.base, cert.local


def log_float(value):
    return -float('inf') if value.is_zero() else float(value.log()/arb(2).log())


class GammaDirectedTests(unittest.TestCase):
    def test_exact_upward_packets_at_tiny_mark(self):
        z = F(1, 1<<200)
        with ctx.workprec(256):
            for state in (24, 30):
                table = cert.packet_table(z, 2, state=state)
                for j, values in table.items():
                    for actual, exact in zip(values, local.exact_packet(j, z, state=state)):
                        self.assertTrue(actual >= base.rational(exact))
                self.assertTrue(all(value > 0 for value in table[1]))
                self.assertLess(log_float(table[1][2]), -128*state)

    def test_directed_coefficients_match_gamma_float_and_raw_order(self):
        z, alpha, state = F(1, 2), F(1, 4), 24
        floating = local.AliasCoefficients(float(z), float(alpha), state=state)
        with ctx.workprec(256):
            singles, pairs = cert.directed_parts(z, alpha, 2, state, 80, base.Budget(10))
            for (s, f), values in singles.items():
                np.testing.assert_allclose([log_float(v) for v in values], floating.single[:, s, f], atol=1e-10, rtol=0)
            for (s, f), values in pairs.items():
                np.testing.assert_allclose([log_float(v) for v in values], floating.composed[:, s, f], atol=1e-10, rtol=0)
            table = cert.packet_table(z, 2, state=state)
            expected = [arb(0) for unused in range(4)]
            for left, right in ((table[0], table[1]), (table[1], table[0])):
                for entry, value in enumerate(base.matrix_product(left, right)):
                    if not value.is_zero():
                        expected[entry] += 32*value**base.rational(alpha)
            for actual, oracle in zip(pairs[1, 0], expected):
                self.assertTrue(actual.overlaps(oracle))

    def test_region_alpha1_and_reusable_prefix_marker_cancellation(self):
        with ctx.workprec(256):
            for state in (24, 30):
                paired = cert.regional_directed(F(7, 8), F(1), 2, state=state)
                unpaired = cert.regional_directed(F(7, 8), F(1), 2, state=state, paired=False)
                self.assertEqual(str(paired), str(unpaired))
            small = cert.regional_directed(F(1, 2), F(1, 4), 1)
            wide = cert.regional_directed(F(1, 2), F(1, 4), 2)
            for key in small:
                np.testing.assert_allclose([log_float(v) for v in small[key]],
                                           [log_float(v) for v in wide[key]], atol=1e-10, rtol=0)

    def test_mocked_receipt_metadata_ceiling_and_closure(self):
        # No real regional coefficients or selected-weight gate in this test.
        artificial = {(1, 0): tuple(map(arb, (1, 0, 0, 1)))}
        with patch.object(cert, 'regional_directed', return_value=artificial):
            row = cert.weight_receipt(1, F(1, 2), F(1, 4), state=30)
        self.assertEqual(row['construction'], local.construction(30))
        self.assertEqual(row['feedback'], 'gamma A transpose')
        self.assertEqual(row['gamma_extension_degree'], 16)
        self.assertEqual(row['profiles'], 1)
        self.assertTrue(row['covers_all_supports_at_weight'])
        self.assertTrue(row['directed'])
        self.assertFalse(row['whole_code_certificate'])
        self.assertEqual(row['projective_divisor'], str((1<<128)-1))
        with ctx.workprec(256):
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())
        for name in ('hd_field128_full_rs_gamma_certificate.py', 'hd_field128_full_rs_gamma_screen.py',
                     'hd_field128_full_rs_state_screen.py', 'hd_field128_lift_certificate.py'):
            self.assertIn(name, row['source_sha256'])

    def test_scope_guards_and_source_conflict(self):
        for kwargs in (dict(state=32), dict(bits=8), dict(paired=1)):
            with self.assertRaises(ValueError):
                cert.regional_directed(F(1, 2), F(1, 4), 2, **kwargs)
        with self.assertRaises(ValueError):
            cert.regional_directed(.5, F(1, 4), 2)
        with self.assertRaises(ValueError):
            cert.regional_directed(F(1, 2), F(1, 4), 65)
        with self.assertRaises(TimeoutError):
            cert.packet_table(F(1, 2), 1, budget=base.Budget(1e-12))
        artificial = {(1, 0): tuple(map(arb, (1, 0, 0, 1)))}
        with patch.object(cert, 'regional_directed', return_value=artificial):
            with self.assertRaises(RuntimeError):
                cert.weight_receipt(1, F(1, 2), F(1, 4), _pins={})


if __name__ == '__main__':
    unittest.main()
