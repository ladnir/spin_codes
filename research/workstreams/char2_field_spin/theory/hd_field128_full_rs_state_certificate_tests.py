"""Tiny exact/directed full-RS state checks; no construction-size gate."""

from fractions import Fraction as F
from math import comb
import unittest

from flint import arb, ctx
import numpy as np

import hd_field128_full_rs_state_certificate as certificate
import hd_field128_full_rs_screen as fixed16
from hd_field128_full_rs_screen_tests import toy_moments


base, local = certificate.base, certificate.local


def log_float(value):
    return -float('inf') if value.is_zero() else float(value.log() / arb(2).log())


class DirectedStateTests(unittest.TestCase):
    def test_exact_local_oracle_toy_and_state16_control(self):
        a = local.lagrange(range(8), range(8, 11), 16)
        for support in ((), (0,), (0, 1), (1, 6)):
            quiet, actual_return = toy_moments(a, support, q=16)
            k, birth, ret, live = local.exact_moments(len(support), F(1, 2), 16, 8, 3)
            self.assertEqual(k, quiet)
            self.assertGreaterEqual(ret, actual_return)
        for joint in (False, True):
            for occupancy in (0, 1, 16, 33, 64):
                self.assertEqual(local.exact_packet(occupancy, F(7, 8), 256, state=16, joint=joint),
                                 fixed16.exact_packet(occupancy, F(7, 8), 256, joint=joint))

    def test_table_rounds_up_and_retains_tiny_fullfield_values(self):
        z, q = F(1, 1 << 200), 1 << 128
        with ctx.workprec(256):
            for state in (24, 32):
                table = certificate.packet_table(z, 1, q=q, state=state, budget=base.Budget(10))
                for occupancy, actual in table.items():
                    exact = local.exact_packet(occupancy, z, q, state=state)
                    for stored, value in zip(actual, exact):
                        self.assertTrue(stored >= base.rational(value))
                self.assertTrue(all(value > 0 for value in table[1]))
                self.assertLess(log_float(table[1][2]), -128*state)

    def test_swap_geometry_and_gg_before_alpha(self):
        z, alpha, maximum, q, state = F(7, 8), F(1, 4), 2, 256, 24
        for s, f, count in local.fixed16.geometry(maximum):
            self.assertEqual(count, comb(32, s) * comb(32 - s, f))
        with ctx.workprec(256):
            single, composed = certificate.directed_parts(z, alpha, maximum, q, state, 80, base.Budget(10))
            table = certificate.packet_table(z, maximum, q=q, state=state, budget=base.Budget(10))
            # Degree(1,0) consists of an ordered empty/one-singleton pair
            # in either order, each with multiplicity32 outside alpha.
            expected = [arb(0) for _ in range(4)]
            for left, right in ((table[0], table[1]), (table[1], table[0])):
                raw = base.matrix_product(left, right)
                for entry, value in enumerate(raw):
                    expected[entry] += 32 * value**base.rational(alpha) if not value.is_zero() else 0
            for actual, oracle in zip(composed[1, 0], expected):
                self.assertTrue(actual.overlaps(oracle))
            float_single = local.single_coefficients(float(z), float(alpha), maximum, q, state, local.old.Budget(10))
            float_pair = local.good_pair_coefficients(float(z), float(alpha), maximum, q, state, local.old.Budget(10))
            for (s, f), values in single.items():
                np.testing.assert_allclose([log_float(x) for x in values], float_single[:, s, f], rtol=0, atol=1e-10)
            for (s, f), values in composed.items():
                np.testing.assert_allclose([log_float(x) for x in values], float_pair[:, s, f], rtol=0, atol=1e-10)

    def test_alpha1_reusable_region_and_fractional_float_prefix(self):
        with ctx.workprec(256):
            for state in local.STATES:
                for alpha in (F(1), F(1, 4)):
                    region = certificate.regional_directed(F(7, 8), alpha, 2, q=256, state=state,
                                                           budget=base.Budget(10))
                    control = local.regional_screen(7/8, float(alpha), 2, q=256, state=state, max_seconds=10)
                    for key, values in region.items():
                        np.testing.assert_allclose([log_float(x) for x in values], np.array(control[key]).ravel(),
                                                   rtol=0, atol=1e-9)
                    if alpha == 1:
                        unpaired = certificate.regional_directed(F(7, 8), alpha, 2, q=256, state=state,
                                                                 paired=False, budget=base.Budget(10))
                        self.assertEqual(str(region), str(unpaired))
            # Compatibility is with single-RS state16, never the old four-block h4 map.
            self.assertNotEqual(local.construction(16), base.CONSTRUCTION)

    def test_tiny_selected_receipt_all_profiles_scope_and_float_agreement(self):
        # q256/r1 is a tiny harness check, not a full-field distance gate.
        row = certificate.weight_receipt(1, F(7, 8), F(1, 4), q=256, state=24, distance=F(2, 5), max_seconds=20)
        control = local.weight_screen(1, 7/8, .25, q=256, state=24, distance=F(2, 5), max_seconds=10)
        self.assertLess(abs(float(F(row['log2_upper_rational'])) - control['log2_upper_screen']), 1e-7)
        self.assertEqual(row['profiles'], 1)
        self.assertEqual(row['construction'], local.construction(24))
        self.assertEqual(row['state_dimension'], 24)
        self.assertEqual(row['emission_minimum'], 41)
        self.assertEqual(row['projective_divisor'], '255')
        self.assertTrue(row['covers_all_supports_at_weight'])
        self.assertTrue(row['covers_all_even_odd_splits'])
        self.assertTrue(row['directed'])
        self.assertFalse(row['whole_code_certificate'])
        with ctx.workprec(256):
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())

    def test_frozen_polynomial_helpers_and_guards(self):
        with ctx.workprec(256):
            source = {(2, 0): tuple(map(arb, (1, 0, 0, 1)))}
            encoded = base.encode(source, 2, F(1))
            squared = base.polynomial_product(encoded, encoded, 2, base.Budget(10))
            self.assertTrue(all(poly[index].is_zero() for poly in squared for index in range(6)))
        for kwargs in (dict(state=20), dict(q=16), dict(bits=8), dict(paired=1)):
            with self.assertRaises(ValueError):
                certificate.regional_directed(F(7, 8), F(1, 4), 2, **kwargs)
        with self.assertRaises(ValueError):
            certificate.regional_directed(F(7, 8), F(1, 4), 65)
        with self.assertRaises(ValueError):
            certificate.regional_directed(.875, F(1, 4), 2)
        with self.assertRaises(TimeoutError):
            certificate.packet_table(F(7, 8), 1, budget=base.Budget(1e-12))
        pins = certificate.source_pins()
        for name in ('hd_field128_full_rs_state_certificate.py', 'hd_field128_full_rs_state_screen.py',
                     'hd_field128_full_rs_screen.py', 'hd_field128_lift_certificate.py', 'hd_shifted_pairs_count.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
