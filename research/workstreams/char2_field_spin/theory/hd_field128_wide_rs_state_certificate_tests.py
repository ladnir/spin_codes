"""Tiny wide directed checks; no construction-size numerical gate."""

from fractions import Fraction as F
from math import comb
import unittest

from flint import arb, ctx
import numpy as np

import hd_field128_wide_rs_state_certificate as certificate
import hd_field128_full_rs_state_certificate as narrow
from hd_field128_full_rs_screen_tests import toy_moments


base, local = certificate.base, certificate.local


def log_float(value):
    return -float('inf') if value.is_zero() else float(value.log() / arb(2).log())


class DirectedWideTests(unittest.TestCase):
    def test_exact_caps_toy_and_narrow_compatibility(self):
        a = local.frozen.lagrange(range(8), range(8, 11), 16)
        for support in ((), (0,), (0, 1), (1, 6)):
            quiet, actual_return = toy_moments(a, support, q=16)
            cap = local.exact_moments(len(support), F(1, 2), 16, 8, 3)
            self.assertEqual(cap[0], quiet)
            self.assertGreaterEqual(cap[2], actual_return)
        with ctx.workprec(256):
            for state in narrow.STATES:
                wide_table = certificate.packet_table(F(7, 8), 2, q=256, width=64, state=state)
                old_table = narrow.packet_table(F(7, 8), 2, q=256, state=state)
                self.assertEqual(str(wide_table), str(old_table))
                for alpha in (F(1), F(1, 4)):
                    wide = certificate.regional_directed(F(7, 8), alpha, 2, q=256, width=64, state=state)
                    old = narrow.regional_directed(F(7, 8), alpha, 2, q=256, state=state)
                    self.assertEqual(str(wide), str(old))

    def test_fullfield_upward_rounding_including_state48(self):
        z, q = F(1, 1 << 200), 1 << 128
        with ctx.workprec(256):
            for state in (1, 48, 64):
                table = certificate.packet_table(z, 1, q=q, width=128, state=state, budget=base.Budget(10))
                for occupancy, actual in table.items():
                    exact = local.exact_packet(occupancy, z, q, 128, state)
                    for stored, value in zip(actual, exact):
                        self.assertTrue(stored >= base.rational(value))
                self.assertTrue(all(value > 0 for value in table[1]))
                self.assertLess(log_float(table[1][2]), -128*state)

    def test_width_specific_integer_region_normalizer(self):
        maximum = 4
        for width in certificate.WIDTHS:
            route = local.routing(width)
            self.assertEqual(route['packets_per_region'] * width, 8192)
            self.assertEqual(route['composed_pairs_per_region'] * 2, route['packets_per_region'])
            counts = {(0, 0): 1}
            for unused in range(route['packets_per_region']):
                next_counts = {}
                for (s, f), count in counts.items():
                    for t, g, other in local.geometry(maximum, width=width):
                        self.assertEqual(other, comb(width//2, t) * comb(width//2-t, g))
                        if s+t+2*(f+g) <= maximum:
                            key = s+t, f+g
                            next_counts[key] = next_counts.get(key, 0) + count*other
                counts = next_counts
            self.assertEqual(counts, {(s, f): comb(4096, s)*comb(4096-s, f)
                                      for s, f in base.valid_degrees(maximum)})
        self.assertEqual(local.routing(128)['pair_slots_per_packet'], 64)
        self.assertEqual(local.routing(128)['packets_per_region'], 64)
        self.assertEqual(local.routing(128)['composed_pairs_per_region'], 32)

    def test_ordered_raw_gg_and_unpowered_multiplicities(self):
        z, alpha, maximum, q, width, state = F(7, 8), F(1, 4), 2, 256, 128, 48
        with ctx.workprec(256):
            single, composed = certificate.directed_parts(z, alpha, maximum, q, width, state, 80, base.Budget(10))
            table = certificate.packet_table(z, maximum, q=q, width=width, state=state)
            expected = [arb(0) for unused in range(4)]
            for left, right in ((table[0], table[1]), (table[1], table[0])):
                for entry, value in enumerate(base.matrix_product(left, right)):
                    if not value.is_zero():
                        expected[entry] += 64 * value**base.rational(alpha)
            for actual, oracle in zip(composed[1, 0], expected):
                self.assertTrue(actual.overlaps(oracle))
            float_single = local.single_coefficients(float(z), float(alpha), maximum, q, width, state, local.old.Budget(10))
            float_pair = local.good_pair_coefficients(float(z), float(alpha), maximum, q, width, state, local.old.Budget(10))
            for (s, f), values in single.items():
                np.testing.assert_allclose([log_float(x) for x in values], float_single[:, s, f], rtol=0, atol=1e-10)
            for (s, f), values in composed.items():
                np.testing.assert_allclose([log_float(x) for x in values], float_pair[:, s, f], rtol=0, atol=1e-10)

    def test_wide_regional_float_and_alpha1_replay(self):
        with ctx.workprec(256):
            for state in (1, 48, 64):
                for alpha in (F(1), F(1, 4)):
                    region = certificate.regional_directed(F(7, 8), alpha, 2, q=256, width=128, state=state)
                    control = local.regional_screen(7/8, float(alpha), 2, q=256, width=128, state=state, max_seconds=10)
                    for key, values in region.items():
                        np.testing.assert_allclose([log_float(x) for x in values], np.array(control[key]).ravel(),
                                                   rtol=0, atol=1e-9)
                    if alpha == 1:
                        unpaired = certificate.regional_directed(F(7, 8), alpha, 2, q=256, width=128,
                                                                 state=state, paired=False)
                        self.assertEqual(str(region), str(unpaired))

    def test_tiny_receipt_profiles_metadata_and_scalar_ceiling(self):
        # q256/r2 is a tiny harness check, not a GF128 selected-weight gate.
        row = certificate.weight_receipt(2, F(7, 8), F(1, 4), q=256, width=128, state=48,
                                         distance=F(12, 25), max_seconds=20)
        control = local.weight_screen(2, 7/8, .25, q=256, width=128, state=48,
                                      distance=F(12, 25), max_seconds=10)
        self.assertLess(abs(float(F(row['log2_upper_rational'])) - control['log2_upper_screen']), 1e-7)
        self.assertEqual(row['profiles'], sum(bool(base.profile_count_upper(4096, 2, a, b))
                                               for a in range(2) for b in range(2)))
        for key, value in dict(inner_width=128, state_dimension=48, emission_minimum=81,
                               pair_slots_per_packet=64, pair_slots_per_region=4096,
                               packets_per_region=64, composed_pairs_per_region=32, regions=256,
                               fixed_map_field_order=256, projective_divisor='255').items():
            self.assertEqual(row[key], value)
        self.assertEqual(row['construction'], local.construction(128, 48))
        self.assertNotEqual(row['construction'], local.construction(128, 64))
        self.assertTrue(row['covers_all_supports_at_weight'])
        self.assertTrue(row['covers_all_even_odd_splits'])
        self.assertTrue(row['directed'])
        self.assertFalse(row['whole_code_certificate'])
        with ctx.workprec(256):
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())

    def test_carryfree_helpers_guards_and_complete_pins(self):
        with ctx.workprec(256):
            encoded = base.encode({(2, 0): tuple(map(arb, (1, 0, 0, 1)))}, 2, F(1))
            squared = base.polynomial_product(encoded, encoded, 2, base.Budget(10))
            self.assertTrue(all(poly[index].is_zero() for poly in squared for index in range(6)))
        for kwargs in (dict(width=32), dict(width=64, state=33), dict(state=65), dict(state=0),
                       dict(q=16), dict(bits=8), dict(paired=1), dict(joint=1)):
            with self.assertRaises(ValueError):
                certificate.regional_directed(F(7, 8), F(1, 4), 2, **kwargs)
        for maximum in (0, 65):
            with self.assertRaises(ValueError):
                certificate.regional_directed(F(7, 8), F(1, 4), maximum)
        with self.assertRaises(ValueError):
            certificate.regional_directed(.875, F(1, 4), 2)
        with self.assertRaises(TimeoutError):
            certificate.packet_table(F(7, 8), 1, budget=base.Budget(1e-12))
        with self.assertRaises(RuntimeError):
            certificate.weight_receipt(1, F(7, 8), F(1, 4), q=256, _pins={})
        pins = certificate.source_pins()
        for name in ('hd_field128_wide_rs_state_certificate.py', 'hd_field128_wide_rs_state_screen.py',
                     'hd_field128_full_rs_state_screen.py', 'hd_field128_lift_certificate.py',
                     'hd_shifted_pairs_count.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
