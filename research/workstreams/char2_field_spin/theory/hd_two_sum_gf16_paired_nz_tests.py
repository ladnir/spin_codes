"""Small explicit-table tests; no production NZ kernel or regional bank."""
from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_paired_nz as target
from hd_shifted_swap_sparse import power_screen, product_screen, sparse_geometry
from hd_two_sum_gf16_fractional_state import normalize_regional_coefficients


def toy_table(maximum, constant=False):
    result = {}
    for _, _, _, occupancy in sparse_geometry(maximum):
        if constant:
            result[occupancy] = np.array([1., 0., 0., 1.])
        else:
            degree = sum(occupancy)
            result[occupancy] = (np.array([1., 0., 0., .25]) if not degree else
                np.array([1/(degree+2), (1+occupancy[0])/16, (1+occupancy[-1])/32, 1/(degree+4)]))
    return result


class PairedNZTests(unittest.TestCase):
    def test_good_geometry_matches_swap_averaged_single_counts_without_extra_divisor(self):
        maximum = 4
        old = {}
        for s, f, count, occupancy in sparse_geometry(maximum):
            if s <= 1:
                key = s, f, occupancy
                old[key] = old.get(key, F(0))+F(count, 2**s)
        new = {}
        for s, f, count, occupancy in target.relevant_good_geometry(maximum):
            key = s, f, occupancy
            new[key] = new.get(key, 0)+count
        self.assertEqual(old, new)
        self.assertEqual(sum(n for s, f, n, _ in target.relevant_good_geometry(1) if s == 1), 32)

    def test_identity_table_counts_every_two_packet_geometry(self):
        maximum = 4
        pieces = target.paired_components_from_table(toy_table(maximum, True), F(1, 3), maximum, marker=1.)
        for s in range(maximum+1):
            for f in range(maximum//2+1):
                expected = comb(64, s)*comb(64-s, f) if s+2*f <= maximum else 0
                self.assertEqual(pieces['paired'][0, s, f], expected)
                self.assertEqual(pieces['paired'][3, s, f], expected)
                self.assertEqual(pieces['paired'][1, s, f], 0.)
                self.assertEqual(pieces['paired'][2, s, f], 0.)
        # Each good packet allows only0/1 singleton: 32*32, no extra1/4.
        self.assertEqual(pieces['composed'][0, 2, 0], 1024.)

    def test_alpha1_composition_equals_packetwise_coefficients_and_region(self):
        maximum, marker = 3, 1/512
        table = toy_table(maximum)
        pieces = target.paired_components_from_table(table, 1, maximum, marker=marker)
        np.testing.assert_allclose(pieces['composed'], product_screen(pieces['good'], pieces['good'], maximum),
                                   rtol=1e-14, atol=1e-20)
        np.testing.assert_allclose(pieces['paired'], product_screen(pieces['single'], pieces['single'], maximum),
                                   rtol=1e-14, atol=1e-20)
        paired = target.regions_from_table(table, 1, maximum, marker=marker)
        reference = normalize_regional_coefficients(power_screen(pieces['single'], 128, maximum), maximum, marker)
        for key in paired:
            np.testing.assert_allclose(paired[key], reference[key], atol=1e-11)

    def test_fractional_good_pair_and_region_never_exceed_packetwise_bound(self):
        maximum, marker = 3, 1/512
        for alpha in (F(1, 2), F(1, 5)):
            table = toy_table(maximum)
            pieces = target.paired_components_from_table(table, alpha, maximum, marker=marker)
            before = product_screen(pieces['single'], pieces['single'], maximum)
            self.assertTrue(np.all(pieces['paired'] <= before*(1+2e-14)+1e-30))
            self.assertLess(pieces['paired'][0, 2, 0], before[0, 2, 0])
            paired = target.regions_from_table(table, alpha, maximum, marker=marker)
            reference = normalize_regional_coefficients(power_screen(pieces['single'], 128, maximum), maximum, marker)
            for key in paired:
                self.assertTrue(np.all(paired[key] <= reference[key]+1e-10))

    def test_exact_two_state_concavity_and_noncommuting_composition(self):
        # Integer square roots make the alpha=1/2 right-hand side exact.
        left_root = np.array([[1, 2], [3, 4]], dtype=object)
        right_root = np.array([[5, 6], [7, 8]], dtype=object)
        left, right = left_root**2, right_root**2
        composed = left@right
        packetwise_squared = (left_root@right_root)**2
        self.assertTrue(np.all(composed <= packetwise_squared))
        self.assertFalse(np.array_equal(left@right, right@left))

    def test_explicit_ordered_good_products_against_reference_sum(self):
        maximum, alpha, marker = 3, F(1, 2), .25
        table = toy_table(maximum)
        expected = np.zeros((4, maximum+1, maximum//2+1))
        rows = target.relevant_good_geometry(maximum)
        for s, f, n, occ in rows:
            for u, v, m, other in rows:
                if s+u+2*(f+v) <= maximum:
                    product = table[occ].reshape(2, 2)@table[other].reshape(2, 2)
                    expected[:, s+u, f+v] += n*m*marker**(s+u+f+v)*np.sqrt(product).ravel()
        got = target.good_composition_from_table(table, alpha, maximum, marker)
        np.testing.assert_allclose(got, expected, rtol=1e-14)

    def test_guards_and_metadata_without_hidden_kernel_call(self):
        for maximum in (0, 65, True):
            with patch.object(target, 'nonzero_matrices') as kernel, self.assertRaises(ValueError):
                target.paired_nz_regions(F(99, 100), F(1, 5), maximum)
            kernel.assert_not_called()
        table = toy_table(2)
        table[(0, 0, 0, 0)] = np.array([1., -1., 0., 1.])
        with self.assertRaises(ValueError):
            target.paired_components_from_table(table, F(1, 5), 2)
        metadata = target.provider_metadata()
        self.assertEqual(metadata['paired_blocks']*metadata['composition_steps'], 128)
        self.assertFalse(metadata['extra_good_orientation_divisor'])
        self.assertFalse(metadata['route_denominator_powered'])
        self.assertFalse(metadata['directed'])


if __name__ == '__main__':
    unittest.main()
