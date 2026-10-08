"""Independent exact tiny compositions and pooled-state interface checks."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

import mask_regional_screen as screen


class MaskRegionalTests(unittest.TestCase):
    def setUp(self):
        self.local = np.asarray([[[F(1), F(0)], [F(0), F(2, 3)]],
                                 [[F(1, 7), F(2, 5)], [F(1, 11), F(3, 4)]],
                                 [[F(0), F(4, 5)], [F(2, 7), F(1, 2)]]], dtype=object)

    def direct_block(self, packets):
        slots = 2*packets
        sums = np.zeros((slots+1, 2, 2), dtype=object)
        for bits in product((0, 1), repeat=slots):
            value = np.eye(2, dtype=object)
            for packet in range(packets):
                value = value @ self.local[sum(bits[2*packet:2*packet+2])]
            sums[sum(bits)] += value
        return np.asarray([value/F(comb(slots, h)) for h, value in enumerate(sums)], dtype=object)

    def test_exact_hypergeometric_distribution(self):
        for nleft, nright in ((0, 5), (2, 3), (17, 29), (100000, 100001)):
            for row in screen.hypergeometric_weights(nleft, nright, 6):
                self.assertEqual(sum(row), 1)
                self.assertTrue(all(0 <= value <= 1 for value in row))

    def test_power_against_independent_slot_subsets(self):
        expected = self.direct_block(3)
        actual = screen.normalized_power(self.local, 2, 3, 6, exact=True)
        np.testing.assert_array_equal(actual, expected)
        truncated = screen.normalized_power(self.local, 2, 3, 3, exact=True)
        np.testing.assert_array_equal(truncated, expected[:4])
        floating = screen.normalized_power(np.asarray(self.local, dtype=float), 2, 3, 6)
        np.testing.assert_allclose(floating, np.asarray(expected, dtype=float), rtol=3e-15)

    def test_noncommuting_block_order_and_raw_coefficients(self):
        right = self.local[:, ::-1, ::-1]
        actual = screen.normalized_product(self.local, 2, right, 2, 4, exact=True)
        raw = np.zeros((5, 2, 2), dtype=object)
        for i in range(3):
            for j in range(3):
                raw[i+j] += comb(2, i)*comb(2, j)*(self.local[i] @ right[j])
        expected = np.asarray([row/F(comb(4, h)) for h, row in enumerate(raw)], dtype=object)
        np.testing.assert_array_equal(actual, expected)
        reversed_order = screen.normalized_product(right, 2, self.local, 2, 4, exact=True)
        self.assertTrue((actual != reversed_order).any())

    def test_exact_selected_row_relaxation(self):
        # Select A regional positions, then independently sum every active
        # subset of those positions. This is NOT Bernoulli slot placement.
        slots, active, tau = 6, 3, F(7, 3)
        total = np.zeros((2, 2), dtype=object)
        for selected in combinations(range(slots), active):
            for bits in product((0, 1), repeat=active):
                support = {i for i, bit in zip(selected, bits) if bit}
                value = np.eye(2, dtype=object)
                for packet in range(3):
                    weight = len(support & {2*packet, 2*packet+1})
                    value = value @ self.local[weight]
                total += tau**len(support)*value/F(comb(slots, active))
        coefficients = screen.normalized_power(self.local, 2, 3, active, exact=True)
        expected = sum((comb(active, h)*tau**h*coefficients[h] for h in range(active+1)),
                       np.zeros((2, 2), dtype=object))
        np.testing.assert_array_equal(total, expected)

    @staticmethod
    def fake_mask(weight, y, tau, phase):
        value = np.eye(16)*(phase+1)
        if weight:
            value[0] = 0
            value[0, 15] = y**weight
        return value

    def test_pooled_fallback_and_supported_routing(self):
        with patch.object(screen, 'local_mask_transfer', self.fake_mask), \
             patch.object(screen, 'local_rank_transfer', return_value=[[2., 3.], [5., 7.]]), \
             patch.object(screen, 'rank_return_moment', return_value=4.):
            supported = screen.packet_matrix(1, .9, 35., cutoff=1)
            self.assertEqual(supported[0, 15], .9)
            self.assertTrue((supported[:16, 16] == 0).all())
            self.assertEqual(supported[16, 0], 4.)
            self.assertEqual(supported[16, 16], 7.)
            self.assertTrue((supported[16, 1:16] == 0).all())
            fallback = screen.packet_matrix(2, .9, 35., cutoff=1)
            self.assertEqual(fallback[0, 0], 2.)
            self.assertEqual(fallback[0, 16], 3.)
            self.assertTrue((fallback[:, 1:16] == 0).all())
            self.assertTrue((fallback[1:16, 0] == 4.).all())
            self.assertTrue((fallback[1:16, 16] == 7.).all())

    def test_empty_packet_and_phase_max(self):
        value = screen.packet_matrix(0, .99, 35.)
        self.assertTrue((value-np.diag(np.diag(value)) == 0).all())
        self.assertEqual(value[0, 0], 1.)
        with patch.object(screen, 'local_mask_transfer', self.fake_mask):
            maximum = screen.packet_coefficients(0, .99, 35.)[0]
            self.assertEqual(maximum[1, 1], 4.)
            for phase in range(4):
                self.assertTrue((maximum >= screen.packet_matrix(0, .99, 35., phase=phase)).all())

    def test_invalid_and_nonfinite_parameters(self):
        for y, tau in ((float('nan'), 35.), (.9, float('inf')), (0., 35.), (.9, 256.)):
            with self.assertRaises(ValueError):
                screen.packet_matrix(1, y, tau)
        with self.assertRaises(ValueError):
            screen.region_coefficients(3, .9, 35., 17)
        with self.assertRaises(ValueError):
            screen.normalized_product(self.local[:1], 2, self.local, 2, 4)


if __name__ == '__main__':
    unittest.main()
