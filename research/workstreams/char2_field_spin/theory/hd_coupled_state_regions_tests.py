"""Exact tiny coefficient and placement tests; no production coupled bank."""
from itertools import product
from math import comb
from unittest.mock import patch
import unittest

import numpy as np

import hd_coupled_state_regions as coupled


def packet(terms, maximum, full_max=None, marker=1.):
    full_max = coupled.checked_scope(maximum, full_max)
    result = np.zeros((5, 5, maximum+1, full_max+1))
    for single, full, matrix in terms:
        if full <= full_max:
            result[:, :, single, full] += marker**(single+full)*matrix
    return result


def integer_placements(terms, count, maximum):
    """Independent exact object-integer enumeration of chronological words."""
    result = {}
    for word in product(terms, repeat=count):
        single, full = sum(x[0] for x in word), sum(x[1] for x in word)
        if single+2*full > maximum:
            continue
        value = np.eye(5, dtype=object)
        for _, _, matrix in word:
            value = value@np.asarray(matrix, dtype=object)
        key = single, full
        result[key] = result.get(key, np.zeros((5, 5), dtype=object))+value
    return result


class CoupledRegionsTests(unittest.TestCase):
    def setUp(self):
        self.empty = np.eye(5, dtype=int)
        self.empty[0, 1] = 1
        self.single = np.eye(5, dtype=int)
        self.single[1, 0], self.single[2, 4] = 2, 1
        self.full = np.eye(5, dtype=int)
        self.full[4, 2], self.full[3, 1] = 3, 1
        self.terms = [(0, 0, self.empty), (1, 0, self.single), (0, 1, self.full)]

    def test_noncommuting_exact_placements_and_global_normalization(self):
        self.assertFalse(np.array_equal(self.empty@self.single, self.single@self.empty))
        maximum, packets, marker = 4, 4, .25
        actual = coupled.regions_from_packet(packet(self.terms, maximum, marker=marker),
            maximum, packets=packets, pairs_per_packet=1, marker=marker)['regions']
        expected = integer_placements(self.terms, packets, maximum)
        for (single, full), coefficient in expected.items():
            denominator = comb(packets, single)*comb(packets-single, full)
            np.testing.assert_allclose(np.exp2(actual[single+2*full, full]),
                np.asarray(coefficient, dtype=float)/denominator, rtol=2e-14, atol=1e-14)

    def test_generic_product_preserves_order_and_zero_blocks(self):
        left = packet([(0, 0, self.empty), (1, 0, self.single)], 2)
        right = packet([(0, 0, self.full), (1, 0, self.empty)], 2)
        actual = coupled.product_coefficients(left, right, 2)
        np.testing.assert_array_equal(actual[:, :, 0, 0], self.empty@self.full)
        np.testing.assert_array_equal(actual[:, :, 1, 0], self.empty@self.empty+self.single@self.full)
        np.testing.assert_array_equal(actual[:, :, 2, 0], self.single@self.empty)
        reverse = coupled.product_coefficients(right, left, 2)
        self.assertFalse(np.array_equal(actual, reverse))
        zero = np.zeros_like(left)
        np.testing.assert_array_equal(coupled.product_coefficients(left, zero, 2), zero)

    def test_full_pair_axis_truncation_is_exact(self):
        maximum = 4
        wide = coupled.power_coefficients(packet(self.terms, maximum), 4, maximum)
        narrow = coupled.power_coefficients(packet(self.terms, maximum, 1), 4, maximum, full_max=1)
        np.testing.assert_array_equal(narrow, wide[:, :, :, :2])
        singles = coupled.power_coefficients(packet(self.terms, maximum, 0), 4, maximum, full_max=0)
        np.testing.assert_array_equal(singles, wide[:, :, :, :1])

    def test_orientation_and_marker_once(self):
        table = {'e': self.empty, 'left': self.single, 'right': self.full}
        geometry = [(0, 0, 1, 'e'), (1, 0, 1, 'left'), (1, 0, 1, 'right'), (0, 1, 1, 'left')]
        result = coupled.packet_from_table(table, 2, marker=.5, geometry=geometry)
        np.testing.assert_array_equal(result[:, :, 0, 0], self.empty)
        np.testing.assert_array_equal(result[:, :, 1, 0], (self.single+self.full)/4)
        np.testing.assert_array_equal(result[:, :, 0, 1], self.single/2)

    def test_multinomial_identity_marker_cancellation(self):
        maximum, slots, packets = 4, 2, 4
        for marker in (1., .125):
            terms = [(s, f, comb(slots, s)*comb(slots-s, f)*np.eye(5))
                     for f in range(3) for s in range(maximum-2*f+1) if s+f <= slots]
            result = coupled.regions_from_packet(packet(terms, maximum, marker=marker), maximum,
                packets=packets, pairs_per_packet=slots, marker=marker)['regions']
            for matrix in result.values():
                np.testing.assert_allclose(np.exp2(matrix), np.eye(5), atol=2e-14)

    def test_identity_and_generic_dimension(self):
        base = np.zeros((3, 3, 3, 2))
        base[:, :, 0, 0] = [[1, 2, 0], [0, 1, 3], [4, 0, 1]]
        identity = coupled.power_coefficients(base, 0, 2)
        np.testing.assert_array_equal(identity[:, :, 0, 0], np.eye(3))
        np.testing.assert_array_equal(coupled.power_coefficients(base, 1, 2), base)
        squared = coupled.power_coefficients(base, 2, 2)
        np.testing.assert_array_equal(squared[:, :, 0, 0], base[:, :, 0, 0]@base[:, :, 0, 0])

    def test_public_provider_builds_one_table_and_retains_scope(self):
        geometry = [(0, 0, 1, 'e'), (1, 0, 2, 's'), (0, 1, 1, 'f')]
        table = {'e': np.eye(5), 's': np.eye(5), 'f': np.eye(5)}
        with patch.object(coupled, '_packet_table', return_value=table) as build, \
                patch.object(coupled, 'sparse_geometry', return_value=geometry):
            result = coupled.raw_region_bank(.9, 2, full_max=0)
        build.assert_called_once_with(.9)
        self.assertEqual(set(result['regions']), {(0, 0), (1, 0), (2, 0)})
        self.assertEqual(result['metadata']['state_order'], ['Z', 'L', 'R', 'B', 'P'])
        self.assertEqual(result['metadata']['global_normalization_count'], 1)
        self.assertEqual(result['metadata']['total_packets'], 128)
        self.assertEqual(result['metadata']['total_pairs'], 4096)
        self.assertFalse(result['metadata']['tensor_factorization'])

    def test_guards_fail_closed(self):
        for maximum in (-1, 65, True):
            with self.assertRaises(ValueError):
                coupled.checked_scope(maximum)
        for full_max in (-1, 3, True):
            with self.assertRaises(ValueError):
                coupled.checked_scope(4, full_max)
        base = packet(self.terms, 2)
        for exponent in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                coupled.power_coefficients(base, exponent, 2)
        for bad in (-1., np.nan, np.inf):
            malformed = base.copy()
            malformed[0, 0, 0, 0] = bad
            with self.assertRaises(ValueError):
                coupled.power_coefficients(malformed, 2, 2)
        outside = base.copy()
        outside[0, 0, 2, 1] = 1
        with self.assertRaises(ValueError):
            coupled.power_coefficients(outside, 2, 2)


if __name__ == '__main__':
    unittest.main()
