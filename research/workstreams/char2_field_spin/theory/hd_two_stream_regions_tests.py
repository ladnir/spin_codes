"""Exact tiny noncommuting placement controls; no production NZ bank."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest.mock import patch
import unittest

import numpy as np

import hd_two_stream_regions as stream


def mm(a, b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) for i in range(2) for j in range(2))


def kron(a, b):
    return tuple(a[2*(i//2)+j//2]*b[2*(i%2)+j%2] for i in range(4) for j in range(4))


def exact_placements(terms, packets, maximum):
    one, two = {}, {}
    identity = (F(1), F(0), F(0), F(1))
    for choices in product(terms, repeat=packets):
        single = sum(term[0] for term in choices)
        full = sum(term[1] for term in choices)
        if single+2*full > maximum:
            continue
        sequential, halves = identity, [identity, identity]
        for index, (_, _, matrix) in enumerate(choices):
            sequential = mm(sequential, matrix)
            halves[index % 2] = mm(halves[index % 2], matrix)
        key = single, full
        one[key] = tuple(a+b for a, b in zip(one.get(key, (F(0),)*4), sequential))
        joint = kron(*halves)
        two[key] = tuple(a+b for a, b in zip(two.get(key, (F(0),)*16), joint))
    return one, two


def make_packet(terms, maximum, marker=F(1)):
    result = np.zeros((4, maximum+1, maximum//2+1))
    for single, full, matrix in terms:
        result[:, single, full] += np.asarray(matrix, dtype=float)*float(marker**(single+full))
    return result


class TwoStreamTests(unittest.TestCase):
    def test_noncommuting_chronological_packet_placement(self):
        empty, single, full = ((1, 1, 0, 1), (1, 0, 2, 1), (2, 1, 1, 1))
        self.assertNotEqual(mm(empty, single), mm(single, empty))
        terms = [(0, 0, empty), (1, 0, single), (0, 1, full)]
        maximum, packets = 4, 4
        expected_one, expected_two = exact_placements(terms, packets, maximum)
        packet = make_packet(terms, maximum, F(1, 4))
        actual = stream.controls_from_packet(packet, maximum, packets=packets,
                    pairs_per_packet=1, marker=.25)
        for (s, f), values in expected_one.items():
            denominator = comb(packets, s)*comb(packets-s, f)
            if not denominator:
                continue
            key = s+2*f, f
            expected = np.asarray(values, dtype=float).reshape(2, 2)/denominator
            np.testing.assert_allclose(np.exp2(actual['one_stream'][key]), expected, rtol=2e-14)
            expected = np.asarray(expected_two[s, f], dtype=float).reshape(4, 4)/denominator
            np.testing.assert_allclose(np.exp2(actual['two_stream'][key]), expected, rtol=2e-14)
        control = stream.power_screen(packet, packets, maximum)
        expected = stream.normalize_coefficients(control, maximum, pairs=packets, marker=.25)
        for key in expected:
            np.testing.assert_allclose(actual['one_stream'][key], expected[key], atol=1e-13)

    def test_singleton_orientation_average_exactly_once(self):
        table = {'e': (1, 0, 0, 1), 'l': (1, 2, 3, 4), 'r': (2, 1, 4, 3), 'f': (5, 6, 7, 8)}
        geometry = ((0, 0, 1, 'e'), (1, 0, 1, 'l'), (1, 0, 1, 'r'), (0, 1, 1, 'f'))
        packet = stream.packet_from_table(table, 2, marker=.5, geometry=geometry)
        np.testing.assert_array_equal(packet[:, 1, 0], (np.array(table['l'])+table['r'])/4)
        np.testing.assert_array_equal(packet[:, 0, 1], np.array(table['f'])/2)
        self.assertEqual(packet[:, 0, 0].tolist(), list(table['e']))

    def test_global_normalizer_identity_and_marker_cancellation(self):
        maximum, packets, slots = 4, 4, 2
        identity = np.array([1., 0., 0., 1.])
        baseline = None
        for marker in (1., .125):
            packet = np.zeros((4, maximum+1, maximum//2+1))
            for single, full in stream.valid_degrees(maximum):
                if single+full <= slots:
                    packet[:, single, full] = (comb(slots, single)*comb(slots-single, full)
                        *marker**(single+full)*identity)
            result = stream.controls_from_packet(packet, maximum, packets=packets,
                                                   pairs_per_packet=slots, marker=marker)
            for matrix in result['two_stream'].values():
                np.testing.assert_allclose(np.exp2(matrix), np.eye(4), atol=2e-15)
            if baseline is not None:
                for key in baseline:
                    np.testing.assert_allclose(result['two_stream'][key], baseline[key], atol=2e-14)
            baseline = result['two_stream']

    def test_joint_half_allocation_is_hypergeometric_not_resampled(self):
        #A nonempty packet changes0 to1; the state remembers whether this
        #stream received any occupied slot. Each packet has one pair slot.
        empty, active = (1, 0, 0, 1), (0, 1, 0, 1)
        terms = [(0, 0, empty), (1, 0, active), (0, 1, active)]
        maximum = 4
        result = stream.controls_from_packet(make_packet(terms, maximum), maximum,
                                               packets=4, pairs_per_packet=1)
        for (weight, full), matrix in result['two_stream'].items():
            occupied = weight-full
            if not occupied:
                continue
            one_half = F(comb(2, occupied), comb(4, occupied)) if occupied <= 2 else F(0)
            row = np.exp2(matrix[0])
            self.assertAlmostEqual(row[1]+row[2], float(2*one_half), places=14)
            self.assertAlmostEqual(sum(row), 1., places=14)
        #One full pair stays inside one stream; two singleton slots need not.
        self.assertEqual(np.exp2(result['two_stream'][2, 1][0, 3]), 0.)
        self.assertAlmostEqual(np.exp2(result['two_stream'][2, 0][0, 3]), 2/3)

    def test_exact_swap_symmetry_and_noncommuting_lumped_products(self):
        terms1 = [(0, 0, (1, 2, 0, 1)), (1, 0, (2, 0, 1, 1)), (0, 1, (1, 1, 3, 2))]
        terms2 = [(0, 0, (2, 1, 1, 1)), (1, 0, (1, 3, 2, 1)), (0, 1, (3, 1, 0, 2))]
        controls = [stream.controls_from_packet(make_packet(t, 2), 2, packets=2, pairs_per_packet=1)
                    for t in (terms1, terms2)]
        matrices = [np.exp2(c['two_stream'][2, 0]) for c in controls]
        lumped = [np.exp2(c['two_stream_lumped'][2, 0]) for c in controls]
        self.assertFalse(np.allclose(matrices[0]@matrices[1], matrices[1]@matrices[0]))
        for order in ((0, 1), (1, 0), (0, 1, 0)):
            full, reduced = np.eye(4), np.eye(3)
            for index in order:
                full = full@matrices[index]
                reduced = reduced@lumped[index]
            self.assertAlmostEqual(sum(full[0]), sum(reduced[0]), places=9)
        for control in controls:
            for key, matrix in control['two_stream'].items():
                swap = [0, 2, 1, 3]
                np.testing.assert_array_equal(matrix, matrix[np.ix_(swap, swap)])
                raw, quotient = np.exp2(matrix), np.exp2(control['two_stream_lumped'][key])
                average = (raw[1]+raw[2])/2
                np.testing.assert_allclose(quotient[1], [average[0], average[1]+average[2], average[3]])

    def test_provider_builds_one_packet_and_full_filter_only_changes_outputs(self):
        maximum = 2
        packet = make_packet([(0, 0, (1, 0, 0, 1)), (1, 0, (1, 1, 0, 1)),
                              (0, 1, (1, 0, 1, 1))], maximum, F(maximum, 4096))
        with patch.object(stream, 'raw_packet', return_value=packet) as builder:
            result = stream.raw_region_controls(.9, maximum, full_max=0)
        builder.assert_called_once_with(.9, maximum, marker=maximum/4096)
        self.assertEqual(set(result['two_stream']), {(0, 0), (1, 0), (2, 0)})
        self.assertEqual(result['metadata']['state_order'], ['00', '01', '10', '11'])
        self.assertEqual(result['metadata']['packets_per_stream'], 64)
        self.assertEqual(result['metadata']['global_normalization_count'], 1)
        self.assertFalse(result['metadata']['independent_half_normalization'])

    def test_guards(self):
        for maximum in (0, 65, True):
            with self.assertRaises(ValueError):
                stream.raw_region_controls(.9, maximum)
        for mark in (0., 1., np.nan, True):
            with self.assertRaises(ValueError):
                stream.raw_packet(mark, 2)
        with self.assertRaises(ValueError):
            stream.lump_stream_symmetry(np.log2(np.arange(1, 17).reshape(4, 4)))
        with self.assertRaises(ValueError):
            stream.controls_from_packet(np.zeros((4, 3, 2)), 2, packets=3)


if __name__ == '__main__':
    unittest.main()
