"""Independent small exhaustive checks, plus the candidate's finite algebra."""
from itertools import product
from math import comb
import unittest

import numpy as np

import maps
import screen


def tiny():
    return maps.prepare_maps((15, 10, 12), (1, 2, 4, 3), packet_bits=2)[0]


class FiniteTests(unittest.TestCase):
    def test_weighted_syndromes_and_emissions(self):
        data = tiny()
        for z in (.39, .83, 1.):
            weighted, emission, _ = screen.moments(data, z)
            direct_w, direct_m = np.zeros((8, 3)), np.zeros((8, 3))
            for x in range(16):
                j = int(bool(x & 3))+int(bool(x >> 2))
                denominator = comb(2, j)*3**j
                direct_w[maps.apply(data['columns'], x), j] += z**x.bit_count()/denominator
                for a in range(8):
                    direct_m[a, j] += z**(x ^ int(data['images'][a])).bit_count()/denominator
            np.testing.assert_allclose(weighted, direct_w, atol=2e-15, rtol=2e-14)
            np.testing.assert_allclose(emission, direct_m, atol=2e-15, rtol=2e-14)

    def test_pointwise_transition_domination(self):
        data, z = tiny(), .71
        weighted, emission, _ = screen.moments(data, z)
        local = screen.operators(weighted, emission)
        born = weighted.copy()
        born[0, :] = 0
        sources = [np.eye(8)[0], np.r_[0., np.ones(7)/7],
                   born[:, 1]/born[:, 1].sum(), born[:, 2]/born[:, 2].sum()]
        for j in range(3):
            for index, source in enumerate(sources):
                direct = np.zeros(8)
                for x in range(16):
                    if int(bool(x & 3))+int(bool(x >> 2)) != j:
                        continue
                    syndrome = maps.apply(data['columns'], x)
                    for a, mass in enumerate(source):
                        value = mass*z**(x ^ int(data['images'][a])).bit_count()/(comb(2, j)*3**j)
                        if a == 0:
                            direct[syndrome] += value
                        else:
                            for ma in range(1, 8):
                                direct[ma ^ syndrome] += value/7
                represented = sum(local[j, index, k]*mu for k, mu in enumerate(sources))
                self.assertTrue(np.all(direct <= represented+2e-14))
        self.assertTrue(np.all(local[0, 1:, 0] == 0))

    def test_ordered_placement(self):
        local = np.array([[[1., 2.], [0., 1.]], [[1., 0.], [3., 1.]], [[2., 1.], [1., 3.]]])
        direct = np.zeros((7, 2, 2))
        for occupancies in product(range(3), repeat=3):
            value = np.eye(2)
            weight = 1
            for j in occupancies:
                value = value @ local[j]
                weight *= comb(2, j)
            total = sum(occupancies)
            direct[total] += weight/comb(6, total)*value
        fast, _ = screen.placement(local, 6, epochs=3, windows=2)
        slow, _ = screen.placement(local, 6, epochs=3, windows=2, force_log=True)
        np.testing.assert_allclose(np.exp(fast), direct, atol=2e-13, rtol=2e-13)
        np.testing.assert_allclose(np.exp(slow), direct, atol=2e-13, rtol=2e-13)

    def test_uniform_packet_mixture(self):
        regional = np.arange(1., 13.).reshape(3, 2, 2)
        for bits in (2, 8):
            p = 1-2.**(-bits)
            direct = sum(comb(2, j)*p**j*(1-p)**(2-j)*regional[j] for j in range(3))
            actual = np.exp(screen.uniform_mixture(np.log(regional), 2, packet_bits=bits))
            np.testing.assert_allclose(actual, direct, atol=1e-13, rtol=1e-13)

    def test_q1_ordered_region_supports(self):
        regional = np.array([[[.8, .2], [0., .4]], [[.1, .6], [.3, .2]]])
        expected = np.zeros(4)
        for hits in product(range(2), repeat=3):
            row = np.array([1., 0.])
            for hit in hits:
                row = row @ regional[hit]
            expected[sum(hits)] += row.sum()/comb(3, sum(hits))
        actual = np.exp(screen.q1_support_logs(np.log(regional, where=regional > 0,
            out=np.full_like(regional, -np.inf)), regions=3, tilt=.1, cutoff=0))
        np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=1e-13)


class CandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.record = maps.prepare()

    def test_literal_basis_and_ranks(self):
        data, record = self.data, self.record
        self.assertTrue(record['CA_zero'])
        self.assertEqual(record['packet_restriction_rank_counts'][1], {8: 8})
        self.assertEqual(record['packet_restriction_rank_counts'][2], {16: 28})
        for bit in range(16):
            state = 1 << bit
            for h in range(8):
                value = (state & 255) ^ maps.multiply(h, state >> 8)
                self.assertEqual((int(data['images'][state]) >> (8*h)) & 255, value)
        self.assertEqual(record['expansion_spectrum'][0], 1)
        self.assertEqual(record['expansion_spectrum'][8], 8)
        self.assertEqual(sum(record['expansion_spectrum'].values()), 65536)

    def test_candidate_mass_and_structural_zeros(self):
        weighted, emission, _ = screen.moments(self.data, .89)
        packet = ((1+.89)**8-1)/255
        np.testing.assert_allclose(weighted.sum(axis=0), packet**np.arange(9), atol=1e-14, rtol=1e-11)
        self.assertEqual(weighted[0, 1], 0)
        self.assertLess(abs(weighted[0, 2]), 1e-15)
        local = screen.operators(weighted, emission)
        self.assertEqual(local.shape, (9, 10, 10))
        self.assertTrue(np.isfinite(local).all())

    def test_outer_geometry(self):
        counts = screen.rs_outer.expected_group_support_counts(16, 8, 8, 2)
        self.assertEqual(len(counts), 33)
        self.assertEqual(sum(counts), (1 << 128)-1)
        self.assertEqual(counts[:9], (0,)*9)
        envelope = screen.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
        self.assertEqual(envelope.regions, 32)
        self.assertEqual(envelope.beta, screen.Fraction(2**256, 65535**8))
        self.assertTrue(all(x <= y for x, y in zip(counts, envelope.shell_caps())))

    def test_A_only_scale_is_a_distinct_map(self):
        scaled, record = maps.prepare('byte_native_A_scaled')
        self.assertEqual(scaled['columns'], self.data['columns'])
        self.assertNotEqual(scaled['rows'], self.data['rows'])
        self.assertFalse(record['CA_zero'])
        self.assertEqual(min(w for w in record['expansion_spectrum'] if w), 17)
        self.assertEqual(record['expansion_spectrum'][17], 3)
        for i, scale in enumerate(maps.SCALES):
            self.assertEqual(scale, maps.multiply(maps.SCALES[i-1], 3) if i else 1)
        for bit in range(16):
            state = 1 << bit
            for h, scale in enumerate(maps.SCALES):
                expected = maps.multiply(scale, (state & 255) ^ maps.multiply(h, state >> 8))
                self.assertEqual((int(scaled['images'][state]) >> (8*h)) & 255, expected)

    def test_t32_changes_the_actual_maps_and_geometry(self):
        data, record = maps.prepare('byte_native_t32')
        self.assertEqual(len(data['columns']), 32)
        self.assertEqual(data['windows'], 4)
        self.assertEqual(record['packet_restriction_rank_counts'][1], {8: 4})
        self.assertEqual(record['packet_restriction_rank_counts'][2], {16: 6})
        self.assertTrue(record['CA_zero'])
        self.assertTrue(all(int(x) < 1 << 32 for x in data['images']))
        shape = screen.geometry(data)
        self.assertEqual(shape['physical_steps_per_region'], 128)
        self.assertEqual(shape['physical_t'], 32)
        self.assertEqual(32*128*32, shape['N'])
        w, m, _ = screen.moments(data, .85)
        self.assertEqual(screen.operators(w, m).shape, (5, 6, 6))
        self.assertEqual(w[0, 2], 0)

    def test_zero_path_is_a_subcontribution(self):
        w, m, _ = screen.moments(tiny(), .73)
        local = screen.operators(w, m)
        full, _ = screen.placement(local, 6, epochs=3, windows=2)
        zero, _ = screen.placement(local[:, :1, :1], 6, epochs=3, windows=2)
        for q in range(1, 7):
            full_mix = screen.uniform_mixture(full, q, packet_bits=2)
            zero_mix = float(screen.uniform_mixture(zero, q, packet_bits=2)[0, 0])
            actual = screen.logarithmic.log_power_matrix(full_mix, 4)
            self.assertGreaterEqual(actual+1e-12, 4*zero_mix)


if __name__ == '__main__':
    unittest.main()
