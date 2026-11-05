"""Small exact regressions; no bank, coefficient table, file write or replay.

Run remotely: OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -B -m unittest hd_field128_translated_packet_v50_tests -v
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import hd_field128_translated_packet_v50 as candidate


class Tests(unittest.TestCase):
    def test_original_denominator_by_independent_counting(self):
        # Syndrome = the first u free coordinates; its first b coordinates
        # vanish. Every active coordinate is a fixed function of the free
        # vector. Enumerate ORIGINAL inputs, including the forced positions.
        for q, j, c, u, b in ((2, 4, 1, 2, 1), (4, 4, 2, 2, 0),
                              (2, 5, 2, 2, 1), (4, 3, 1, 2, 2),
                              (2, 3, 0, 2, 1)):
            accepted = 0
            for vector in product(range(q), repeat=j):
                active, free = vector[:c], vector[c:]
                syndrome = free[:u]
                forced = tuple(free[i % len(free)] for i in range(c))
                if (active == forced and any(syndrome)
                        and all(x == 0 for x in syndrome[:b])):
                    accepted += 1
            self.assertEqual(candidate.probability_from_ranks(j, c, u, b, q),
                             F(accepted, q**j))
        self.assertEqual(candidate.probability_from_ranks(4, 1, 2, 1, 2), F(1, 8))

    def test_rank_floor_and_reversal(self):
        # A small deterministic spread, not all field elements or occupancies.
        cases = ((0, 0, 0), (1, 0, 31), (1, 0, 32), (16, 8, 16),
                 (24, 16, 24), (32, 31, 31), (32, 16, 16),
                 (48, 24, 16), (63, 62, 1), (64, 63, 0), (64, 64, 0))
        for j, c, e in cases:
            t, v, u = j-c, 64-j-e, min(32, j-c)
            if e >= 32:
                expected = u
            else:
                h = max(0, e+u-32, min((e+1)//2, t), min(e, (t+1)//2))
                loss = min(c, v, 32-h)//16 if h else 0
                expected = max(0, e+u-32, h-loss)
            b = candidate.rank_floor(j, c, e)
            self.assertEqual(b, expected)
            self.assertEqual(b, candidate.rank_floor(t+v, v, e))
            self.assertTrue(0 <= b <= u)
            self.assertEqual(candidate.return_probability(j, c, e),
                             candidate.reversal.probability(j, c, e, 1 << 128))
        self.assertEqual(candidate.return_probability(32, 31, 31), 0)
        self.assertEqual(candidate.return_probability(32, 32, 0), 0)
        self.assertEqual(candidate.return_probability(1, 0, 32), 0)

    def test_positive_expansion_matches_independent_nested_sum(self):
        q, z = 1 << 128, F(1, 4)
        for j in (0, 1, 2):
            expected = F(0)
            for c in range(j):
                u = min(32, j-c)
                for e in range(min(32, 65-j)):
                    b = candidate.rank_floor(j, c, e)
                    # Equivalent q^(-c-b)*(1-q^(-(u-b))) form is independent
                    # of the numerator/preimage implementation under test.
                    probability = F(1, q**(c+b))*(1-F(1, q**(u-b)))
                    expected += (comb(j, c)*comb(64-j, e)*z**(64-c-e)
                                 *(1-z)**(c+e)*probability)
            self.assertEqual(candidate.return_cap(j, z), expected)

    def test_zero_occupancy_and_sampled_real_moments(self):
        z, q = F(1, 2), 1 << 128
        for j in (0, 1, 2):
            baseline = candidate.generic.exact_moments(j, z, q, 64, 32, False)
            actual = candidate.exact_moments(j, z)
            self.assertEqual((actual[0], actual[1], actual[3]),
                             (baseline[0], baseline[1], baseline[3]))
            self.assertEqual(actual[2], min(baseline[2], candidate.return_cap(j, z)))
            self.assertTrue(all(isinstance(x, F) and x >= 0 for x in actual))
        self.assertEqual(candidate.exact_moments(0, z)[:3], (F(1), F(0), F(0)))
        self.assertEqual(candidate.exact_packet(0, z)[:3], (F(1), F(0), F(0)))

    def test_packet_normalizations_independently(self):
        q, d, z = 1 << 128, (1 << 4096)-1, F(1, 2)
        density = F(d, (q-1)**32)
        local = (F(1, 3), F(1, 5), F(1, 7), F(1, 11))
        with patch.object(candidate, 'exact_moments', return_value=local):
            self.assertEqual(candidate.exact_packet(1, z),
                             (local[0], density*local[1], local[2]/d, density*local[3]))
            self.assertEqual(candidate.exact_packet(0, z)[3], z**33)
        # Also exercise the other branch of min at j=0.
        with patch.object(candidate, 'exact_moments', return_value=(F(1), F(0), F(0), F(0))):
            self.assertEqual(candidate.exact_packet(0, z)[3], 0)

    def test_setup_risk_and_scope(self):
        q = 1 << 128
        risk = candidate.principal_bad_gamma_upper()
        self.assertEqual(risk, F(8*comb(64, 32), q-(1 << 64)))
        self.assertLess(risk, F(1, 1 << 64))
        with patch.object(candidate, 'source_pins', return_value={'test': 'pin'}):
            meta = candidate.metadata()
        self.assertEqual(meta['emission_nodes'], [128, 159])
        self.assertEqual(meta['feedback_nodes'], [64, 95])
        self.assertEqual(meta['state_dimension'], 32)
        self.assertEqual(meta['principal_bad_gamma_upper'], str(risk))
        self.assertFalse(meta['setup_risk_included_in_packet'])
        self.assertFalse(meta['additional_gamma_exclusions'])
        self.assertFalse(meta['distance_threshold_included'])
        self.assertFalse(meta['directed_certificate'])
        self.assertFalse(meta['whole_code_certificate'])
        self.assertTrue(meta['conditional_pointwise_in_every_good_gamma'])

    def test_frozen_source_identity(self):
        pins = candidate.source_pins()
        self.assertIn(candidate.NAME, pins)
        for name, digest in candidate.FROZEN_SHA256.items():
            self.assertEqual(pins[name], digest)
        candidate.check_source_pins(pins)
        with patch.dict(candidate.FROZEN_SHA256,
                        {'hd_field128_reversed_return_v44.py': '0'*64}):
            with self.assertRaises(RuntimeError):
                candidate.verify_frozen()
        with patch.object(candidate.generic, 'source_pins',
                          return_value={'hd_field128_reversed_return_v44.py': '0'*64}):
            with self.assertRaises(RuntimeError):
                candidate.source_pins()

    def test_reject_inexact_inputs_and_invalid_geometry(self):
        for j, z in ((True, F(1, 2)), (-1, F(1, 2)), (65, F(1, 2)),
                     (1, .5), (1, True), (1, 0), (1, 1)):
            with self.assertRaises(ValueError):
                candidate.exact_moments(j, z)
        for args in ((2, -1, 0), (2, 3, 0), (2, 0, 63), (2, 0, True)):
            with self.assertRaises(ValueError):
                candidate.rank_floor(*args)
        for args in ((3, 1, 3, 0, 2), (3, 1, 2, 3, 2), (3, 1, 2, 0, 3)):
            with self.assertRaises(ValueError):
                candidate.probability_from_ranks(*args)


if __name__ == '__main__':
    unittest.main()
