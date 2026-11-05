"""Small independent exact tests. No artifacts, coefficient bank or replay.

Root may run remotely with OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -B -m unittest hd_field128_return_tail_moments_v57_tests -v
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest
from unittest.mock import patch

import hd_field128_return_tail_moments_v57 as candidate


def moments(distribution):
    n = len(distribution)-1
    return tuple(sum((mass*comb(k, l) for k, mass in enumerate(distribution)
                      if k >= l), F(0)) for l in range(n+1))


class Tests(unittest.TestCase):
    def test_independent_finite_distributions(self):
        # Direct probability enumeration, not the implementation's expansion.
        cases = ((F(1),), (F(0), F(0), F(0)),
                 (F(1, 4), F(1, 8), F(1, 16), F(0)),
                 (F(1, 7), F(2, 7), F(0), F(3, 7), F(1, 7)),
                 (F(0), F(0), F(1), F(0), F(0), F(0)))
        for distribution in cases:
            bound = moments(distribution)
            maximum, clipped, tails = candidate.factorial_envelope(bound)
            self.assertEqual(clipped, bound)
            for k in range(maximum+1):
                self.assertGreaterEqual(tails[k], sum(distribution[k:], F(0)))
            for z in (F(1, 4), F(1, 2), F(3, 4)):
                actual = sum((p*z**(len(distribution)-1-k)
                              for k, p in enumerate(distribution)), F(0))
                caps = candidate.evaluate_envelope(bound, z)
                self.assertEqual(caps['original_bernstein'], actual)
                self.assertEqual(caps['refined'], actual)
                self.assertGreaterEqual(caps['tail'], actual)

    def test_strict_clipping_and_exact_zero_cutoff(self):
        # Valid upper moments for the deterministic law K=2, with loose H1.
        bound = (F(1), F(10), F(1), F(0), F(0))
        maximum, clipped, tails = candidate.factorial_envelope(bound)
        self.assertEqual(maximum, 2)
        self.assertEqual(clipped, (F(1), F(2), F(1), F(0), F(0)))
        self.assertEqual(tails, (F(1), F(1), F(1)))
        caps = candidate.evaluate_envelope(bound, F(1, 2))
        self.assertEqual(caps['original_bernstein'], F(3, 4))
        self.assertEqual(caps['clipped_bernstein'], F(1, 4))
        self.assertEqual(caps['tail'], F(1, 4))
        # No epsilon-based truncation of a positive exact moment is allowed.
        tiny = F(1, 1 << 256)
        self.assertEqual(candidate.factorial_envelope((F(1), tiny))[0], 1)

    def test_rank_plateau_clips_later_moment(self):
        bound = (F(1), F(1, 8), F(1, 8), F(0))
        maximum, clipped, unused = candidate.factorial_envelope(bound)
        self.assertEqual(maximum, 2)
        self.assertEqual(clipped[2], F(1, 16))
        distribution = (F(15, 16), F(0), F(1, 16), F(0))
        for z in (F(1, 3), F(2, 3)):
            actual = sum((p*z**(3-k) for k, p in enumerate(distribution)), F(0))
            caps = candidate.evaluate_envelope(bound, z)
            self.assertEqual(caps['refined'], actual)
            self.assertLess(caps['refined'], caps['original_bernstein'])

    def test_tiny_binary_linear_oracle_original_denominator(self):
        # Uniform X on exactly three potential coordinates, denominator 2^3.
        # Syndrome rank2; output is a fixed injective linear map. No MDS or
        # generic-construction assumption is needed for the moment identity.
        rows = []
        for x in product((0, 1), repeat=3):
            syndrome = (x[0]^x[1], x[1]^x[2])
            y = (x[0]^x[2], x[1], x[2], x[0]^x[1]^x[2])
            rows.append((bool(any(syndrome)), y))
        distribution = [F(0)]*5
        for nonzero, y in rows:
            if nonzero:
                distribution[y.count(0)] += F(1, 8)
        grouped = [F(0)]*5
        for c in range(4):
            for e in range(2):
                for active in combinations(range(3), c):
                    for inactive in combinations((3,), e):
                        fixed = active+inactive
                        grouped[c+e] += F(sum(nonzero and all(y[k] == 0 for k in fixed)
                                                   for nonzero, y in rows), 8)
        self.assertEqual(tuple(grouped), moments(distribution))
        self.assertEqual(grouped[0], F(3, 4))
        for z in (F(1, 2), F(1, 8)):
            direct = sum((F(1, 8)*z**sum(bool(v) for v in y)
                          for nonzero, y in rows if nonzero), F(0))
            self.assertEqual(candidate.evaluate_envelope(grouped, z)['refined'], direct)

    def test_small_production_projection_and_zero_mass(self):
        # Only j0/j1; this does not build 65 occupancies or any coefficients.
        z, q = F(1, 2), 1 << 128
        for j in (0, 1):
            old = candidate.provider.exact_moments(j, z)
            actual = candidate.exact_moments(j, z)
            h = candidate.factorial_bounds(j)
            self.assertEqual(h[0], F(0) if j == 0 else 1-F(1, q))
            self.assertEqual((actual[0], actual[1], actual[3]), (old[0], old[1], old[3]))
            self.assertLessEqual(actual[2], old[2])
            self.assertGreaterEqual(actual[2], z**64*h[0])
            caps = candidate.cap_components(j, z)
            self.assertEqual(caps['selected'], actual[2])
        self.assertEqual(candidate.return_cap(0, z), 0)
        self.assertEqual(candidate.exact_packet(0, z)[:3], (F(1), F(0), F(0)))

    def test_packet_normalization_and_setup_scope(self):
        q, d, z = 1 << 128, (1 << 4096)-1, F(1, 2)
        density = F(d, (q-1)**32)
        local = (F(1, 3), F(1, 5), F(1, 7), F(1, 11))
        with patch.object(candidate, 'exact_moments', return_value=local):
            self.assertEqual(candidate.exact_packet(1, z),
                             (local[0], density*local[1], local[2]/d, density*local[3]))
            self.assertEqual(candidate.exact_packet(0, z)[3], z**33)
        with patch.object(candidate, 'source_pins', return_value={'mock': 'pin'}):
            meta = candidate.metadata()
        self.assertEqual(meta['construction'], candidate.provider.CONSTRUCTION)
        self.assertFalse(meta['construction_changed'])
        self.assertFalse(meta['additional_gamma_exclusions'])
        self.assertFalse(meta['setup_risk_included_in_packet'])
        self.assertFalse(meta['numerical_improvement_measured'])
        self.assertFalse(meta['whole_code_certificate'])
        self.assertEqual(candidate.principal_bad_gamma_upper(),
                         candidate.provider.principal_bad_gamma_upper())

    def test_source_closure_and_guards(self):
        pins = candidate.source_pins()
        self.assertEqual(pins[candidate.provider.NAME], candidate.PROVIDER_SHA256)
        self.assertIn(candidate.NAME, pins)
        candidate.check_source_pins(pins)
        with patch.object(candidate, 'PROVIDER_SHA256', '0'*64):
            with self.assertRaises(RuntimeError):
                candidate.verify_frozen()
        for bounds in ((), (F(1),)*66, (F(-1),), (True,), (.5,)):
            with self.assertRaises(ValueError):
                candidate.factorial_envelope(bounds)
        for j, z in ((True, F(1, 2)), (-1, F(1, 2)), (65, F(1, 2)),
                     (1, .5), (1, True), (1, 0), (1, 1)):
            with self.assertRaises(ValueError):
                candidate.exact_moments(j, z)


if __name__ == '__main__':
    unittest.main()
