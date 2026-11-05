import unittest
from math import log2

import numpy as np

import hd_field128_gamma_near_face_counts as near
from hd_shifted_pairs_count import exact_profiles, profile_count_upper


class Expired:
    def check(self):
        raise TimeoutError('test budget')


class NearFaceTests(unittest.TestCase):
    def test_exhaustive_small_cycles_and_scalar_formula(self):
        for pairs in range(1, 8):
            actual = exact_profiles(pairs)
            for weight in range(2*pairs+1):
                counts = near.NearFaceCounts(pairs, weight, max_singletons=2*pairs)
                for first in counts.full:
                    for second in counts.full:
                        value = counts.count(first, second)
                        old = profile_count_upper(pairs, weight, first, second)
                        runs = weight-first-second
                        both = [near.near_face_count_upper(pairs, weight, f, runs) for f in (first, second)]
                        self.assertEqual(value, min(old, *both))
                        self.assertGreaterEqual(value, actual.get((weight, first, second), 0))
                        self.assertLessEqual(value, old)
                        self.assertEqual(value, counts.count(second, first))
                        if 2*first == weight or 2*second == weight:
                            self.assertEqual(value, actual.get((weight, first, second), 0))

    def test_constants_empty_remnant_invalid_and_cutoff(self):
        self.assertEqual(near.near_face_count_upper(4, 0, 0, 0), 1)
        self.assertEqual(near.near_face_count_upper(4, 8, 4, 0), 1)
        self.assertEqual(near.near_face_count_upper(4, 4, 0, 4), 16)
        self.assertEqual(near.near_face_count_upper(4, 4, 0, 5), 0)
        self.assertEqual(near.near_face_count_upper(4, 2, 3, 0), 0)
        for pairs in range(1, 6):
            for weight in range(2*pairs+1):
                counts = near.NearFaceCounts(pairs, weight, max_singletons=0)
                for first in counts.full:
                    for second in counts.full:
                        self.assertEqual(counts.count(first, second), profile_count_upper(pairs, weight, first, second))
                self.assertEqual(counts.count(-1, 0), 0)
        for args in ((0, 0), (4, -1), (4, 9)):
            with self.assertRaises(ValueError):
                near.NearFaceCounts(*args)

    def test_vector_logs_are_logs_of_exact_counts(self):
        for weight in range(13):
            counts = near.NearFaceCounts(6, weight, max_singletons=6)
            full = np.array(list(counts.full))
            logs = counts.log_counts(full[:, None], full[None])
            for i, first in enumerate(full):
                for j, second in enumerate(full):
                    value = counts.count(int(first), int(second))
                    self.assertAlmostEqual(logs[i, j], log2(value) if value else -np.inf)
                    scalar = counts.log_counts(int(first), int(second))
                    self.assertAlmostEqual(float(scalar), log2(value) if value else -np.inf)

    def test_transfer_matches_independent_noncommuting_matrix_oracle(self):
        rng = np.random.default_rng(73129)
        for weight in (2, 3, 4, 6):
            counts = near.NearFaceCounts(4, weight, max_singletons=4)
            matrices = rng.uniform(.001, .2, (len(counts.full), 2, 2))
            total = 0.
            for i, first in enumerate(counts.full):
                for j, second in enumerate(counts.full):
                    power = np.linalg.matrix_power(matrices[i] @ matrices[j], 4)
                    total += counts.count(first, second)*sum(power[0])
            result = near.support_sum(4, weight, np.log2(matrices), repetitions=8, max_singletons=4, chunk=2)
            old = near.provider.profiles.support_sum(4, weight, np.log2(matrices), repetitions=8)
            self.assertAlmostEqual(result['log2_support_sum_screen'], log2(total), places=11)
            self.assertLessEqual(result['log2_support_sum_screen'], old['log2_support_sum_screen']+1e-10)
            self.assertEqual(result['profile_pairs_evaluated'], len(counts.full)**2)
            self.assertFalse(result['directed'])

    def test_timeout_and_source_closure(self):
        with self.assertRaises(TimeoutError):
            near.NearFaceCounts(4, 4, budget=Expired())
        with self.assertRaises(TimeoutError):
            near.support_sum(4, 4, np.zeros((3, 2, 2)), budget=Expired())
        with self.assertRaises(ValueError):
            near.support_sum(4, 4, np.full((3, 2, 2), np.nan))
        pins = near.source_pins()
        for name in ('hd_field128_gamma_near_face_counts.py', 'hd_field128_gamma_profile_certificate.py',
                     'hd_shifted_pairs_count.py', 'hd_field128_gamma_all_pair_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
