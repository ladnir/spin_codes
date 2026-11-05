import itertools
import unittest
from collections import Counter, defaultdict
from math import comb
from hd_shifted_pairs_count import (exact_profiles, one_pairing_count,
                                    circular_run_count, profile_count_upper)


class ShiftedPairCountTests(unittest.TestCase):
    def test_transfer_matches_exhaustive(self):
        for pairs in range(1, 7):
            expected = Counter()
            for bits in itertools.product((0, 1), repeat=2 * pairs):
                a = sum(bits[i] * bits[i+1] for i in range(0, 2*pairs, 2))
                b = sum(bits[i] * bits[(i+1) % (2*pairs)] for i in range(1, 2*pairs, 2))
                expected[sum(bits), a, b] += 1
            self.assertEqual(dict(expected), exact_profiles(pairs))

    def test_marginals_runs_and_upper(self):
        for pairs in range(1, 11):
            profile = exact_profiles(pairs)
            even, by_run, by_weight = defaultdict(int), defaultdict(int), defaultdict(int)
            for (r, a, b), count in profile.items():
                self.assertLessEqual(count, profile_count_upper(pairs, r, a, b))
                self.assertEqual(count, profile[r, b, a])
                even[r, a] += count
                by_run[r, r-a-b] += count
                by_weight[r] += count
            for (r, a), count in even.items():
                self.assertEqual(count, one_pairing_count(pairs, r, a))
            for (r, runs), count in by_run.items():
                self.assertEqual(count, circular_run_count(2*pairs, r, runs))
            for r, count in by_weight.items():
                self.assertEqual(count, comb(2*pairs, r))
            self.assertEqual(sum(profile.values()), 1 << (2*pairs))

    def test_constant_and_impossible_profiles(self):
        self.assertEqual(profile_count_upper(4096, 0, 0, 0), 1)
        self.assertEqual(profile_count_upper(4096, 8192, 4096, 4096), 1)
        self.assertEqual(profile_count_upper(4096, 0, 1, 0), 0)
        self.assertEqual(profile_count_upper(4096, 1, 1, 0), 0)
        self.assertEqual(profile_count_upper(4096, 1024, 512, 512), 0)
        self.assertEqual(profile_count_upper(4096, 1024, 512, 511), 4096)


if __name__ == '__main__':
    unittest.main()
