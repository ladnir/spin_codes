import unittest
from math import log2
import numpy as np

import hd_field128_gamma_profile_diagnostic as diagnostic
from hd_shifted_pairs_count import exact_profiles, profile_count_upper


class ProfileTests(unittest.TestCase):
    def test_every_small_count_matches_existing_integer_bound(self):
        for pairs in range(1, 10):
            actual = exact_profiles(pairs)
            for weight in range(2*pairs+1):
                full = diagnostic.transfer.tail.full_counts(pairs, weight)
                logs = diagnostic.profile_log_counts(pairs, weight, full[:, None], full[None])
                for i, first in enumerate(full):
                    for j, second in enumerate(full):
                        integer = profile_count_upper(pairs, weight, int(first), int(second))
                        self.assertGreaterEqual(integer, actual.get((weight, int(first), int(second)), 0))
                        if integer:
                            self.assertAlmostEqual(logs[i, j], log2(integer), places=10)
                        else:
                            self.assertEqual(logs[i, j], -np.inf)

    def test_sharpened_sum_is_no_larger_than_run_only(self):
        rng = np.random.default_rng(87128)
        for weight in (1, 2, 3, 4, 7, 8):
            full = diagnostic.transfer.tail.full_counts(4, weight)
            regions = -rng.uniform(0, 20, (len(full), 2, 2))
            old = diagnostic.transfer.support_sum(4, weight, regions, repetitions=8)
            new = diagnostic.support_sum(4, weight, regions, repetitions=8)
            self.assertLessEqual(new['log2_support_sum_screen'], old['log2_support_sum_screen']+1e-10)


if __name__ == '__main__':
    unittest.main()
