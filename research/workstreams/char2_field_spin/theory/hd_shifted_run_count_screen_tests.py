import unittest

import numpy as np

from hd_shifted_run_count_screen import log_run_counts_float, run_log_screen_float
from hd_shifted_run_potential import log_run_counts, run_log_screen
from hd_shifted_swap_potential import full_counts


class RunCountScreenTests(unittest.TestCase):
    def test_all_small_weights_against_exact_integer_counts(self):
        for pairs in range(1, 13):
            for weight in range(2*pairs+1):
                actual = log_run_counts_float(pairs, weight)
                expected = log_run_counts(pairs, weight)
                np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)
                self.assertFalse(actual.flags.writeable)

    def test_production_count_tables_against_exact_integer_logs(self):
        for weight in (1280, 1536, 2048, 3072):
            actual = log_run_counts_float(4096, weight)
            expected = log_run_counts(4096, weight)
            np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-9)

    def test_noncommuting_regional_fixtures_preserve_screen(self):
        for weight in range(1, 13):
            regions = {f: np.log2([[2+f, 1+2*f], [3+weight-f, 4+f]])
                       for f in full_counts(8, weight)}
            for ratio in (-7., 0., 2.):
                exact = run_log_screen(8, weight, regions, ratio, repetitions=8)
                actual = run_log_screen_float(8, weight, regions, ratio, repetitions=8)
                self.assertAlmostEqual(actual, exact, places=10)

    def test_guards(self):
        for pairs, weight in ((0, 0), (4, -1), (4, 9), (4.5, 1)):
            with self.assertRaises(ValueError):
                log_run_counts_float(pairs, weight)
        with self.assertRaises(ValueError):
            run_log_screen_float(4, 2, {}, 0, repetitions=3)


if __name__ == '__main__':
    unittest.main()
