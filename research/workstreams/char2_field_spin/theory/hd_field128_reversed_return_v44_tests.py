"""Independent tiny-field rank oracle inherited, plus reversal checks."""
import unittest
import hd_field128_dense_feedback_return_v42_tests as oracle
import hd_field128_reversed_return_v44 as candidate
import hd_field128_dense_gamma_exclusions_v43 as exclusions


class Tests(oracle.Tests):
    def test_toy_actual_ranks(self):
        # Substitute only the candidate API in the existing independent
        # actual-matrix rank oracle; restore it even on an assertion failure.
        previous = oracle.candidate
        try:
            oracle.candidate = candidate
            super().test_toy_actual_ranks()
        finally:
            oracle.candidate = previous

    def test_all_count_reversals(self):
        for j in range(65):
            for c in range(j+1):
                for e in range(65-j):
                    t, v = j-c, 64-j-e
                    new = candidate.rank_floor(j, c, e)
                    self.assertEqual(new, candidate.rank_floor(t+v, v, e))
                    self.assertGreaterEqual(new, candidate.prior.rank_floor(j, c, e))
        for row in exclusions.exclusion_count(58)['shapes']:
            j, c, e = row[:3]
            self.assertEqual(candidate.probability(j, c, e, 1 << 128), 0)


if __name__ == '__main__':
    unittest.main()
