import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from flint import arb, ctx

from prime_field_biregular_ec_d22_m12_certificate import (
    BLOCKS,
    CUTOFF,
    DEGREE,
    K,
    MEMORY,
    N,
    REGION_LENGTH,
    SPARSE_LIMIT,
    VerificationCache,
    sparse_structural_term_arb,
)
from prime_field_biregular_ec_certificate import VerificationResult
from prime_field_biregular_ec_diagnostic import log_binom
from prime_field_regular_ea_trace_diagnostic import (
    balanced_occupancy_size_distribution,
)


class Degree22Memory12CertificateTests(unittest.TestCase):
    def test_support_blocks_are_contiguous(self) -> None:
        expected = 1401
        for lo, hi in BLOCKS:
            self.assertEqual(lo, expected)
            self.assertGreaterEqual(hi, lo)
            expected = hi + 1
        self.assertEqual(expected, K)

    def test_sparse_structural_arb_matches_float_diagnostic(self) -> None:
        ctx.prec = 160
        actual = sparse_structural_term_arb(
            k=K,
            n=N,
            cutoff=CUTOFF,
            region_count=DEGREE,
            memory=MEMORY,
            support_limit=SPARSE_LIMIT,
        )
        logs = []
        for r in range(1, SPARSE_LIMIT + 1):
            law = balanced_occupancy_size_distribution(
                REGION_LENGTH, K // REGION_LENGTH, r
            )
            bad_probability = sum(law[:r // 2 + 1])
            zero_intervals = 1 + (r - 1) // MEMORY
            needed_regions = (
                (N - CUTOFF - (r - 1) + REGION_LENGTH - 1)
                // REGION_LENGTH
                - 2 * zero_intervals
            )
            if bad_probability:
                logs.append(
                    log_binom(K, r)
                    + log_binom(DEGREE, needed_regions)
                    + needed_regions * math.log(bad_probability)
                )
        expected_log2 = math.log2(sum(math.exp(x - max(logs)) for x in logs))
        expected_log2 += max(logs) / math.log(2.0)
        actual_log2 = float(actual.log() / type(actual)(2).log())
        self.assertAlmostEqual(actual_log2, expected_log2, places=8)

    def test_verification_cache_is_content_addressed(self) -> None:
        certificate = {"verification": {"precision_bits": 128}, "value": 7}
        with TemporaryDirectory() as directory:
            cache = VerificationCache(
                certificate=certificate, directory=Path(directory)
            )
            self.assertIsNone(cache.get_bound("exact-1-2"))
            cache.put_bound("exact-1-2", arb("0.125"))
            self.assertEqual(cache.get_bound("exact-1-2"), arb("0.125"))
            other = VerificationCache(
                certificate={**certificate, "value": 8},
                directory=Path(directory),
            )
            self.assertIsNone(other.get_bound("exact-1-2"))
            refresh = VerificationCache(
                certificate=certificate,
                directory=Path(directory),
                refresh=True,
            )
            self.assertIsNone(refresh.get_bound("exact-1-2"))

    def test_verification_cache_round_trips_complete_result(self) -> None:
        certificate = {"verification": {"precision_bits": 128}}
        result = VerificationResult(
            success=True,
            total_bound=arb("0.01"),
            security_bits=arb("6.5"),
            exact_bound=arb("0.009"),
            block_bound=arb("0.0009"),
            full_support_bound=arb("0.0001"),
            largest_exact_band=(3, 4),
            largest_exact_term=arb("0.008"),
            largest_block=(5, 6),
            largest_block_term=arb("0.0008"),
        )
        with TemporaryDirectory() as directory:
            cache = VerificationCache(
                certificate=certificate, directory=Path(directory)
            )
            cache.put_result(result)
            restored = cache.get_result()
            self.assertIsNotNone(restored)
            self.assertTrue(restored.total_bound.contains(result.total_bound))
            self.assertEqual(restored.largest_exact_band, (3, 4))


if __name__ == "__main__":
    unittest.main()
