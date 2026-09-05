#!/usr/bin/env python3

import math
import unittest

from flint import arb, ctx

from expander_bounds import log_binom, log_weight_mgf
from prime_field_biregular_ec_certificate import (
    DEFAULT_BLOCKS,
    DEFAULT_EXACT_BANDS,
    DEFAULT_K,
    balanced_occupancy_prefix_arb,
    block_term_arb,
    exact_band_term_arb,
    nonwrapping_trace_matrices_arb,
    singleton_exact_band_term_arb,
    singleton_block_term_arb,
    singleton_region_polynomial_arb,
    validate_coverage,
)
from prime_field_biregular_ec_diagnostic import (
    balanced_constraint_slot_trace_matrix,
    balanced_occupancy_prefix,
    constraint_balanced_ec_block_logbound,
    nonwrapping_constraint_trace_matrices,
    singleton_constraint_uniform_support_trace_transfers,
    singleton_constraint_balanced_ec_block_logbound,
)


class PrimeFieldBiregularECCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 160

    def test_nonwrapping_arb_matrix_matches_float(self) -> None:
        expected = nonwrapping_constraint_trace_matrices(
            memory=2, output_marker=0.4, equation_marker=0.3
        )
        actual = nonwrapping_trace_matrices_arb(
            memory=2,
            output_marker=arb("0.4"),
            equation_marker=arb("0.3"),
        )
        for expected_matrix, actual_matrix in zip(expected, actual):
            for row in range(3):
                for column in range(3):
                    self.assertAlmostEqual(
                        float(actual_matrix[row, column]),
                        expected_matrix[row, column],
                        places=14,
                    )

    def test_occupancy_prefix_arb_matches_float(self) -> None:
        expected = balanced_occupancy_prefix(
            region_length=5, group_size=2, limit=4
        )
        actual = balanced_occupancy_prefix_arb(
            region_length=5, group_size=2, limit=4
        )
        for expected_law, actual_law in zip(expected, actual):
            for expected_value, actual_value in zip(expected_law, actual_law):
                self.assertAlmostEqual(
                    float(actual_value), expected_value, places=14
                )

    def test_exact_band_term_is_positive(self) -> None:
        term = exact_band_term_arb(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            lo=1,
            hi=3,
            markers={
                "structural": ["0.7", "0.2"],
                "field": ["0.7", "0.2"],
            },
        )
        self.assertTrue(term > 0)

    def test_singleton_region_polynomial_matches_normalized_float_slices(self) -> None:
        region_length, group_size, maximum = 3, 2, 4
        polynomial = singleton_region_polynomial_arb(
            region_length=region_length,
            group_size=group_size,
            max_weight=maximum,
            memory=1,
            output_marker=arb("0.6"),
            equation_marker=arb("0.2"),
        )
        slices = singleton_constraint_uniform_support_trace_transfers(
            region_length=region_length,
            group_size=group_size,
            max_weight=maximum,
            memory=1,
            output_marker=0.6,
            equation_marker=0.2,
        )
        for weight in range(maximum + 1):
            normalizer = math.comb(region_length * group_size, weight)
            for row in range(2):
                for column in range(2):
                    self.assertAlmostEqual(
                        float(polynomial[row][column][weight] / normalizer),
                        slices[weight][row, column],
                        places=13,
                    )

    def test_singleton_exact_band_term_is_positive(self) -> None:
        term = singleton_exact_band_term_arb(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            lo=1,
            hi=3,
            markers={
                "structural": ["0.7", "0.2"],
                "field": ["0.7", "0.2"],
            },
        )
        self.assertTrue(term > 0)

    def test_uniform_block_dominates_singletons(self) -> None:
        prime, k, n, cutoff, regions, memory = 101, 30, 60, 20, 6, 1
        lo, hi = 3, 5
        bound = constraint_balanced_ec_block_logbound(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=regions,
            memory=memory,
            support_start=lo,
            support_limit=hi,
        )
        markers = {
            "structural": [
                format(value, ".17g") for value in bound.structural_markers
            ],
            "field": [
                format(value, ".17g") for value in bound.field_markers
            ],
        }
        block = block_term_arb(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=regions,
            memory=memory,
            lo=lo,
            hi=hi,
            markers=markers,
        )
        total = 0.0
        group_size = k // (n // regions)
        for r in range(lo, hi + 1):
            for name, field in (("structural", False), ("field", True)):
                z, v, x = map(float, markers[name])
                value = (
                    (1 - regions) * log_binom(k, r)
                    - regions * r * math.log(x)
                    + log_weight_mgf(
                        balanced_constraint_slot_trace_matrix(
                            memory=memory,
                            group_size=group_size,
                            slot_marker=x,
                            output_marker=z,
                            equation_marker=v,
                        ),
                        n,
                        memory,
                    )
                    - cutoff * math.log(z)
                    - (r - (0 if field else 1)) * math.log(v)
                    - (math.log(prime - 1) if field else 0.0)
                )
                total += math.exp(value)
        self.assertGreater(float(block), total)

    def test_singleton_block_arb_matches_float_diagnostic(self) -> None:
        parameters = dict(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            memory=2,
            support_start=3,
            support_limit=5,
        )
        bound = singleton_constraint_balanced_ec_block_logbound(**parameters)
        markers = {
            "structural": [
                format(value, ".17g") for value in bound.structural_markers
            ],
            "field": [
                format(value, ".17g") for value in bound.field_markers
            ],
        }
        interval = singleton_block_term_arb(
            prime=parameters["prime"],
            k=parameters["k"],
            n=parameters["n"],
            cutoff=parameters["cutoff"],
            region_count=parameters["region_count"],
            memory=parameters["memory"],
            lo=parameters["support_start"],
            hi=parameters["support_limit"],
            markers=markers,
        )
        certified_log2 = float(interval.log() / arb(2).log())
        # The interval verifier uses the sum of the two endpoints; the float
        # selector uses their maximum.  This costs at most one bit.
        self.assertGreaterEqual(certified_log2, bound.total_log2 - 1e-8)
        self.assertLessEqual(certified_log2, bound.total_log2 + 1.0 + 1e-8)

    def test_default_partition_covers_every_support(self) -> None:
        certificate = {
            "parameters": {"k": DEFAULT_K},
            "exact_bands": [
                {"lo": lo, "hi": hi} for lo, hi in DEFAULT_EXACT_BANDS
            ],
            "blocks": [{"lo": lo, "hi": hi} for lo, hi in DEFAULT_BLOCKS],
            "full_support": {"r": DEFAULT_K},
        }
        validate_coverage(certificate)


if __name__ == "__main__":
    unittest.main()
