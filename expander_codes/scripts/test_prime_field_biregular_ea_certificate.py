#!/usr/bin/env python3

import math
import unittest

import numpy as np
from flint import arb, ctx

from expander_bounds import log_binom, log_weight_mgf
from prime_field_biregular_ea_certificate import (
    balanced_occupancy_size_distribution_arb,
    block_markers_float,
    block_term_arb,
    exact_region_trace_matrix_arb,
    exact_marker_schedule_float,
)
from prime_field_regular_ea_trace_diagnostic import (
    balanced_occupancy_size_distribution,
    balanced_slot_trace_matrix,
    uniform_occupancy_trace_transfers,
)


class PrimeFieldBiregularEACertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 160

    def test_balanced_occupancy_arb_matches_float(self) -> None:
        expected = balanced_occupancy_size_distribution(4, 3, 5)
        actual = balanced_occupancy_size_distribution_arb(4, 3, 5)
        for interval, value in zip(actual, expected):
            self.assertAlmostEqual(float(interval), value, places=14)

    def test_exact_region_arb_contains_float(self) -> None:
        length, group_size, r = 4, 3, 3
        z, v = 0.4, 0.3
        occupancy = balanced_occupancy_size_distribution(length, group_size, r)
        slices = uniform_occupancy_trace_transfers(
            region_length=length, max_occupied=r, z=z, v=v
        )
        expected = sum(
            (probability * slices[occupied]
             for occupied, probability in enumerate(occupancy)),
            np.zeros((2, 2)),
        )
        actual = exact_region_trace_matrix_arb(
            region_length=length,
            group_size=group_size,
            message_weight=r,
            z=arb(str(z)),
            v=arb(str(v)),
        )
        for row in range(2):
            for column in range(2):
                self.assertAlmostEqual(
                    float(actual[row, column]), expected[row, column], places=14
                )

    def test_uniform_block_dominates_its_singletons(self) -> None:
        prime, k, n, cutoff, regions = 101, 30, 60, 20, 6
        lo, hi = 3, 5
        markers, _ = block_markers_float(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=regions,
            lo=lo,
            hi=hi,
        )
        block = block_term_arb(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=regions,
            lo=lo,
            hi=hi,
            markers=markers,
        )
        total = 0.0
        group_size = k // (n // regions)
        for r in range(lo, hi + 1):
            for name, field in (("structural", False), ("field", True)):
                x, z, v = map(float, markers[name])
                value = (
                    (1 - regions) * log_binom(k, r)
                    - regions * r * math.log(x)
                    + log_weight_mgf(
                        balanced_slot_trace_matrix(group_size, x, z, v), n, 0
                    )
                    - cutoff * math.log(z)
                    - (r - (0 if field else 1)) * math.log(v)
                    - (math.log(prime - 1) if field else 0.0)
                )
                total += math.exp(value)
        self.assertGreater(float(block), total)

    def test_exact_marker_schedule_covers_every_weight(self) -> None:
        schedule = exact_marker_schedule_float(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            limit=5,
        )
        self.assertEqual([item["r"] for item in schedule], [1, 2, 3, 4, 5])


if __name__ == "__main__":
    unittest.main()
