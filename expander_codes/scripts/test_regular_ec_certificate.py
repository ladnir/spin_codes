#!/usr/bin/env python3

from __future__ import annotations

import unittest

from flint import arb, ctx

from regular_ec_certificate import (
    exact_regular_tail_arb,
    poisson_penalty_arb,
    regional_shell_distribution_arb,
)
from regular_ec_diagnostic import exact_regular_ec_logterm


class RegularECCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 128

    def test_shell_distribution_is_normalized(self) -> None:
        distribution = regional_shell_distribution_arb(11, 5)
        self.assertTrue(sum(distribution.values(), arb(0)).overlaps(arb(1)))

    def test_poisson_penalty_is_greater_than_one(self) -> None:
        for r in range(1, 8):
            self.assertTrue(poisson_penalty_arb(r) > 1)

    def test_exact_arb_tail_contains_float_bound(self) -> None:
        selected = exact_regular_ec_logterm(
            k=16,
            region_length=8,
            region_count=4,
            cutoff=3,
            memory=3,
            message_weight=2,
        )
        tail = exact_regular_tail_arb(
            region_length=8,
            region_count=4,
            cutoff=3,
            memory=3,
            message_weight=2,
            output_marker=format(selected.output_marker, ".17g"),
        )
        self.assertTrue(tail > 0)
        self.assertTrue(tail <= 1)


if __name__ == "__main__":
    unittest.main()
