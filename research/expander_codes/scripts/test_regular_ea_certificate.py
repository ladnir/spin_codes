#!/usr/bin/env python3

from __future__ import annotations

import unittest

from flint import arb, ctx

from regular_ea_certificate import (
    exact_regular_ea_tail_arb,
    regional_shell_distribution_arb,
)
from regular_ea_diagnostic import exact_regular_ea_logterm


class RegularEACertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 128

    def test_shell_distribution_is_normalized(self) -> None:
        distribution = regional_shell_distribution_arb(11, 5)
        self.assertTrue(sum(distribution.values(), arb(0)).overlaps(arb(1)))

    def test_exact_tail_is_a_probability_bound(self) -> None:
        selected = exact_regular_ea_logterm(
            k=16,
            region_length=8,
            region_count=4,
            cutoff=3,
            message_weight=2,
        )
        tail = exact_regular_ea_tail_arb(
            region_length=8,
            region_count=4,
            cutoff=3,
            message_weight=2,
            output_marker=format(selected.output_marker, ".17g"),
        )
        self.assertTrue(tail > 0)
        self.assertTrue(tail <= 1)


if __name__ == "__main__":
    unittest.main()
