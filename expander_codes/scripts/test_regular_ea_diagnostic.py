#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest

import numpy as np

from regular_ea_diagnostic import (
    accumulator_input_matrices,
    accumulator_uniform_slice_transfers,
    exact_regular_ea_logterm,
)


class RegularEADiagnosticTests(unittest.TestCase):
    def test_accumulator_slices_match_direct_products(self) -> None:
        zero, one = accumulator_input_matrices(0.7)
        slices = accumulator_uniform_slice_transfers(
            length=2, max_weight=2, output_marker=0.7
        )
        np.testing.assert_allclose(slices[0], zero @ zero)
        np.testing.assert_allclose(slices[1], (zero @ one + one @ zero) / 2.0)
        np.testing.assert_allclose(slices[2], one @ one)

    def test_exact_small_term_is_finite(self) -> None:
        term = exact_regular_ea_logterm(
            k=16,
            region_length=8,
            region_count=4,
            cutoff=3,
            message_weight=2,
        )
        self.assertTrue(math.isfinite(term.log2_bound))


if __name__ == "__main__":
    unittest.main()
