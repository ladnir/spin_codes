#!/usr/bin/env python3
"""Tests for scaled binary biregular EC diagnostics."""

from __future__ import annotations

import math
import unittest

import numpy as np

from binary_biregular_diagnostic import (
    _exact_region_matrix,
    biregular_parity_shell_distribution,
)
from binary_biregular_hybrid_diagnostic import (
    scaled_exact_fixed_marker_block,
    scaled_exact_region_logs,
    scaled_uniform_slice_transfers,
)
from regular_ec_diagnostic import uniform_slice_transfer_matrices


class ScaledUniformSliceTest(unittest.TestCase):
    def test_scaled_slices_match_ordinary_arithmetic(self) -> None:
        ordinary = uniform_slice_transfer_matrices(
            length=11,
            max_weight=5,
            output_marker=0.73,
            memory=3,
        )
        scaled = scaled_uniform_slice_transfers(
            length=11,
            max_weight=5,
            output_marker=0.73,
            memory=3,
        )
        for expected, (matrix, log_scale) in zip(ordinary, scaled, strict=True):
            np.testing.assert_allclose(
                math.exp(log_scale) * matrix,
                expected,
                rtol=2e-13,
                atol=2e-15,
            )

    def test_scaled_shell_mixture_matches_ordinary_arithmetic(self) -> None:
        k = 15
        right_degree = 3
        support = 5
        memory = 3
        marker = 0.81
        shell = biregular_parity_shell_distribution(
            left_vertices=k,
            right_degree=right_degree,
            support_size=support,
        )
        ordinary = _exact_region_matrix(
            code="ec",
            region_length=k // right_degree,
            support_size=support,
            shell=shell,
            output_marker=marker,
            memory=memory,
        )
        logs = scaled_exact_region_logs(
            k=k,
            right_degree=right_degree,
            memory=memory,
            message_weight=support,
            output_marker=marker,
        )
        np.testing.assert_allclose(
            np.exp(logs), ordinary, rtol=3e-13, atol=2e-15
        )

    def test_scaled_block_matches_individual_region_evaluation(self) -> None:
        summed, worst_support, worst, terms = scaled_exact_fixed_marker_block(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            support_start=2,
            support_limit=4,
            output_marker=0.8,
        )
        self.assertEqual(len(terms), 3)
        self.assertEqual(worst_support, 2 + int(np.argmax(terms)))
        self.assertEqual(worst, max(terms))
        self.assertGreaterEqual(summed, worst)


if __name__ == "__main__":
    unittest.main()
