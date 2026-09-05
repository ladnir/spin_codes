#!/usr/bin/env python3

from __future__ import annotations

import unittest

from binary_biregular_activation_diagnostic import (
    activation_band_split_logterm,
    activation_bands,
    exact_logterm,
)


class ActivationBandDiagnosticTests(unittest.TestCase):
    def test_bands_cover_each_feasible_count_once(self) -> None:
        bands = activation_bands(10, 4)
        self.assertEqual(bands, [(0, 0), (1, 4), (5, 8), (9, 10)])

    def test_band_bound_dominates_exact_term(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=3,
            output_marker=0.8,
        )
        exact = exact_logterm(**parameters)
        bound, bands, _ = activation_band_split_logterm(
            **parameters,
            band_width=2,
            log_input_min=-5.0,
            log_input_max=2.0,
            input_steps=15,
            log_activation_limit=5.0,
            activation_steps=15,
        )
        self.assertEqual(bands, [(0, 0), (1, 2), (3, 3)])
        self.assertGreaterEqual(bound + 1e-10, exact)


if __name__ == "__main__":
    unittest.main()
