#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest

from flint import arb, ctx

from exact_row_ec_certificate import (
    wrapping_enumerator_matrix_arb,
    wrapping_slice_tail_bound_arb,
)
from expander_bounds import ec_enumerator_matrix, log_weight_mgf


class ExactRowECCertificateTests(unittest.TestCase):
    def test_arb_wrapping_matrix_contains_float_evaluation(self) -> None:
        ctx.prec = 128
        n, m, x_value, z_value = 7, 3, 0.031, 0.73
        x, z = arb("0.031"), arb("0.73")
        power = wrapping_enumerator_matrix_arb(x, z, m) ** n
        arb_value = sum((power[m, j] for j in range(m + 1)), arb(0))
        float_value = math.exp(
            log_weight_mgf(ec_enumerator_matrix(x_value, z_value, m, True), n, m)
        )
        # The Arb path encloses the exact decimal markers, whereas binary64
        # evaluates nearby dyadic values. Their midpoints should nevertheless
        # agree to ordinary floating-point accuracy.
        self.assertAlmostEqual(float(arb_value), float_value, places=13)

    def test_slice_bound_accepts_fixed_decimal_markers(self) -> None:
        ctx.prec = 128
        bound = wrapping_slice_tail_bound_arb(
            n=7,
            cutoff=2,
            shell_weight=2,
            memory=2,
            input_marker="0.2",
            output_marker="0.7",
        )
        self.assertTrue(bound > 0)


if __name__ == "__main__":
    unittest.main()
