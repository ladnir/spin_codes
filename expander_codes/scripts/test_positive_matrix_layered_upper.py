#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest

import numpy as np
from flint import arb, arb_mat, ctx

from positive_matrix_layered_upper import (
    layered_from_matrix,
    multiply_layered,
    multiply_scale_layered_rational,
    power_layered,
    rational_mantissa_exponent_upper,
    row_sum_layers,
    scale_layered_rational,
    sum_layered,
)


def arb_from_binary64(value: float) -> arb:
    numerator, denominator = value.as_integer_ratio()
    return arb(numerator) / denominator


def layered_entry(value, row: int, column: int) -> arb:
    return sum((
        arb_from_binary64(float(layer.matrix[row, column]))
        * arb(2) ** layer.exponent
        for layer in value.layers
    ), arb(0))


class LayeredUpperMatrixTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 160

    def test_independent_entry_scales_survive(self) -> None:
        high = layered_from_matrix(np.array([[1.0, 0.0], [0.0, 1.0]]))
        low = layered_from_matrix(
            np.array([[0.0, 0.75], [0.0, 0.0]]), exponent=-1600
        )
        total = sum_layered([high, low])
        self.assertEqual(len(total.layers), 2)
        self.assertTrue(layered_entry(total, 0, 1) > 0)
        self.assertTrue(layered_entry(total, 0, 1) < arb(2) ** -1599)

    def test_rational_scale_does_not_underflow(self) -> None:
        numerator = 1
        denominator = 1 << 2000
        mantissa, exponent = rational_mantissa_exponent_upper(
            numerator, denominator
        )
        self.assertGreaterEqual(mantissa, 0.5)
        self.assertEqual(exponent, -1999)
        identity = layered_from_matrix(np.eye(2))
        scaled = scale_layered_rational(identity, numerator, denominator)
        self.assertTrue(layered_entry(scaled, 0, 0) >= arb(2) ** -2000)
        self.assertTrue(layered_entry(scaled, 0, 0) < arb(2) ** -1999)

    def test_product_dominates_exact_arb(self) -> None:
        left_array = np.array(((0.75, 0.125), (0.5, 0.25)))
        right_array = np.array(((0.2, 0.4), (0.3, 0.1)))
        left = layered_from_matrix(left_array)
        right = layered_from_matrix(right_array)
        product = multiply_layered(left, right)
        exact = arb_mat(left_array.tolist()) * arb_mat(right_array.tolist())
        for row in range(2):
            for column in range(2):
                self.assertTrue(
                    layered_entry(product, row, column)
                    >= exact[row, column]
                )

    def test_fused_product_scale_dominates_exact_arb(self) -> None:
        left_array = np.array(((0.75, 0.125), (0.5, 0.25)))
        right_array = np.array(((0.2, 0.4), (0.3, 0.1)))
        product = multiply_scale_layered_rational(
            layered_from_matrix(left_array),
            layered_from_matrix(right_array),
            1,
            1 << 1400,
        )
        exact = (arb_mat(left_array.tolist()) * arb_mat(right_array.tolist()))
        exact *= arb(2) ** -1400
        for row in range(2):
            for column in range(2):
                self.assertTrue(
                    layered_entry(product, row, column)
                    >= exact[row, column]
                )

    def test_power_and_row_sum_dominate_exact_arb(self) -> None:
        matrix_array = np.array(((0.8, 0.1), (0.25, 0.6)))
        matrix = layered_from_matrix(matrix_array)
        powered = power_layered(matrix, 7)
        exact = arb_mat(matrix_array.tolist()) ** 7
        components = row_sum_layers(powered, 1)
        upper = sum((
            arb_from_binary64(mantissa) * arb(2) ** exponent
            for mantissa, exponent in components
        ), arb(0))
        exact_sum = exact[1, 0] + exact[1, 1]
        self.assertTrue(upper >= exact_sum)
        self.assertTrue(upper / exact_sum < arb("1.0000000001"))


if __name__ == "__main__":
    unittest.main()
