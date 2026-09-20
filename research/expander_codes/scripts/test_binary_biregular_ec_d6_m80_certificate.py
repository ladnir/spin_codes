#!/usr/bin/env python3

from __future__ import annotations

import unittest
from copy import deepcopy
from unittest import mock

import numpy as np
from flint import arb, ctx

import binary_biregular_ec_d6_m80_certificate as certificate_module
from binary_biregular_ec_d6_m80_certificate import (
    degree_three_dense_low_half_block_arb,
    degree_three_odd_variance_arb,
    exact_block_arb,
    exact_block_factor_upper,
    exact_block_layered_upper,
    exact_block_rows_upper,
    exact_block_upper,
    verify_certificate,
)
from check_binary_biregular_ec_d6_m80_certificate import (
    SCHEMA,
    validate_certificate_structure,
)
from positive_matrix_upper import matmul_upper


def small_certificate() -> dict:
    return {
        "schema": SCHEMA,
        "parameters": {
            "k": 303,
            "n": 606,
            "cutoff": 60,
            "left_degree": 6,
            "right_degree": 3,
            "memory": 3,
            "target_bits": 0,
        },
        "verification": {"precision_bits": 128},
        "exact_blocks": [{"lo": 1, "hi": 2, "z": "0.9"}],
        "outer_blocks": [
            {"lo": 3, "hi": 119, "x": "0.3", "z": "0.7"},
            {"lo": 184, "hi": 302, "x": "3.0", "z": "0.7"},
        ],
        "central_blocks": [
            {"lo": 120, "hi": 151, "include_complement": True},
        ],
        "full_support": {"r": 303, "z": "0.7"},
    }


class DegreeSixECCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 128

    def test_degree_three_variance_at_half(self) -> None:
        value = degree_three_odd_variance_arb(arb(1) / 2)
        self.assertTrue(value.contains(arb(3) / 16))

    def test_directed_product_preserves_structural_zero(self) -> None:
        left = np.array([[1.0, 0.0], [0.0, 1.0]])
        right = np.array([[1.0, 0.0], [1.0, 0.0]])
        upper = matmul_upper(left, right)
        self.assertEqual(upper[0, 1], 0.0)
        self.assertEqual(upper[1, 1], 0.0)

    def test_exact_block_is_positive(self) -> None:
        value = exact_block_arb(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            lo=2,
            hi=4,
            output_marker="0.8",
        )
        self.assertTrue(value > 0)

    def test_directed_block_dominates_direct_arb(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            lo=2,
            hi=4,
            output_marker="0.8",
        )
        exact = exact_block_arb(**parameters)
        upper = exact_block_upper(**parameters)
        self.assertTrue(upper > exact)
        self.assertTrue(upper / exact < arb("1.000000001"))

    def test_factor_block_dominates_direct_arb(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            lo=2,
            hi=4,
            output_marker="0.8",
        )
        exact = exact_block_arb(**parameters)
        upper = exact_block_factor_upper(**parameters)
        self.assertTrue(upper > exact)
        self.assertTrue(upper / exact < arb("1.000000001"))

    def test_row_scaled_block_dominates_direct_arb(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            lo=2,
            hi=4,
            output_marker="0.8",
        )
        exact = exact_block_arb(**parameters)
        upper = exact_block_rows_upper(**parameters)
        self.assertTrue(upper > exact)
        self.assertTrue(upper / exact < arb("1.000000001"))

    def test_layered_block_dominates_direct_arb(self) -> None:
        parameters = dict(
            k=15,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            lo=2,
            hi=4,
            output_marker="0.8",
        )
        exact = exact_block_arb(**parameters)
        upper = exact_block_layered_upper(**parameters)
        self.assertTrue(upper > exact)
        self.assertTrue(upper / exact < arb("1.000000001"))

    def test_degree_three_dense_block_is_positive(self) -> None:
        value = degree_three_dense_low_half_block_arb(
            k=303,
            left_degree=6,
            cutoff=60,
            lo=120,
            hi=130,
        )
        self.assertTrue(value > 0)

    def test_checker_accepts_complete_partition(self) -> None:
        validate_certificate_structure(small_certificate())

    def test_checker_rejects_gap(self) -> None:
        certificate = deepcopy(small_certificate())
        certificate["outer_blocks"][0]["lo"] = 4
        with self.assertRaisesRegex(ValueError, "gap"):
            validate_certificate_structure(certificate)

    def test_checker_rejects_wrong_degrees(self) -> None:
        certificate = deepcopy(small_certificate())
        certificate["parameters"]["left_degree"] = 10
        with self.assertRaisesRegex(ValueError, "degrees 6/3"):
            validate_certificate_structure(certificate)

    def test_verification_does_not_call_marker_optimizers(self) -> None:
        with (
            mock.patch.object(
                certificate_module,
                "saddle_biregular_block_logterm",
                side_effect=AssertionError("optimizer was called"),
            ),
            mock.patch.object(
                certificate_module,
                "full_support_marker",
                side_effect=AssertionError("optimizer was called"),
            ),
        ):
            result = verify_certificate(small_certificate())
        self.assertTrue(result.total_bound > 0)


if __name__ == "__main__":
    unittest.main()
