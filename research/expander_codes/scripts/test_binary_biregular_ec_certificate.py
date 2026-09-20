#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest
from copy import deepcopy
from unittest import mock

from flint import arb, ctx

import binary_biregular_ec_certificate as certificate_module
from binary_biregular_ec_certificate import (
    coefficient_block_arb,
    dense_low_half_block_arb,
    exact_term_arb,
    full_support_term_arb,
    odd_degree_variance_lower_arb,
    odd_half_count_variance_arb,
    verify_certificate,
)
from check_binary_biregular_ec_certificate import (
    ODD_DEGREE_SCHEMA,
    SCHEMA,
    validate_certificate_structure,
)


def small_certificate() -> dict:
    return {
        "schema": SCHEMA,
        "parameters": {
            "k": 25,
            "n": 50,
            "cutoff": 4,
            "left_degree": 10,
            "right_degree": 5,
            "memory": 3,
            "target_bits": 0,
        },
        "verification": {"precision_bits": 128},
        "exact_weights": [{"r": 1, "z": "0.7"}],
        "outer_blocks": [
            {"lo": 2, "hi": 4, "x": "0.3", "z": "0.7"},
            {"lo": 21, "hi": 24, "x": "3.0", "z": "0.7"},
        ],
        "central_blocks": [
            {"lo": 5, "hi": 12, "include_complement": True},
        ],
        "full_support": {"r": 25, "z": "0.7"},
    }


class BinaryBiregularECCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 128

    def test_odd_variance_at_half(self) -> None:
        self.assertTrue(odd_half_count_variance_arb(arb(1) / 2).contains(arb(5) / 16))

    def test_generic_variance_lower_bound(self) -> None:
        for degree in (3, 5, 7, 9):
            p_lo = arb(1) / degree
            lower = odd_degree_variance_lower_arb(p_lo, degree)
            self.assertTrue(lower > 0)
            for p in (p_lo, arb(1) / 4, arb(2) / 5, arb(1) / 2):
                if p < p_lo:
                    continue
                for parity in (0, 1):
                    masses = [
                        arb(math.comb(degree, count))
                        * p**count * (1 - p) ** (degree - count)
                        for count in range(parity, degree + 1, 2)
                    ]
                    total = sum(masses, arb(0))
                    mean = sum(
                        mass * index for index, mass in enumerate(masses)
                    ) / total
                    variance = sum(
                        mass * (index - mean) ** 2
                        for index, mass in enumerate(masses)
                    ) / total
                    self.assertFalse(lower > variance)

    def test_exact_term_is_positive(self) -> None:
        value = exact_term_arb(
            k=15,
            left_degree=10,
            right_degree=5,
            cutoff=4,
            memory=3,
            message_weight=2,
            output_marker="0.7",
        )
        self.assertTrue(value > 0)

    def test_coefficient_block_is_positive(self) -> None:
        value = coefficient_block_arb(
            k=15,
            left_degree=10,
            right_degree=5,
            cutoff=4,
            memory=3,
            lo=2,
            hi=3,
            input_marker="0.2",
            output_marker="0.7",
        )
        self.assertTrue(value > 0)

    def test_dense_block_is_positive(self) -> None:
        value = dense_low_half_block_arb(
            k=25,
            left_degree=10,
            cutoff=4,
            lo=5,
            hi=8,
        )
        self.assertTrue(value > 0)

    def test_generic_odd_degree_dense_block_is_positive(self) -> None:
        value = dense_low_half_block_arb(
            k=3507,
            left_degree=14,
            right_degree=7,
            cutoff=500,
            lo=501,
            hi=800,
        )
        self.assertTrue(value > 0)

    def test_checker_accepts_generic_odd_degree(self) -> None:
        certificate = small_certificate()
        certificate["schema"] = ODD_DEGREE_SCHEMA
        certificate["parameters"].update({
            "k": 35,
            "n": 70,
            "cutoff": 5,
            "left_degree": 14,
            "right_degree": 7,
        })
        certificate["outer_blocks"] = [
            {"lo": 2, "hi": 4, "x": "0.3", "z": "0.7"},
            {"lo": 31, "hi": 34, "x": "3.0", "z": "0.7"},
        ]
        certificate["central_blocks"] = [
            {"lo": 5, "hi": 17, "include_complement": True},
        ]
        certificate["full_support"] = {"r": 35, "z": "0.7"}
        validate_certificate_structure(certificate)

    def test_full_support_term_is_positive(self) -> None:
        value = full_support_term_arb(
            n=30,
            cutoff=4,
            memory=3,
            output_marker="0.7",
        )
        self.assertTrue(value > 0)

    def test_independent_checker_accepts_complete_partition(self) -> None:
        validate_certificate_structure(small_certificate())

    def test_independent_checker_rejects_support_gap(self) -> None:
        certificate = deepcopy(small_certificate())
        certificate["outer_blocks"][0]["lo"] = 3
        with self.assertRaisesRegex(ValueError, "gap"):
            validate_certificate_structure(certificate)

    def test_independent_checker_rejects_nondecimal_marker(self) -> None:
        certificate = deepcopy(small_certificate())
        certificate["exact_weights"][0]["z"] = "NaN"
        with self.assertRaisesRegex(ValueError, "finite and positive"):
            validate_certificate_structure(certificate)

    def test_independent_checker_rejects_reordered_blocks(self) -> None:
        certificate = deepcopy(small_certificate())
        certificate["outer_blocks"].reverse()
        with self.assertRaisesRegex(ValueError, "ordered and disjoint"):
            validate_certificate_structure(certificate)

    def test_verification_uses_only_frozen_markers(self) -> None:
        with (
            mock.patch.object(
                certificate_module,
                "exact_biregular_logterm",
                side_effect=AssertionError("optimizer was called"),
            ),
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
