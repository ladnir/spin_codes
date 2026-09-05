#!/usr/bin/env python3
"""Evaluate the frozen random-linear 512-bit proof-gym model.

Each data-group position receives an independent uniform injective linear map
from 256 to 512 bits. The maps are sampled during setup and then frozen.
Independence makes the exact ensemble first moment factor across data groups.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln

import analyze_riffle_rm512_randomstepconv_g4_sigma20 as base
from analyze_riffle_randomstepconv_g4_sigma20_goal03 import (
    DEFAULT_INNER_RECEIPT,
    load_inner_caps,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/"
    "receipts/frozenrandom512_p2_distance09.json"
)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * base.LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


class RandomLinearMoments(base.ExactRMMoments):
    """Exact ensemble-average spectrum and exact packetization law."""

    def __init__(self) -> None:
        log_nonzero_ratio = (
            log_two_power_minus_one(base.LOCAL_INPUT_BITS)
            - log_two_power_minus_one(base.LOCAL_BITS)
        )

        self.weights = np.arange(base.LOCAL_BITS + 1, dtype=np.int64)
        self.log_multiplicities = np.full(
            base.LOCAL_BITS + 1, -math.inf, dtype=np.float64
        )
        self.log_multiplicities[0] = 0.0
        for weight in range(1, base.LOCAL_BITS + 1):
            self.log_multiplicities[weight] = (
                log_nonzero_ratio
                + float(gammaln(base.LOCAL_BITS + 1))
                - float(gammaln(weight + 1))
                - float(gammaln(base.LOCAL_BITS - weight + 1))
            )

        self.nonzero = self.weights != 0
        self.supports = np.arange(base.PACKETS_PER_LOCAL + 1, dtype=np.float64)
        self.log_support_probabilities = np.full(
            (len(self.weights), base.PACKETS_PER_LOCAL + 1),
            -math.inf,
            dtype=np.float64,
        )
        support_polynomials = base.packet_support_coefficients()
        for weight in range(base.LOCAL_BITS + 1):
            if weight == 0:
                self.log_support_probabilities[weight, 0] = 0.0
                continue
            denominator = math.comb(base.LOCAL_BITS, weight)
            observed = 0
            for support in range(1, base.PACKETS_PER_LOCAL + 1):
                polynomial = support_polynomials[support]
                if weight >= len(polynomial):
                    continue
                coefficient = polynomial[weight]
                if coefficient == 0:
                    continue
                count = math.comb(base.PACKETS_PER_LOCAL, support) * coefficient
                observed += count
                self.log_support_probabilities[weight, support] = (
                    math.log(count) - math.log(denominator)
                )
            if observed != denominator:
                raise AssertionError(
                    f"packet-support law failed at weight {weight}"
                )

        self.spectrum_receipt = {
            "kind": "EXACT_RANDOM_LINEAR_ENSEMBLE_EXPECTATION",
            "length": base.LOCAL_BITS,
            "dimension": base.LOCAL_INPUT_BITS,
            "expected_weight_formula": "((2^256-1)/(2^512-1))*binom(512,w)",
            "expected_packet_support_formula": "((2^256-1)/(2^512-1))*binom(128,s)*15^s",
            "typical_minimum_distance": "about 58; no deterministic local floor is imposed",
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inner-receipt", type=Path, default=DEFAULT_INNER_RECEIPT)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    grouped = load_inner_caps(args.inner_receipt)
    if args.relative_distance not in grouped:
        raise ValueError("inner receipt does not contain the requested distance")
    model = RandomLinearMoments()
    result = base.evaluate(model, grouped[args.relative_distance])
    result["canonical_late_placement_test"]["validity"] = (
        "The event charges the maximum 64 packets from the two encoded field "
        "parities. The data coefficient is exact in expectation over the "
        "random-code setup experiment."
    )
    result["canonical_late_placement_test"]["inner_event"]["scope"] = (
        "This per-word event charges the maximum parity support for every data "
        "configuration at the representative support."
    )

    payload = {
        "schema": "riffle-frozenrandom512-p2-randomstepconv-v2",
        "model": "Riffle FrozenRandom512-P2-RandomStepConv g=4 sigma=20",
        "parameters": {
            "data_local_code": [base.LOCAL_BITS, base.LOCAL_INPUT_BITS],
            "data_groups": base.DATA_GROUPS,
            "independent_setup_maps": True,
            "field_parity_symbols": 2,
            "parity_symbol_code": [128, 64, 22],
            "packet_positions": base.TARGET_PACKETS,
            "packet_bits": base.PACKET_BITS,
            "relative_binary_weight": args.relative_distance,
        },
        "spectrum": model.spectrum_receipt,
        "evidence": {
            "data_local_weight_spectrum": "EXACT_ENSEMBLE_EXPECTATION",
            "packet_support_law_given_weight": "EXACT_INTEGER_VALIDATED",
            "data_configuration_sum": "EXACT_ENSEMBLE_FIRST_MOMENT",
            "parity_block": "RIGOROUS_WORST_CASE_MOMENT",
            "inner_caps": "OUTWARD_ROUNDED_FROM_ONE_LAP_GOAL02",
            "optimization_fft_and_final_sum": "FLOATING_DIAGNOSTIC",
        },
        "result": result,
        "scope": (
            "The result concerns independently sampled random linear local "
            "encoders that are frozen after setup. It isolates the structural "
            "outer spectrum and does not claim an efficient dense encoder."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"total_log2,{result['total_log2_first_moment_upper']:.6f}")
    print(f"dominant_support,{result['dominant_support_interval']}")
    print(f"dominant_mode,{result['dominant_occupation_mode']:.0f}")
    print(
        "pointwise_upper_log2,"
        f"{result['representative_support_audit']['pointwise_first_moment_upper_log2']:.6f}"
    )
    print(
        "explicit_event_expected_log2,"
        f"{result['canonical_late_placement_test']['modeled_expected_count_log2']:.6f}"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
