#!/usr/bin/env python3
"""Evaluate the minimum-weight gate for a Hölder product-envelope proof.

The exact EBCH [128,64,22] spectrum supplies A_22.  If a product envelope
must cover Q(22,22,22), Hölder's inequality charges at least one A_22-sized
marginal shell per data position.  This script compares that charge with the
current three-block bound and with the scalar full-bit-permutation benchmark.
It prints JSON and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

from analyze_riffle_bchblockperm_parallelacc_g4_goal01 import (  # noqa: E402
    scalar_accumulator_bad_probability,
)


DEFAULT_SPECTRUM = SCRIPT_DIRECTORY / "ebch128_64_spectrum.csv"
DATA_POSITIONS = 1 << 14
GLOBAL_BINARY_LENGTH = 2_097_408
DISTANCE = 188_766
MINIMUM_TRIPLE_WEIGHT = 66


def load_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        rows = csv.DictReader(source)
        spectrum = {int(row["weight"]): int(row["count"]) for row in rows}
    if sum(spectrum.values()) != 1 << 64:
        raise RuntimeError("EBCH spectrum failed total-mass validation")
    if spectrum.get(0) != 1 or spectrum.get(128) != 1:
        raise RuntimeError("EBCH spectrum failed endpoint validation")
    if any(spectrum.get(weight, 0) != spectrum.get(128 - weight, 0) for weight in spectrum):
        raise RuntimeError("EBCH spectrum failed complement symmetry")
    if min(weight for weight, count in spectrum.items() if weight and count) != 22:
        raise RuntimeError("EBCH spectrum failed minimum-distance validation")
    return spectrum


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument(
        "--current-profile-log2-bound", type=float, default=-36.4472162884044
    )
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    minimum_count = spectrum[22]
    shell_log2 = math.log2(DATA_POSITIONS) + math.log2(minimum_count)
    current_total_log2 = shell_log2 + args.current_profile_log2_bound

    scalar_probability = scalar_accumulator_bad_probability(
        GLOBAL_BINARY_LENGTH, MINIMUM_TRIPLE_WEIGHT, DISTANCE
    )
    scalar_log2 = math.log2(scalar_probability)
    scalar_total_log2 = shell_log2 + scalar_log2

    targets = {}
    for target_bits in (20, 40):
        required_profile_log2 = -target_bits - shell_log2
        targets[str(target_bits)] = {
            "required_profile_log2_bound": required_profile_log2,
            "current_shortfall_bits": (
                args.current_profile_log2_bound - required_profile_log2
            ),
            "full_bit_permutation_shortfall_bits": (
                scalar_log2 - required_profile_log2
            ),
            "current_passes": current_total_log2 <= -target_bits,
            "full_bit_permutation_passes": scalar_total_log2 <= -target_bits,
        }

    payload = {
        "schema": "riffle-bch-holder-minimum-shell-gate-v1",
        "evidence_label": "EXACT_SPECTRUM_ARITHMETIC_WITH_SCALAR_IDEAL_BENCHMARK",
        "parameters": {
            "data_positions": DATA_POSITIONS,
            "minimum_bch_weight": 22,
            "minimum_three_block_weight": MINIMUM_TRIPLE_WEIGHT,
            "distance_threshold": DISTANCE,
            "global_binary_length": GLOBAL_BINARY_LENGTH,
        },
        "exact_bch_spectrum": {
            "A_22": minimum_count,
            "log2_A_22": math.log2(minimum_count),
            "validation": "PASS",
            "source": str(args.spectrum),
        },
        "holder_minimum_shell": {
            "position_and_message_log2_mass": shell_log2,
            "inequality": "sum_x a_x^2 a_(alpha*x) <= sum_x a_x^3",
            "current_profile_log2_bound": args.current_profile_log2_bound,
            "current_shell_log2_bound": current_total_log2,
            "interpretation": (
                "An envelope justified only by the current certified "
                "Q(22,22,22) upper bound inherits this minimum-shell charge."
            ),
        },
        "full_global_bit_permutation_comparison": {
            "fixed_weight_66_log2_probability": scalar_log2,
            "minimum_shell_log2_bound": scalar_total_log2,
            "role": "optimistic scalar-accumulator comparison, not this construction",
        },
        "targets": targets,
        "scope": (
            "This gate concerns the one-data-symbol minimum-weight shell under "
            "a correlation-free product envelope. It neither bounds other "
            "profiles nor uses the shifted alpha spectrum."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
