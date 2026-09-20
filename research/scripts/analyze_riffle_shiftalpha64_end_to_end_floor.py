#!/usr/bin/env python3
"""Prove a finite end-to-end distance floor for ShiftAlpha64.

The inner sweep gives a uniform fixed-word bound rho^w.  The outer sum treats
one-active-data messages through the exact shifted minimum weight 72.  For
messages with at least two active data blocks, it drops both parity blocks and
uses the exact extended-BCH weight enumerator on the systematic blocks.
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

from analyze_riffle_denseouter_parallelacc_g4_goal02 import (  # noqa: E402
    exact_small_support_profile,
    optimize_profile,
)


DATA_BLOCKS = 1 << 14
FIELD_NONZERO = (1 << 64) - 1
PACKET_POSITIONS = 524_352
PACKET_WIDTH = 4
GLOBAL_BINARY_LENGTH = PACKET_POSITIONS * PACKET_WIDTH
ONE_DATA_MINIMUM_WEIGHT = 72
OUTER_FIELD_BLOCK_DISTANCE = 3
BCH_MINIMUM_WEIGHT = 22
OUTER_BINARY_MINIMUM_WEIGHT = OUTER_FIELD_BLOCK_DISTANCE * BCH_MINIMUM_WEIGHT
DETERMINISTIC_OUTPUT_MINIMUM_WEIGHT = (OUTER_BINARY_MINIMUM_WEIGHT + 1) // 2
DEFAULT_SPECTRUM = SCRIPT_DIRECTORY / "ebch128_64_spectrum.csv"


def load_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        spectrum = {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(source)
        }
    if sum(spectrum.values()) != 1 << 64:
        raise RuntimeError("EBCH spectrum failed total-mass validation")
    if spectrum.get(0) != 1 or spectrum.get(22) != 243_840:
        raise RuntimeError("EBCH spectrum failed endpoint or minimum-shell validation")
    if any(spectrum.get(weight, 0) != spectrum.get(128 - weight, 0) for weight in spectrum):
        raise RuntimeError("EBCH spectrum failed complement symmetry")
    return spectrum


def inner_sweep(distance: int) -> dict[str, object]:
    profiles = []
    worst = None
    for input_weight in range(1, 2 * distance + 1):
        support_minimum = (input_weight + PACKET_WIDTH - 1) // PACKET_WIDTH
        for support in range(support_minimum, input_weight + 1):
            if support <= 2:
                profile = exact_small_support_profile(
                    PACKET_POSITIONS,
                    PACKET_WIDTH,
                    distance,
                    input_weight,
                    support,
                )
            else:
                profile = optimize_profile(
                    PACKET_POSITIONS,
                    PACKET_WIDTH,
                    distance,
                    input_weight,
                    support,
                )
            profiles.append(profile)
            if worst is None or profile["root"] > worst["root"]:
                worst = profile
    if worst is None:
        raise RuntimeError("inner profile sweep is empty")
    return {
        "uniform_rho": worst["root"],
        "uniform_rho_log2": math.log2(worst["root"]),
        "profile_count": len(profiles),
        "failed_optimizers": sum(not profile["optimizer_success"] for profile in profiles),
        "worst_profile": worst,
        "justification": (
            "Profiles with input weight at most 2D are swept; larger input "
            "weight cannot produce output weight at most D. Every returned "
            "numerical point is a valid moment bound regardless of optimizer status."
        ),
    }


def outer_sum(rho: float, spectrum: dict[int, int]) -> dict[str, object]:
    block_partition = sum(
        count * rho**weight for weight, count in spectrum.items()
    )
    active_mass = block_partition - 1.0
    all_systematic_nonzero = math.expm1(DATA_BLOCKS * math.log1p(active_mass))
    one_systematic_term = DATA_BLOCKS * active_mass
    multi_data_bound = all_systematic_nonzero - one_systematic_term
    one_data_bound = (
        DATA_BLOCKS * FIELD_NONZERO * rho**ONE_DATA_MINIMUM_WEIGHT
    )
    total = multi_data_bound + one_data_bound
    return {
        "exact_bch_block_partition": block_partition,
        "nonzero_bch_block_mass": active_mass,
        "one_data_bound": one_data_bound,
        "multi_data_bound_after_dropping_parities": multi_data_bound,
        "total_expected_bad_codewords_bound": total,
        "failure_bits": -math.log2(total),
        "success_probability_lower_bound": max(0.0, 1.0 - total),
        "identity": (
            "(A(rho)^B - 1 - B(A(rho)-1)) + B(2^64-1)rho^72"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    args = parser.parse_args()
    if args.distance < 1:
        raise ValueError("distance threshold must be positive")

    spectrum = load_spectrum(args.spectrum)
    inner = inner_sweep(args.distance)
    outer = outer_sum(float(inner["uniform_rho"]), spectrum)
    deterministic_impossibility = args.distance < DETERMINISTIC_OUTPUT_MINIMUM_WEIGHT
    first_moment_passes = (
        outer["total_expected_bad_codewords_bound"] < 1.0
    )
    if deterministic_impossibility:
        failure_probability = 0.0
        failure_bits: float | str = "infinity"
        proof_source = "deterministic outer-distance/accumulator inequality"
        proved_minimum_distance = DETERMINISTIC_OUTPUT_MINIMUM_WEIGHT
    elif first_moment_passes:
        failure_probability = outer["total_expected_bad_codewords_bound"]
        failure_bits = outer["failure_bits"]
        proof_source = "first moment over all nonzero messages"
        proved_minimum_distance = args.distance + 1
    else:
        failure_probability = 1.0
        failure_bits = 0.0
        proof_source = "deterministic baseline; first moment is trivial"
        proved_minimum_distance = DETERMINISTIC_OUTPUT_MINIMUM_WEIGHT
    payload = {
        "schema": "riffle-shiftalpha64-end-to-end-floor-v1",
        "evidence_label": "RIGOROUS_FINITE_DISTANCE_BOUND",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "parameters": {
            "data_blocks": DATA_BLOCKS,
            "global_binary_length": GLOBAL_BINARY_LENGTH,
            "bad_output_weight_inclusive": args.distance,
            "proved_minimum_distance_lower_bound": proved_minimum_distance,
            "one_data_outer_minimum_weight": ONE_DATA_MINIMUM_WEIGHT,
            "outer_field_block_distance": OUTER_FIELD_BLOCK_DISTANCE,
            "bch_minimum_weight": BCH_MINIMUM_WEIGHT,
            "outer_binary_minimum_weight": OUTER_BINARY_MINIMUM_WEIGHT,
            "deterministic_output_minimum_weight": (
                DETERMINISTIC_OUTPUT_MINIMUM_WEIGHT
            ),
        },
        "inner_fixed_word_contraction": inner,
        "outer_first_moment": outer,
        "theorem": {
            "failure_event": f"minimum distance is at most {args.distance}",
            "failure_probability_upper_bound": failure_probability,
            "failure_bits": failure_bits,
            "proof_source": proof_source,
            "deterministic_event_impossible": deterministic_impossibility,
            "passes_positive_margin": failure_probability < 1.0,
            "passes_0_1_bit_margin": (
                deterministic_impossibility
                or float(failure_bits) >= 0.1
            ),
            "passes_20_bit_margin": (
                deterministic_impossibility
                or float(failure_bits) >= 20.0
            ),
        },
        "validation": {
            "exact_bch_spectrum": "PASS",
            "systematic_support_partition": "EXACT",
            "one_data_replacement_uses_exact_minimum_72": True,
            "multi_data_parity_blocks_are_dropped_for_an_upper_bound": True,
            "outer_distance_times_accumulator_half_distance": (
                "3 * 22 / 2 gives deterministic output distance 33"
            ),
        },
        "scope": (
            "Finite ensemble theorem at the stated length. It proves a distance "
            "floor, not positive relative distance as length grows."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
