#!/usr/bin/env python3
"""Test whether a two-stage sparse-EA mean spectrum can reuse frozen caps.

This is a nearest-binary64 design diagnostic.  It computes the complete
one-word spectrum of ``A E1 A E0`` and compares it with the exact integer caps
from the closed degree-33 theorem.  It also reports the largest hypothetical
uniform variance factor compatible with selected setup and end-to-end
failure margins.  It does not prove a two-stage variance bound.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from sweep_one_stage_extreme_shells import activation_probabilities
from sweep_two_stage_extreme_shells import (
    first_stage_expected_spectrum,
    full_accumulator_shells,
)


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
DEFAULT_CAPS = (
    WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_spectrum_caps_outward.json"
)
DEFAULT_TRANSFER = (
    WORKSTREAM
    / "one_stage_sparse_ea_K256_B512_r33_randomstepconv_M22_d109_outward.json"
)
DEFAULT_OUTPUT = HERE / "two_stage_cap_reuse_B512_diagnostic.json"


def parse_candidate(value: str) -> tuple[int, int]:
    try:
        first, second = (int(part) for part in value.split(","))
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("candidate must have form R0,R1") from error
    if first < 1 or second < 1:
        raise argparse.ArgumentTypeError("degrees must be positive")
    return first, second


def bits(probability: float) -> float:
    if probability <= 0:
        return math.inf
    return -math.log2(probability)


def setup_failure(
    means: np.ndarray,
    caps: list[int],
    variance_factor: float,
    rank_attempts: int,
    lower_support: int,
    upper_support: int,
) -> tuple[float, float, float, float]:
    kernel = float(means[0])
    tail = float(
        np.sum(means[1:lower_support])
        + np.sum(means[upper_support + 1 :])
    )
    cap_failure = 0.0
    for shell in range(lower_support, upper_support + 1):
        mean = float(means[shell])
        deviation = float(caps[shell] + 1) - mean
        if deviation <= 0:
            return math.inf, kernel, tail, math.inf
        variance = variance_factor * mean
        cap_failure += variance / (variance + deviation * deviation)
    if kernel >= 1:
        return math.inf, kernel, tail, cap_failure
    conditioned = (tail + cap_failure) / (1.0 - kernel)
    abort = kernel**rank_attempts
    return conditioned + abort, kernel, tail, cap_failure


def largest_factor(
    means: np.ndarray,
    caps: list[int],
    rank_attempts: int,
    lower_support: int,
    upper_support: int,
    target_margin: float,
    conditional_bad: float,
) -> float:
    def passes(factor: float) -> bool:
        setup, *_ = setup_failure(
            means,
            caps,
            factor,
            rank_attempts,
            lower_support,
            upper_support,
        )
        return bits(setup + conditional_bad) >= target_margin

    low = 0.0
    high = 1.0
    while passes(high) and high < 2.0**40:
        low = high
        high *= 2.0
    for _ in range(100):
        middle = (low + high) / 2.0
        if passes(middle):
            low = middle
        else:
            high = middle
    return low


def evaluate(
    first_degree: int,
    second_degree: int,
    second_transition: np.ndarray,
    caps: list[int],
    rank_attempts: int,
    lower_support: int,
    upper_support: int,
    conditional_bad: float,
) -> dict[str, object]:
    first_spectrum = first_stage_expected_spectrum(
        message_bits=256,
        output_bits=512,
        right_degree=first_degree,
    )
    means = first_spectrum @ second_transition
    expected_mass = float(np.sum(means))
    expected_mass_reference = float((1 << 256) - 1)

    cap_headrooms = []
    random_excesses = []
    for shell in range(lower_support, upper_support + 1):
        mean = float(means[shell])
        cap = caps[shell]
        cap_headrooms.append((math.log2(cap) - math.log2(mean), shell))
        random_mean = math.ldexp(float(math.comb(512, shell)), -256)
        random_excesses.append((math.log2(mean / random_mean), shell))

    setup_rows = []
    for factor in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512):
        setup, _kernel, _tail, cap_failure = setup_failure(
            means,
            caps,
            float(factor),
            rank_attempts,
            lower_support,
            upper_support,
        )
        setup_rows.append(
            {
                "uniform_variance_to_mean_factor": factor,
                "positive_cap_failure_margin_bits": bits(cap_failure),
                "setup_margin_bits": bits(setup),
                "end_to_end_margin_bits_if_frozen_transfer_reused": bits(
                    setup + conditional_bad
                ),
            }
        )

    setup_at_one, kernel, tail, _ = setup_failure(
        means,
        caps,
        1.0,
        rank_attempts,
        lower_support,
        upper_support,
    )
    minimum_headroom, minimum_headroom_shell = min(cap_headrooms)
    maximum_random_excess, maximum_random_excess_shell = max(random_excesses)
    return {
        "first_right_degree": first_degree,
        "second_right_degree": second_degree,
        "xor_count": 512 * (first_degree + second_degree) - 2,
        "expected_spectrum_mass": expected_mass,
        "expected_spectrum_mass_relative_error": (
            expected_mass / expected_mass_reference - 1.0
        ),
        "kernel_expected_count": kernel,
        "kernel_expected_count_margin_bits": bits(kernel),
        "zero_cap_tail_expected_count": tail,
        "zero_cap_tail_margin_bits": bits(tail),
        "all_expected_positive_shells_below_frozen_caps": all(
            means[shell] < caps[shell]
            for shell in range(lower_support, upper_support + 1)
        ),
        "minimum_cap_to_mean_headroom_bits": minimum_headroom,
        "minimum_cap_to_mean_headroom_shell": minimum_headroom_shell,
        "maximum_excess_over_random_mean_bits": maximum_random_excess,
        "maximum_excess_over_random_mean_shell": maximum_random_excess_shell,
        "setup_margin_bits_if_variance_equals_mean": bits(setup_at_one),
        "largest_uniform_variance_factor_for_41_36_bit_end_to_end_margin": (
            largest_factor(
                means,
                caps,
                rank_attempts,
                lower_support,
                upper_support,
                41.36,
                conditional_bad,
            )
        ),
        "largest_uniform_variance_factor_for_40_bit_end_to_end_margin": (
            largest_factor(
                means,
                caps,
                rank_attempts,
                lower_support,
                upper_support,
                40.0,
                conditional_bad,
            )
        ),
        "uniform_variance_factor_diagnostics": setup_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",
        action="append",
        type=parse_candidate,
        default=None,
        help="degree pair R0,R1; repeat for multiple candidates",
    )
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--caps", type=Path, default=DEFAULT_CAPS)
    parser.add_argument("--transfer", type=Path, default=DEFAULT_TRANSFER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    candidates = args.candidate or [(17, 3), (19, 3)]
    if args.rank_attempts < 1:
        parser.error("rank attempts must be positive")

    cap_receipt = json.loads(args.caps.read_text(encoding="utf-8"))
    transfer_receipt = json.loads(args.transfer.read_text(encoding="utf-8"))
    caps = [int(value) for value in cap_receipt["caps"]]
    lower_support, upper_support = cap_receipt["claim"]["positive_cap_support"]
    conditional_bad = 2.0 ** float(
        transfer_receipt["claim"]["conditional_bad_probability_log2_upper"]
    )

    transitions = {}
    for _first_degree, second_degree in candidates:
        if second_degree not in transitions:
            activations = np.concatenate(
                (
                    np.array([0.0]),
                    activation_probabilities(512, second_degree),
                )
            )
            transitions[second_degree] = full_accumulator_shells(activations, 512)

    rows = [
        evaluate(
            first_degree,
            second_degree,
            transitions[second_degree],
            caps,
            args.rank_attempts,
            lower_support,
            upper_support,
            conditional_bad,
        )
        for first_degree, second_degree in candidates
    ]
    payload = {
        "schema": "two-stage-sparse-ea-cap-reuse-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": "C2 = A E1 A E0",
        "probability_space": {
            "E0": "512 independent uniform fixed-weight rows in F_2^256",
            "E1": "512 independent uniform fixed-weight rows in F_2^512, independent of E0",
            "reuse": "one sampled pair (E0,E1) is reused at every SPIN outer position",
        },
        "frozen_cap_receipt": str(args.caps),
        "frozen_transfer_receipt": str(args.transfer),
        "candidates": rows,
        "scope": [
            "The complete one-word spectra are exact in formula and evaluated in nearest binary64.",
            "Mean domination is necessary but not sufficient for the simultaneous cap event.",
            "Uniform variance factors are hypothetical sensitivity tests, not proved bounds.",
            "Reusing the exact frozen caps would make the existing conditional SPIN transfer applicable without change.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "candidates": [
                    {
                        key: row[key]
                        for key in (
                            "first_right_degree",
                            "second_right_degree",
                            "xor_count",
                            "zero_cap_tail_margin_bits",
                            "all_expected_positive_shells_below_frozen_caps",
                            "minimum_cap_to_mean_headroom_bits",
                            "largest_uniform_variance_factor_for_40_bit_end_to_end_margin",
                        )
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
