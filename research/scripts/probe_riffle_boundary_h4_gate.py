#!/usr/bin/env python3
"""Probe a uniform lower-tail gate for weight-four packets.

For a BCH block of binary weight r, the within-block permutation gives an
exact distribution for the number of full four-bit packets.  A one-parameter
Chernoff envelope then bounds the sum across any outer word using only its
total binary weight.  Floating-point optimization makes this a diagnostic;
the underlying one-block counts are exact integers.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from analyze_riffle_parity_ladder_zero import DEFAULT_SPECTRUM, load_spectrum


PACKETS_PER_BLOCK = 32
BLOCK_LENGTH = 128
NONFULL_PACKET_POLYNOMIAL = (1, 4, 6, 4)
DEFAULT_DISTANCES = (76_000, 77_000, 188_743)
DEFAULT_THRESHOLDS = (0, 1, 8, 16, 32, 64, 96, 128)
CERTIFICATE_THRESHOLD = 64
CERTIFICATE_S_NUMERATOR = 27
CERTIFICATE_S_DENOMINATOR = 64
CERTIFICATE_U_BITS = 64


def polynomial_powers() -> list[list[int]]:
    powers = [[1]]
    for _ in range(PACKETS_PER_BLOCK):
        old = powers[-1]
        following = [0] * (len(old) + len(NONFULL_PACKET_POLYNOMIAL) - 1)
        for left, left_value in enumerate(old):
            for right, right_value in enumerate(NONFULL_PACKET_POLYNOMIAL):
                following[left + right] += left_value * right_value
        powers.append(following)
    return powers


def full_packet_counts(weight: int, powers: list[list[int]]) -> list[int]:
    counts = []
    for full_packets in range(PACKETS_PER_BLOCK + 1):
        remainder = weight - 4 * full_packets
        polynomial = powers[PACKETS_PER_BLOCK - full_packets]
        coefficient = polynomial[remainder] if 0 <= remainder < len(polynomial) else 0
        counts.append(math.comb(PACKETS_PER_BLOCK, full_packets) * coefficient)
    if sum(counts) != math.comb(BLOCK_LENGTH, weight):
        raise RuntimeError(f"full-packet distribution failed at weight {weight}")
    return counts


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def log_mgf(log_s: float, counts: list[int]) -> float:
    denominator = math.log(sum(counts))
    terms = [
        math.log(count) + full_packets * log_s
        for full_packets, count in enumerate(counts)
        if count
    ]
    return logsumexp(terms) - denominator


def uniform_envelope(
    log_s: float,
    distributions: dict[int, list[int]],
) -> tuple[float, int]:
    rows = [
        (log_mgf(log_s, counts) / weight, weight)
        for weight, counts in distributions.items()
    ]
    return max(rows)


def threshold_bound(
    input_weight: int,
    threshold: int,
    distributions: dict[int, list[int]],
) -> dict[str, object]:
    def objective(log_s: float) -> float:
        envelope, _weight = uniform_envelope(log_s, distributions)
        return input_weight * envelope - threshold * log_s

    result = minimize_scalar(
        objective,
        method="bounded",
        bounds=(-40.0, -1e-10),
        options={"xatol": 1e-12, "maxiter": 2000},
    )
    envelope, active_weight = uniform_envelope(float(result.x), distributions)
    natural_log_bound = input_weight * envelope - threshold * float(result.x)
    return {
        "threshold_h4_at_most": threshold,
        "optimized_log_s": float(result.x),
        "active_bch_weight_in_envelope": active_weight,
        "natural_log_upper_bound": natural_log_bound,
        "log2_upper_bound": min(0.0, natural_log_bound / math.log(2.0)),
        "optimizer_reported_success": bool(result.success),
    }


def exact_rational_certificate(
    input_weight: int,
    threshold: int,
    distributions: dict[int, list[int]],
) -> dict[str, object]:
    """Certify one Chernoff point with exact integer comparisons."""
    s_numerator = CERTIFICATE_S_NUMERATOR
    s_denominator = CERTIFICATE_S_DENOMINATOR
    u_denominator = 1 << CERTIFICATE_U_BITS
    common_s_denominator = s_denominator**PACKETS_PER_BLOCK
    rational_mgfs = {}
    for weight, counts in distributions.items():
        numerator = sum(
            count
            * s_numerator**full_packets
            * s_denominator ** (PACKETS_PER_BLOCK - full_packets)
            for full_packets, count in enumerate(counts)
        )
        denominator = sum(counts) * common_s_denominator
        rational_mgfs[weight] = (numerator, denominator)

    def envelope_holds(u_numerator: int) -> bool:
        return all(
            numerator * u_denominator**weight
            <= denominator * u_numerator**weight
            for weight, (numerator, denominator) in rational_mgfs.items()
        )

    low = 0
    high = u_denominator
    while low + 1 < high:
        middle = (low + high) // 2
        if envelope_holds(middle):
            high = middle
        else:
            low = middle
    u_numerator = high
    if not envelope_holds(u_numerator) or envelope_holds(u_numerator - 1):
        raise RuntimeError("failed to isolate the least dyadic envelope")

    left = s_denominator**threshold * u_numerator**input_weight
    right = s_numerator**threshold * u_denominator**input_weight
    certified_bits = right.bit_length() - left.bit_length()
    while certified_bits > 0 and (left << certified_bits) > right:
        certified_bits -= 1
    while (left << (certified_bits + 1)) <= right:
        certified_bits += 1
    approximate_log2_bound = (
        threshold * math.log2(s_denominator / s_numerator)
        + input_weight * math.log2(u_numerator / u_denominator)
    )
    return {
        "threshold_h4_at_most": threshold,
        "s_rational": {
            "numerator": s_numerator,
            "denominator": s_denominator,
        },
        "u_dyadic": {
            "numerator": str(u_numerator),
            "denominator_power_of_two": CERTIFICATE_U_BITS,
        },
        "all_bch_weight_mgf_inequalities_exactly_verified": True,
        "exact_final_inequality_verified": (left << certified_bits) <= right,
        "certified_whole_bits": certified_bits,
        "approximate_log2_bound": approximate_log2_bound,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument(
        "--distances", type=int, nargs="+", default=list(DEFAULT_DISTANCES)
    )
    parser.add_argument(
        "--thresholds", type=int, nargs="+", default=list(DEFAULT_THRESHOLDS)
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    active_weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    powers = polynomial_powers()
    distributions = {
        weight: full_packet_counts(weight, powers) for weight in active_weights
    }

    minimum_expectation_per_input_bit = min(
        sum(index * count for index, count in enumerate(counts))
        / sum(counts)
        / weight
        for weight, counts in distributions.items()
    )
    rows = []
    for distance in args.distances:
        input_weight = 2 * distance
        rows.append(
            {
                "distance": distance,
                "boundary_input_weight": input_weight,
                "uniform_minimum_expected_h4": (
                    input_weight * minimum_expectation_per_input_bit
                ),
                "threshold_bounds": [
                    threshold_bound(input_weight, threshold, distributions)
                    for threshold in args.thresholds
                ],
                "exact_rational_h4_at_most_64_certificate": (
                    exact_rational_certificate(
                        input_weight,
                        CERTIFICATE_THRESHOLD,
                        distributions,
                    )
                ),
            }
        )

    payload = {
        "schema": "riffle-boundary-h4-gate-v1",
        "evidence_label": (
            "EXACT_LOCAL_COUNTS_AND_H4_64_RATIONAL_CERTIFICATE_WITH_"
            "ADDITIONAL_NUMERICAL_OPTIMIZATION"
        ),
        "active_bch_weights": active_weights,
        "minimum_expected_h4_per_input_bit": minimum_expectation_per_input_bit,
        "rows": rows,
        "scope": (
            "Conditional on any fixed outer word, the one-block distributions "
            "and the max-over-BCH-weights envelope are valid. The displayed "
            "optimized logs are numerical diagnostics, not outward-rounded "
            "certificates. This gate bounds packetization probability only; "
            "it does not by itself sum the outer words."
        ),
    }
    output = args.output or Path(
        "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/"
        "receipts/goal26_boundary_h4_gate.json"
    )
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(output),
                "minimum_expected_h4_per_input_bit": minimum_expectation_per_input_bit,
                "rows": [
                    {
                        "distance": row["distance"],
                        "minimum_expected_h4": row["uniform_minimum_expected_h4"],
                        "bounds": {
                            str(bound["threshold_h4_at_most"]): bound["log2_upper_bound"]
                            for bound in row["threshold_bounds"]
                        },
                        "certified_h4_at_most_64_whole_bits": row[
                            "exact_rational_h4_at_most_64_certificate"
                        ]["certified_whole_bits"],
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
