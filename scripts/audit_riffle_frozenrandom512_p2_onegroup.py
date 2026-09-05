#!/usr/bin/env python3
"""Audit one active data group with the exact two-parity BCH spectrum."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_frozenrandom512_p2_randomstepconv_g4_sigma20 import (
    RandomLinearMoments,
)
from analyze_riffle_rm512_randomstepconv_g4_sigma20 import DATA_GROUPS, LOG2
from audit_riffle_randomstepconv_episode_structure import (
    termination_upper_profile,
)


EBCH_BITS = 128
EBCH_INPUT_BITS = 64
PACKET_BITS = 4
EBCH_PACKETS = EBCH_BITS // PACKET_BITS
DEFAULT_SPECTRUM = Path(__file__).with_name("EBCH128_64.wd")
DEFAULT_OUTPUT = Path(
    "constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/"
    "receipts/onegroup_exact_p2.json"
)


def read_ebch_spectrum(path: Path) -> dict[int, int]:
    spectrum: dict[int, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        weight_text, count_text = stripped.split()
        spectrum[int(weight_text)] = int(count_text)
    if sum(spectrum.values()) != 1 << EBCH_INPUT_BITS:
        raise ValueError("EBCH spectrum mass is not 2^64")
    if spectrum.get(0) != 1:
        raise ValueError("EBCH spectrum has no unique zero word")
    return spectrum


def packet_support_polynomials(packet_count: int) -> list[list[int]]:
    polynomials: list[list[int]] = [[1]]
    for _support in range(1, packet_count + 1):
        previous = polynomials[-1]
        current = [0] * (len(previous) + PACKET_BITS)
        for degree, coefficient in enumerate(previous):
            current[degree + 1] += 4 * coefficient
            current[degree + 2] += 6 * coefficient
            current[degree + 3] += 4 * coefficient
            current[degree + 4] += coefficient
        polynomials.append(current)
    return polynomials


def ebch_packet_support_distribution(path: Path) -> np.ndarray:
    spectrum = read_ebch_spectrum(path)
    polynomials = packet_support_polynomials(EBCH_PACKETS)
    counts = np.zeros(EBCH_PACKETS + 1, dtype=np.float64)
    for weight, multiplicity in spectrum.items():
        denominator = math.comb(EBCH_BITS, weight)
        observed = 0
        for support in range(EBCH_PACKETS + 1):
            polynomial = polynomials[support]
            if weight >= len(polynomial):
                continue
            coefficient = polynomial[weight]
            if coefficient == 0:
                continue
            placements = math.comb(EBCH_PACKETS, support) * coefficient
            observed += placements
            counts[support] += multiplicity * placements / denominator
        if observed != denominator:
            raise AssertionError(f"EBCH packet law failed at weight {weight}")
    probabilities = counts / float(1 << EBCH_INPUT_BITS)
    probabilities /= probabilities.sum()
    return probabilities


def parity_pair_distribution(single: np.ndarray) -> np.ndarray:
    pair = np.convolve(single, single)
    total_messages = float(1 << 256)
    pair[0] = (total_messages * pair[0] - 1.0) / (total_messages - 1.0)
    pair[1:] *= total_messages / (total_messages - 1.0)
    pair /= pair.sum()
    return pair


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ebch-spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = RandomLinearMoments()
    log_data_by_support = logsumexp(
        data.log_multiplicities[data.nonzero, None]
        + data.log_support_probabilities[data.nonzero],
        axis=0,
    )
    single_parity = ebch_packet_support_distribution(args.ebch_spectrum)
    pair_parity = parity_pair_distribution(single_parity)
    log_pair_parity = np.full_like(pair_parity, -math.inf)
    positive = pair_parity > 0.0
    log_pair_parity[positive] = np.log(pair_parity[positive])

    minimum_total = int(np.flatnonzero(np.isfinite(log_data_by_support))[0])
    maximum_total = len(log_data_by_support) - 1 + len(pair_parity) - 1
    inner = {
        support: termination_upper_profile(support)
        for support in range(minimum_total, maximum_total + 1)
    }

    rows: list[dict[str, object]] = []
    group_choice_log2 = math.log2(DATA_GROUPS)
    for data_support in np.flatnonzero(np.isfinite(log_data_by_support)):
        if data_support == 0:
            continue
        for parity_support in np.flatnonzero(positive):
            total_support = int(data_support + parity_support)
            outer_log2 = (
                float(log_data_by_support[data_support])
                + float(log_pair_parity[parity_support])
            ) / LOG2 + group_choice_log2
            inner_log2 = float(
                inner[total_support]["termination_separated_log2_upper"]
            )
            rows.append(
                {
                    "data_packet_support": int(data_support),
                    "parity_packet_support": int(parity_support),
                    "total_packet_support": total_support,
                    "outer_coefficient_log2": outer_log2,
                    "inner_probability_upper_log2": inner_log2,
                    "expected_count_log2": outer_log2 + inner_log2,
                }
            )

    total_log2 = float(
        logsumexp([float(row["expected_count_log2"]) * LOG2 for row in rows])
    ) / LOG2
    dominant = sorted(
        rows, key=lambda row: float(row["expected_count_log2"]), reverse=True
    )
    payload = {
        "schema": "riffle-frozenrandom512-p2-onegroup-v1",
        "model": "Riffle FrozenRandom512-P2-RandomStepConv g=4 sigma=20",
        "outer_occupation": 1,
        "field_parity_symbols": 2,
        "data_groups": DATA_GROUPS,
        "ebch_spectrum": str(args.ebch_spectrum),
        "total_log2_first_moment_upper": total_log2,
        "dominant_profiles": dominant[:20],
        "zero_parity_profile_max": max(
            (
                row
                for row in rows
                if int(row["parity_packet_support"]) == 0
            ),
            key=lambda row: float(row["expected_count_log2"]),
        ),
        "evidence": {
            "data": "exact random-linear ensemble expectation",
            "parity": "exact EBCH [128,64,22] weight spectrum and exact packetization average",
            "inner": "pointwise termination-separated floating upper",
        },
        "scope": (
            "This calculation covers exactly one active 256-bit data group. "
            "The two parity equations are assumed independent on that group, "
            "so their output pair is uniform as the group input ranges over "
            "all messages."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"onegroup_total_log2,{total_log2:.6f}")
    print(
        "dominant_profile,"
        f"data={dominant[0]['data_packet_support']},"
        f"parity={dominant[0]['parity_packet_support']},"
        f"total={dominant[0]['total_packet_support']},"
        f"log2={dominant[0]['expected_count_log2']:.6f}"
    )
    print(
        "zero_parity_max_log2,"
        f"{payload['zero_parity_profile_max']['expected_count_log2']:.6f}"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
