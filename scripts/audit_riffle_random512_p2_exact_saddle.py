#!/usr/bin/env python3
"""Exact tilted outer coefficients near a Random512-P2 support saddle."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_rm512_randomstepconv_g4_sigma20 import LOG2
from audit_riffle_random512_p2_global_pointwise import (
    DEFAULT_SPECTRUM,
    GlobalOuterMoment,
    termination_upper_profile,
)


DEFAULT_POINTWISE = Path(
    "constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/"
    "receipts/global_pointwise_p2_sigma18.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/"
    "receipts/global_exact_saddle_p2_sigma18.json"
)


def tilted_global_distribution(
    outer: GlobalOuterMoment, log_t: float
) -> tuple[np.ndarray, float]:
    block_length = outer.local_output_packets + 1
    supports = np.arange(block_length, dtype=np.float64)
    log_block = np.full(block_length, -math.inf)
    log_block[0] = 0.0
    for support in range(1, block_length):
        log_block[support] = (
            outer.log_output_nonzero_ratio
            + float(gammaln(block_length))
            - float(gammaln(support + 1))
            - float(gammaln(block_length - support))
            + support * math.log(15.0)
        )
    log_block_tilted = log_block + supports * log_t
    log_block_normalizer = float(logsumexp(log_block_tilted))
    block = np.exp(log_block_tilted - log_block_normalizer)

    parity_supports = outer.uniform_pair_supports.astype(np.float64)
    log_parity_tilted = outer.log_uniform_pair + parity_supports * log_t
    log_parity_normalizer = float(logsumexp(log_parity_tilted))
    parity = np.exp(log_parity_tilted - log_parity_normalizer)

    maximum_degree = (
        outer.data_groups * outer.local_output_packets + len(parity) - 1
    )
    transform_size = 1 << maximum_degree.bit_length()
    block_transform = np.fft.rfft(
        np.pad(block, (0, transform_size - len(block)))
    )
    parity_transform = np.fft.rfft(
        np.pad(parity, (0, transform_size - len(parity)))
    )
    distribution = np.fft.irfft(
        block_transform**outer.data_groups * parity_transform,
        transform_size,
    )[: maximum_degree + 1]
    minimum = float(distribution.min())
    if minimum < -1e-11:
        raise ArithmeticError("FFT distribution has material negative mass")
    distribution = np.maximum(distribution, 0.0)
    distribution /= distribution.sum()
    log_moment = (
        outer.data_groups * log_block_normalizer + log_parity_normalizer
    )
    return distribution, log_moment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pointwise-receipt", type=Path, default=DEFAULT_POINTWISE)
    parser.add_argument("--ebch-spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--sigma", type=int, default=18)
    parser.add_argument("--center-support", type=int, default=800)
    parser.add_argument("--left", type=int, default=400)
    parser.add_argument("--right", type=int, default=1400)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    pointwise = json.loads(args.pointwise_receipt.read_text(encoding="utf-8"))
    center = min(
        pointwise["rows"],
        key=lambda row: abs(int(row["packet_support"]) - args.center_support),
    )
    log_t = float(center["outer_log_tilt"])
    parity_count = int(pointwise.get("field_parity_symbols", 2))
    outer_bits = int(
        pointwise.get("outer_constituent", {}).get("length_bits", 512)
    )
    outer = GlobalOuterMoment(args.ebch_spectrum, parity_count, outer_bits)
    distribution, log_moment = tilted_global_distribution(outer, log_t)

    rows: list[dict[str, object]] = []
    for support in range(args.left, args.right + 1):
        probability = float(distribution[support])
        if probability <= 0.0:
            continue
        outer_log2 = (
            log_moment + math.log(probability) - support * log_t
        ) / LOG2
        packet_positions = int(pointwise["packet_positions"])
        distance = int(pointwise["distance"])
        inner = termination_upper_profile(
            support, args.sigma, packet_positions, distance
        )
        inner_log2 = float(inner["termination_separated_log2_upper"])
        rows.append(
            {
                "packet_support": support,
                "tilted_probability": probability,
                "outer_exact_coefficient_log2": outer_log2,
                "inner_probability_upper_log2": inner_log2,
                "expected_count_log2": outer_log2 + inner_log2,
                "dominant_termination_count": int(
                    inner["dominant_termination_count"]
                ),
            }
        )
        if support % 100 == 0:
            print(
                f"progress,support,{support},combined,{outer_log2 + inner_log2:.6f}",
                flush=True,
            )

    dominant_rows = sorted(
        rows, key=lambda row: float(row["expected_count_log2"]), reverse=True
    )
    interval_total = float(
        logsumexp([float(row["expected_count_log2"]) * LOG2 for row in rows])
    ) / LOG2
    payload = {
        "schema": "riffle-random512-p2-exact-saddle-v1",
        "state_bits": args.sigma,
        "outer_constituent_bits": outer_bits,
        "support_interval": [args.left, args.right],
        "center_support": int(center["packet_support"]),
        "outer_log_tilt": log_t,
        "tilted_distribution_mass_in_interval": float(
            distribution[args.left : args.right + 1].sum()
        ),
        "interval_log2_first_moment_upper": interval_total,
        "dominant_profiles": dominant_rows[:30],
        "scope": (
            "The FFT computes coefficients of the uniform-parity main term. "
            "The exact global correction is below 2^-512 for occupations at "
            "least two and is omitted here. The interval does not cover all "
            "packet supports."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    top = dominant_rows[0]
    print(
        "dominant_exact,"
        f"support={top['packet_support']},"
        f"combined={top['expected_count_log2']:.6f}"
    )
    print(f"interval_total_log2,{interval_total:.6f}")
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
