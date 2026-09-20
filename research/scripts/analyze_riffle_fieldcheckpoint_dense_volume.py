#!/usr/bin/env python3
"""Dense-subspace Hamming-ball bound for FieldCheckpointAccumulate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_fieldcheckpoint_regular_envelope import (
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import (
    LOG2,
    log_choose,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal14_dense_volume_delta06.json"
)


def hamming_ball_log(length: int, radius: int) -> float:
    weights = np.arange(radius + 1, dtype=np.float64)
    logs = (
        gammaln(length + 1.0)
        - gammaln(weights + 1.0)
        - gammaln(length - weights + 1.0)
    )
    return float(logsumexp(logs))


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta
    regular_point_log_density = regular_log_mass - args.outer_bits * LOG2
    ball_log = hamming_ball_log(output_bits, distance)

    rows = []
    regular_terms = []
    mixed_terms = []
    for regular in range(1, outer_blocks + 1):
        placement_log = log_choose(outer_blocks, regular)
        regular_log = (
            placement_log + ball_log + regular * regular_point_log_density
        )
        # Sum the all-one count b exactly at the combinatorial level:
        # sum_b C(N,a) C(N-a,b) = C(N,a) 2^(N-a).
        mixed_log = regular_log + (outer_blocks - regular) * LOG2
        regular_terms.append(regular_log)
        mixed_terms.append(mixed_log)
        rows.append(
            {
                "regular_active_blocks": regular,
                "regular_pointwise_log2_upper": regular_log / LOG2,
                "any_all_one_count_pointwise_log2_upper": mixed_log / LOG2,
            }
        )

    def first_tail_below(terms: list[float], target_bits: float) -> tuple[int, float]:
        tail = -math.inf
        for index in range(len(terms) - 1, -1, -1):
            tail = float(np.logaddexp(tail, terms[index]))
            if tail >= -target_bits * LOG2:
                return index + 2, float(
                    logsumexp(np.asarray(terms[index + 1 :]))
                )
        return 1, float(logsumexp(np.asarray(terms)))

    regular_start, regular_tail = first_tail_below(
        regular_terms, args.target_margin_bits
    )
    mixed_start, mixed_tail = first_tail_below(
        mixed_terms, args.target_margin_bits
    )
    return {
        "schema": "riffle-fieldcheckpoint-dense-volume-v1",
        "candidate": "Dense-volume bound for any invertible inner",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "target_margin_bits": args.target_margin_bits,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "regular_point_log2_density": regular_point_log_density / LOG2,
        },
        "hamming_ball_log2_upper_float": ball_log / LOG2,
        "regular_tail": {
            "first_active_count": regular_start,
            "log2_upper_float": regular_tail / LOG2,
            "lambda_bits_lower_float": -regular_tail / LOG2,
        },
        "any_all_one_count_tail": {
            "first_regular_active_count": mixed_start,
            "log2_upper_float": mixed_tail / LOG2,
            "lambda_bits_lower_float": -mixed_tail / LOG2,
        },
        "rows": rows,
        "method": (
            "The inner map is invertible.  For fixed active-block positions, "
            "the regular-support envelope is uniform on a 256a-dimensional "
            "coordinate subspace up to its density factor.  Its image meets "
            "a Hamming ball in at most the volume of that ball.  The mixed "
            "bound sums every possible all-one-block placement."
        ),
        "scope": (
            "Nearest-binary64 diagnostic.  The formulas are analytic upper "
            "bounds, but the reported numbers are not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--relative-distance", type=float, default=0.06)
    parser.add_argument("--target-margin-bits", type=float, default=40.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "regular_tail_start,"
        f"{payload['regular_tail']['first_active_count']},"
        "mixed_tail_start,"
        f"{payload['any_all_one_count_tail']['first_regular_active_count']}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
