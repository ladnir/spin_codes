#!/usr/bin/env python3
"""Diagnostic spectrum budget for the unpadded EBCH128 construction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import evaluate_ebch128_randomstepconv_g1 as base


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128_pow2_spectrum_budget.json"


def evaluate() -> dict[str, object]:
    base.configure_outer_rows(base.ACTIVE_ROWS)
    spectrum = base.load_spectrum(base.SPECTRUM)
    body_mass = (1 << base.K) - 2
    even_body_mass = (1 << (base.B - 1)) - 2

    kl_bits = 0.0
    second_moment = 0.0
    maximum_log2_ratio = -math.inf
    maximum_weights: list[int] = []
    for weight, count in spectrum.items():
        if weight in (0, base.B) or count == 0:
            continue
        probability = count / body_mass
        reference = math.comb(base.B, weight) / even_body_mass
        log2_ratio = math.log2(probability / reference)
        kl_bits += probability * log2_ratio
        second_moment += probability * (probability / reference)
        if log2_ratio > maximum_log2_ratio + 1e-14:
            maximum_log2_ratio = log2_ratio
            maximum_weights = [weight]
        elif abs(log2_ratio - maximum_log2_ratio) <= 1e-14:
            maximum_weights.append(weight)

    renyi2_bits = math.log2(second_moment)
    log_term = 0.0
    log_hamming_ball = 0.0
    for weight in range(base.D - 1):
        log_term += math.log(base.N - weight) - math.log(weight + 1)
        log_hamming_ball = (
            max(log_hamming_ball, log_term)
            + math.log1p(math.exp(-abs(log_hamming_ball - log_term)))
        )
    message_bits = base.K * base.ACTIVE_ROWS
    random_margin = (
        base.N
        - message_bits
        - log_hamming_ball / math.log(2.0)
    )
    return {
        "schema": "ebch128-pow2-spectrum-budget-v1",
        "status": "DIAGNOSTIC",
        "parameters": {
            "message_bits": base.K * base.ACTIVE_ROWS,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "outer_rows": base.L,
        },
        "uniform_random_linear_code": {
            "first_moment_margin_bits": random_margin,
            "event": "nonzero codeword weight strictly below the distance cutoff",
        },
        "body_spectrum_vs_uniform_even_body": {
            "excluded_weights": [0, base.B],
            "kl_bits_per_row": kl_bits,
            "kl_bits_at_full_occupation": kl_bits * base.L,
            "renyi2_bits_per_row": renyi2_bits,
            "renyi2_bits_at_full_occupation": renyi2_bits * base.L,
            "max_log2_density_ratio_per_row": maximum_log2_ratio,
            "max_log2_density_ratio_at_full_occupation": maximum_log2_ratio * base.L,
            "maximizing_weights": maximum_weights,
        },
        "interpretation": [
            "The existing dense proof pays the maximum density ratio independently in every active row.",
            "The small KL and Renyi-2 values show that this pointwise payment is not representative of the average body spectrum.",
            "These diagnostics do not themselves imply a distance bound; a shell-sensitive transfer is still required.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = evaluate()
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
