#!/usr/bin/env python3
"""Diagnose a shell-sensitive event for one random [256,128] constituent."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
UNIFORM_RECEIPT = (
    WORKSTREAM / "single_random_constituent_B256_allq_factor8_s22_d109.json"
)
OUTPUT = WORKSTREAM / "single_random_constituent_B256_three_band_event.json"
B = 256
K = 128
FRINGE_CAP = 128
CENTRAL_FACTOR = 1.5


def main() -> None:
    message_count = (1 << K) - 1
    means = np.asarray(
        [
            message_count * math.comb(B, weight) / (1 << B)
            for weight in range(B + 1)
        ]
    )

    extreme_failure = float(
        np.sum(means[1:25]) + np.sum(means[232:257])
    )
    low_mean = float(np.sum(means[25:32]))
    high_mean = float(np.sum(means[225:232]))

    def fringe_failure(mean: float) -> float:
        deviation = (FRINGE_CAP + 1) - mean
        return mean / (mean + deviation * deviation)

    low_failure = fringe_failure(low_mean)
    high_failure = fringe_failure(high_mean)
    central_failures = 1.0 / (
        1.0 + (CENTRAL_FACTOR - 1.0) ** 2 * means[32:225]
    )
    central_failure = float(np.sum(central_failures))
    rank_failure = sum(math.ldexp(1.0, index - B) for index in range(K))
    total_failure = (
        extreme_failure
        + low_failure
        + high_failure
        + central_failure
        + rank_failure
    )

    uniform = json.loads(UNIFORM_RECEIPT.read_text(encoding="utf-8"))
    uniform_factor = float(uniform["parameters"]["spectrum_factor"])
    adjusted = np.asarray(
        [
            float(row["pointwise_log2_upper"])
            - int(row["active_outer_rows"])
            * math.log2(uniform_factor / CENTRAL_FACTOR)
            for row in uniform["occupation_rows"]
        ]
    )
    q2_plus_log2 = float(logsumexp(adjusted[1:] * math.log(2.0)) / math.log(2.0))
    dominant_q2_plus = int(np.argmax(adjusted[1:])) + 2

    payload = {
        "schema": "single-random-constituent-three-band-event-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_code": "one uniform binary [256,128] generator",
            "extreme_weights": [[1, 24], [232, 256]],
            "low_fringe": [25, 31],
            "central_band": [32, 224],
            "high_fringe": [225, 231],
            "fringe_count_cap": FRINGE_CAP,
            "central_shell_factor": CENTRAL_FACTOR,
        },
        "event_probability": {
            "extreme_markov_failure_upper": extreme_failure,
            "low_fringe_mean": low_mean,
            "low_fringe_cantelli_failure_upper": low_failure,
            "high_fringe_mean": high_mean,
            "high_fringe_cantelli_failure_upper": high_failure,
            "central_cantelli_union_upper": central_failure,
            "rank_failure_upper": rank_failure,
            "total_failure_upper": total_failure,
            "success_probability_lower": 1.0 - total_failure,
            "ten_attempt_abort_log2_upper": 10.0 * math.log2(total_failure),
        },
        "pure_central_distance": {
            "source_uniform_factor": uniform_factor,
            "central_factor": CENTRAL_FACTOR,
            "q2_through_L_log2_upper": q2_plus_log2,
            "q2_through_L_margin_bits": -q2_plus_log2,
            "dominant_occupation": dominant_q2_plus,
            "qL_log2_upper": float(adjusted[-1]),
            "qL_margin_bits": -float(adjusted[-1]),
        },
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The event does not yet imply the distance claim for compositions containing fringe rows.",
            "Efficiently testing the event is not supplied.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
