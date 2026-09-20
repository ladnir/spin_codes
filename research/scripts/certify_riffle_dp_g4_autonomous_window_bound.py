#!/usr/bin/env python3
"""Certify a finite-window lower bound for the fixed Riffle DP g=4 map."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEARCH_RECEIPT = (
    ROOT / "explorations" / "riffle_dp_g4_g2_low_output_returns_w8_gap2.json"
)
SEARCH_SOURCE = ROOT / "scripts" / "analyze_riffle_dp_g4_low_output_returns.cpp"
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_window_bound.json"
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lower_bound(length: int, low_maximum: int, short_horizon: int) -> int:
    high_minimum = low_maximum + 1
    long_gap_minimum = short_horizon + 1
    return (
        high_minimum * length
        - low_maximum
        - low_maximum * ((length - 1) // long_gap_minimum)
    )


def main() -> None:
    receipt = json.loads(SEARCH_RECEIPT.read_text())
    if receipt["candidate"] != "Riffle DP g=4@g0-v1":
        raise RuntimeError("autonomous window: candidate mismatch")
    low_maximum = receipt["low_output_weight_maximum"]
    high_minimum = receipt["high_output_weight_minimum"]
    short_horizon = receipt["maximum_gap"]
    if (low_maximum, high_minimum, short_horizon) != (8, 9, 2):
        raise RuntimeError("autonomous window: unexpected search parameters")
    expected_states = sum(math.comb(64, weight) for weight in range(1, 9))
    if receipt["enumerated_low_states"] != expected_states:
        raise RuntimeError("autonomous window: low-state count mismatch")
    first_counts = receipt["first_return_counts_by_gap"]
    if sum(first_counts) + receipt["no_first_return_through_bound"] != expected_states:
        raise RuntimeError("autonomous window: first-return partition mismatch")
    minimum_weights = receipt["minimum_first_return_segment_weight_by_gap"]
    minimum_excesses = receipt[
        "minimum_first_return_excess_over_high_baseline_by_gap"
    ]
    for gap, (count, weight, excess) in enumerate(
        zip(first_counts, minimum_weights, minimum_excesses),
        start=1,
    ):
        if count == 0:
            if weight is not None or excess is not None:
                raise RuntimeError("autonomous window: empty return gap has a minimum")
            continue
        if weight is None or excess != weight - high_minimum * gap:
            raise RuntimeError("autonomous window: short-segment minimum mismatch")
        if excess < 0:
            raise RuntimeError("autonomous window: short segment falls below baseline")

    first_crossing = next(
        length
        for length in range(1, INNER_NODES + 1)
        if lower_bound(length, low_maximum, short_horizon) > DISTANCE_THRESHOLD
    )
    full_bound = lower_bound(INNER_NODES, low_maximum, short_horizon)
    if first_crossing != 29_807 or full_bound != 207_556:
        raise RuntimeError("autonomous window: derived constants changed")
    payload = {
        "schema": "riffle-dp-g4-g2-autonomous-window-bound-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "RIGOROUS_BOUND",
        "search_receipt_sha256": digest(SEARCH_RECEIPT),
        "search_source_sha256": digest(SEARCH_SOURCE),
        "parameters": {
            "low_output_weight_maximum": low_maximum,
            "high_output_weight_minimum": high_minimum,
            "short_return_horizon": short_horizon,
            "long_return_gap_minimum": short_horizon + 1,
            "enumerated_low_states": expected_states,
        },
        "short_return_certificate": {
            "first_return_counts_by_gap": first_counts,
            "minimum_segment_weight_by_gap": minimum_weights,
            "minimum_excess_over_high_baseline_by_gap": minimum_excesses,
        },
        "theorem": {
            "scope": (
                "Every nonzero zero-input autonomous output segment of length L "
                "under the fixed map Q=Acc P."
            ),
            "lower_bound": "delta(L) >= 9*L - 8 - 8*floor((L-1)/3)",
            "proof_summary": [
                "Call an output low when its Hamming weight is at most 8.",
                "The exhaustive receipt covers all nonzero low outputs.",
                "A first return after one node is impossible.",
                "Every first-return segment of length two has weight at least 19, which is at least 9*2.",
                "A later first-return segment has weight at least 1+9*(g-1)=9*g-8.",
                "The initial high prefix costs at least 9 per node, and the final segment loses at most 8.",
                "At most floor((L-1)/3) complete segments can have length at least three.",
            ],
        },
        "distance_consequences": {
            "distance_threshold": DISTANCE_THRESHOLD,
            "first_length_forced_above_threshold": first_crossing,
            "full_chain_length": INNER_NODES,
            "full_chain_lower_bound": full_bound,
            "full_chain_margin": full_bound - DISTANCE_THRESHOLD,
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"enumerated_low_states={expected_states}")
    print("window_bound=9*L-8-8*floor((L-1)/3)")
    print(f"first_length_above_distance={first_crossing}")
    print(f"full_chain_lower_bound={full_bound}")
    print(f"full_chain_margin={full_bound - DISTANCE_THRESHOLD}")
    print(f"output={OUTPUT}")
    print("status=RIGOROUS_G2_AUTONOMOUS_WINDOW_BOUND")


if __name__ == "__main__":
    main()
