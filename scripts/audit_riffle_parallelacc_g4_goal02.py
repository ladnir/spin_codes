#!/usr/bin/env python3
"""Independent fixed-parameter audit of the worst Goal 02 moment row."""

from __future__ import annotations

import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
RECEIPT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal02_support33_moment.json"
OUTPUT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal02_audit.json"
PACKET_COUNT = 524352
DISTANCE = 188766


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    return [
        value
        for packet in range(32)
        if (value := (word >> (4 * packet)) & 15) != 0
    ]


def main() -> None:
    receipt = json.loads(RECEIPT.read_text())
    source = json.loads(SOURCE.read_text())
    family_to_counts = {}
    for outer in source["receipts"]:
        values = []
        for local_word in outer["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        family_to_counts[outer["family_id"]] = Counter(values)

    worst = max(receipt["profiles"], key=lambda row: row["log2_probability_upper_bound"])
    family = worst["outer_families"][0]
    counts = family_to_counts[family]
    values = tuple(sorted(counts))
    initial = tuple(counts[value] for value in values)
    support = sum(initial)
    a = worst["a_equals_minus_log_z"]
    b = worst["b_equals_minus_log_tau"]
    tau = math.exp(-b)
    factors = tuple(
        math.exp(-a * weight) / (1.0 - math.exp(-b - a * weight))
        for weight in range(5)
    )

    @lru_cache(maxsize=None)
    def suffix(remaining: tuple[int, ...], state: int) -> float:
        if not any(remaining):
            return 1.0
        total = 0.0
        for index, value in enumerate(values):
            if not remaining[index]:
                continue
            following = list(remaining)
            following[index] -= 1
            next_state = state ^ value
            total += (
                remaining[index]
                * factors[next_state.bit_count()]
                * suffix(tuple(following), next_state)
            )
        return total

    tilted_order_sum = suffix(initial, 0)
    log_choose = (
        math.lgamma(PACKET_COUNT + 1)
        - math.lgamma(support + 1)
        - math.lgamma(PACKET_COUNT - support + 1)
    )
    objective = (
        DISTANCE * a
        + (PACKET_COUNT - support) * b
        - log_choose
        - math.log1p(-tau)
        - math.lgamma(support + 1)
        + math.log(tilted_order_sum)
    )
    reconstructed = objective / math.log(2)
    recorded = worst["log2_probability_upper_bound"]
    difference = abs(reconstructed - recorded)
    passed = difference < 1e-7
    payload = {
        "schema": "riffle-parallelacc-g4-goal02-audit-v1",
        "method": "independent top-down recursion at the recorded parameters",
        "profile_id": worst["profile_id"],
        "outer_family": family,
        "recorded_log2_bound": recorded,
        "reconstructed_log2_bound": reconstructed,
        "absolute_difference": difference,
        "cached_states": suffix.cache_info().currsize,
        "status": "PASS" if passed else "FAIL",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"profile_id={worst['profile_id']}")
    print(f"reconstructed_log2_bound={reconstructed:.12f}")
    print(f"absolute_difference={difference:.3e}")
    print(f"output={OUTPUT}")
    print(f"status={payload['status']}")
    if not passed:
        raise RuntimeError("Goal 02 audit failed")


if __name__ == "__main__":
    main()
