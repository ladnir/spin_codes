#!/usr/bin/env python3
"""Independent recursive audit of the Riffle ParallelAcc g=4 Goal 01 receipt."""

from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
RECEIPT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal01_reduced_exact.json"
OUTPUT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal01_audit.json"


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    result = []
    for packet in range(32):
        value = (word >> (4 * packet)) & 15
        if value:
            result.append(value)
    return result


def independent_bad_count(length: int, counts: Counter[int], threshold: int) -> int:
    values = tuple(sorted(counts))
    initial = tuple(counts[value] for value in values)

    @lru_cache(maxsize=None)
    def visit(position: int, state: int, remaining: tuple[int, ...], budget: int) -> int:
        if budget < 0:
            return 0
        if position == length:
            return int(not any(remaining))

        active_remaining = sum(remaining)
        zero_remaining = length - position - active_remaining
        ways = 0
        if zero_remaining:
            ways += visit(position + 1, state, remaining, budget - state.bit_count())
        for index, value in enumerate(values):
            if not remaining[index]:
                continue
            following = list(remaining)
            following[index] -= 1
            next_state = state ^ value
            ways += visit(
                position + 1,
                next_state,
                tuple(following),
                budget - next_state.bit_count(),
            )
        return ways

    return visit(0, 0, initial, threshold)


def total_strings(length: int, counts: Counter[int]) -> int:
    result = math.factorial(length) // math.factorial(length - sum(counts.values()))
    for count in counts.values():
        result //= math.factorial(count)
    return result


def main() -> None:
    receipt = json.loads(RECEIPT.read_text())
    source = json.loads(SOURCE.read_text())
    first_word = source["receipts"][0]["local_words"][0]["codeword_hex"]
    profiles = {
        "authenticated_local_support11": Counter(packet_values(first_word)),
        "four_basis_pairs": Counter({1: 2, 2: 2, 4: 2, 8: 2}),
        "four_equal_pairs": Counter({3: 2, 5: 2, 6: 2, 7: 2}),
        "eight_equal_full_packets": Counter({15: 8}),
    }

    length = receipt["experiment"]["packet_positions"]
    threshold = receipt["experiment"]["failure_weight_inclusive"]
    checks = []
    for row in receipt["profiles"]:
        name = row["profile"]
        counts = profiles[name]
        bad = independent_bad_count(length, counts, threshold)
        total = total_strings(length, counts)
        recorded = row["packet_permutation_bad_probability"]
        expected = Fraction(int(recorded["numerator"]), int(recorded["denominator"]))
        actual = Fraction(bad, total)
        checks.append(
            {
                "profile": name,
                "match": actual == expected,
                "bad_strings": str(bad),
                "total_strings": str(total),
            }
        )

    passed = all(check["match"] for check in checks)
    payload = {
        "schema": "riffle-parallelacc-g4-goal01-audit-v1",
        "method": "independent top-down recursion over remaining packet multiplicities",
        "checks": checks,
        "status": "PASS" if passed else "FAIL",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"profiles_checked={len(checks)}")
    print(f"output={OUTPUT}")
    print(f"status={payload['status']}")
    if not passed:
        raise RuntimeError("Goal 01 independent audit failed")


if __name__ == "__main__":
    main()
