#!/usr/bin/env python3
"""Verify the exact DP-2Lap one-data suffix bound."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
POPULATION = EXPLORATIONS / "riffle_dp_g4_one_data_support39.json"
CERTIFICATE = (
    EXPLORATIONS / "riffle_dp_2lap_g4_three_node_block_certificate_w8.json"
)
BOUND = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound_verification.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def probability(total: int, suffix: int, support: int) -> Fraction:
    result = Fraction(1)
    for index in range(support):
        result *= Fraction(suffix - index, total - index)
    return result


def main() -> None:
    population = json.loads(POPULATION.read_text())
    certificate = json.loads(CERTIFICATE.read_text())
    bound = json.loads(BOUND.read_text())
    lower = certificate["three_node_weight_lower_bound"]
    distance = certificate["distance_threshold"]
    minimum_prefix = 3 * (distance // lower + 1)
    if minimum_prefix != bound["minimum_zero_prefix_nodes"]:
        raise RuntimeError("suffix verification: prefix derivation mismatch")
    if (minimum_prefix // 3) * lower <= distance:
        raise RuntimeError("suffix verification: insufficient prefix weight")

    counts = {int(key): int(value) for key, value in population["outer_word_counts"].items()}
    rows = {row["packet_support"]: row for row in bound["rows"]}
    total = Fraction()
    for support, count in counts.items():
        value = probability(
            bound["total_packet_positions"],
            bound["suffix_packet_positions"],
            support,
        )
        row = rows[support]
        recorded = Fraction(
            int(row["suffix_placement_probability_numerator"]),
            int(row["suffix_placement_probability_denominator"]),
        )
        if value != recorded:
            raise RuntimeError("suffix verification: row probability mismatch")
        total += count * value
    recorded_total = Fraction(
        int(bound["aggregate_contribution"]["numerator"]),
        int(bound["aggregate_contribution"]["denominator"]),
    )
    if total != recorded_total or total >= Fraction(1, 1 << 40):
        raise RuntimeError("suffix verification: aggregate mismatch")

    payload = {
        "schema": "riffle-dp-2lap-g4-one-data-suffix-bound-verification-v2",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_INDEPENDENT_VERIFICATION",
        "bound_receipt_sha256": digest(BOUND),
        "verified_outer_word_count": sum(counts.values()),
        "verified_supports": sorted(counts),
        "derived_minimum_zero_prefix_nodes": minimum_prefix,
        "certified_prefix_weight": (minimum_prefix // 3) * lower,
        "aggregate_numerator": str(total.numerator),
        "aggregate_denominator": str(total.denominator),
        "below_2^-40": True,
        "method": (
            "Recompute every hypergeometric probability as a product of falling "
            "ratios, then compare the exact aggregate fraction."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"verified_outer_words={sum(counts.values())}")
    print(f"derived_minimum_prefix_nodes={minimum_prefix}")
    print(f"certified_prefix_weight={payload['certified_prefix_weight']}")
    print("below_2^-40=True")
    print(f"output={OUTPUT}")
    print("status=EXACT_INDEPENDENT_VERIFICATION")


if __name__ == "__main__":
    main()
