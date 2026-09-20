#!/usr/bin/env python3
"""Exact one-data bound for the DP-2Lap certified suffix stratum."""

from __future__ import annotations

import hashlib
import json
import math
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
POPULATION = EXPLORATIONS / "riffle_dp_g4_one_data_support39.json"
BLOCK_CERTIFICATE = (
    EXPLORATIONS / "riffle_dp_2lap_g4_three_node_block_certificate_w8.json"
)
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound.json"
TOTAL_NODES = 32_772
MINIMUM_ZERO_PREFIX_NODES = 22_653
SUFFIX_NODES = TOTAL_NODES - MINIMUM_ZERO_PREFIX_NODES
PACKETS_PER_NODE = 16
TOTAL_PACKET_POSITIONS = TOTAL_NODES * PACKETS_PER_NODE
SUFFIX_PACKET_POSITIONS = SUFFIX_NODES * PACKETS_PER_NODE


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    with localcontext() as context:
        context.prec = 100
        numerator = Decimal(value.numerator)
        denominator = Decimal(value.denominator)
        estimate = (numerator.ln() - denominator.ln()) / Decimal(2).ln()
        unit = Decimal(1).scaleb(-places)
        lower = estimate.quantize(unit, rounding=ROUND_FLOOR)
        upper = estimate.quantize(unit, rounding=ROUND_CEILING)
        return [format(lower, "f"), format(upper, "f")]


def placement_probability(support: int) -> Fraction:
    return Fraction(
        math.comb(SUFFIX_PACKET_POSITIONS, support),
        math.comb(TOTAL_PACKET_POSITIONS, support),
    )


def main() -> None:
    population = json.loads(POPULATION.read_text())
    certificate = json.loads(BLOCK_CERTIFICATE.read_text())
    support_counts = {
        int(support): int(count)
        for support, count in population["outer_word_counts"].items()
    }
    if sum(support_counts.values()) != population["record_count"]:
        raise RuntimeError("DP-2Lap suffix bound: population count mismatch")
    derived_prefix = 3 * (
        certificate["distance_threshold"]
        // certificate["three_node_weight_lower_bound"]
        + 1
    )
    if derived_prefix != MINIMUM_ZERO_PREFIX_NODES:
        raise RuntimeError("DP-2Lap suffix bound: certificate geometry changed")

    rows = []
    total = Fraction()
    for support, count in sorted(support_counts.items()):
        probability = placement_probability(support)
        contribution = count * probability
        total += contribution
        rows.append(
            {
                "packet_support": support,
                "outer_word_count": count,
                "suffix_placement_probability_numerator": str(probability.numerator),
                "suffix_placement_probability_denominator": str(probability.denominator),
                "suffix_placement_log2_interval": log2_interval(probability),
                "contribution_numerator": str(contribution.numerator),
                "contribution_denominator": str(contribution.denominator),
                "contribution_log2_interval": log2_interval(contribution),
            }
        )

    payload = {
        "schema": "riffle-dp-2lap-g4-one-data-suffix-bound-v2",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "RIGOROUS_BOUND",
        "population_receipt_sha256": digest(POPULATION),
        "three_node_certificate_sha256": digest(BLOCK_CERTIFICATE),
        "total_nodes": TOTAL_NODES,
        "minimum_zero_prefix_nodes": MINIMUM_ZERO_PREFIX_NODES,
        "suffix_nodes": SUFFIX_NODES,
        "total_packet_positions": TOTAL_PACKET_POSITIONS,
        "suffix_packet_positions": SUFFIX_PACKET_POSITIONS,
        "event": (
            "All active packets of a fixed outer word occupy the final 10119 nodes."
        ),
        "probability_space": (
            "The construction samples one uniform bijection of all 524352 packet "
            "positions. Thus the active positions of a support-h word form a "
            "uniform h-subset."
        ),
        "argument": (
            "For a nonzero first-lap terminal state, the preceding 22653 zero-input "
            "nodes contribute more than the distance threshold on the second lap. "
            "A bad word in this suffix stratum must therefore have terminal state "
            "zero. Ignoring that additional restriction gives the recorded upper bound."
        ),
        "rows": rows,
        "aggregate_contribution": {
            "numerator": str(total.numerator),
            "denominator": str(total.denominator),
            "log2_interval": log2_interval(total),
            "below_2^-40": total < Fraction(1, 1 << 40),
        },
        "scope_limitation": (
            "This row covers only one-data outer words through packet support 39 "
            "and only the placement stratum whose packets all lie in the final "
            "10119 nodes. It does not bound other placements."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"one_data_words={sum(support_counts.values())}")
    print(f"suffix_nodes={SUFFIX_NODES}")
    print(f"aggregate_log2_interval={payload['aggregate_contribution']['log2_interval']}")
    print(f"below_2^-40={payload['aggregate_contribution']['below_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=RIGOROUS_ONE_DATA_SUFFIX_BOUND")


if __name__ == "__main__":
    main()
