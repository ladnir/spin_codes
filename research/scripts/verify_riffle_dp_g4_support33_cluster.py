#!/usr/bin/env python3
"""Independently replay the support-33 three-node cluster witness."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import apply, systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROBE = ROOT / "explorations" / "riffle_dp_g4_g2_support33_cluster_probe_l5940.json"
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_support33_cluster_family.json"
OUTER_RECEIPTS = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
SAMPLE_SPACE = 524_352
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def receipt_values(receipt: dict) -> list[int]:
    values = []
    for local_word in receipt["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        values.extend(
            (word >> (4 * slot)) & 0xF
            for slot in range(32)
            if ((word >> (4 * slot)) & 0xF) != 0
        )
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, default=DEFAULT_PROBE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    probe = json.loads(args.probe.read_text())
    witness = probe["best_row"]
    if witness["margin_above_distance"] > 0:
        raise RuntimeError("cluster verifier: probe has no sub-threshold witness")
    outer = json.loads(OUTER_RECEIPTS.read_text())
    receipts = {row["family_id"]: row for row in outer["receipts"]}
    expected = Counter(receipt_values(receipts[witness["family_ids"][0]]))
    if sum(expected.values()) != 33:
        raise RuntimeError("cluster verifier: outer profile support mismatch")
    for family_id in witness["family_ids"]:
        if Counter(receipt_values(receipts[family_id])) != expected:
            raise RuntimeError("cluster verifier: family profile mismatch")

    node_words = [int(value, 16) for value in witness["node_words_hex"]]
    actual = Counter(
        (word >> (4 * slot)) & 0xF
        for word in node_words
        for slot in range(16)
        if ((word >> (4 * slot)) & 0xF) != 0
    )
    if actual != expected:
        raise RuntimeError("cluster verifier: node words do not match outer values")

    p_columns = systematic_state_columns()
    state = 0
    total_weight = 0
    prefix_weights = []
    for index in range(probe["window_length"]):
        drive = node_words[index] if index < len(node_words) else 0
        emitted = accumulate(state ^ drive)
        total_weight += emitted.bit_count()
        prefix_weights.append(total_weight)
        state = apply(p_columns, emitted)
    if total_weight != witness["best_weight"]:
        raise RuntimeError("cluster verifier: direct replay weight mismatch")
    if total_weight > DISTANCE_THRESHOLD:
        raise RuntimeError("cluster verifier: witness exceeds the distance threshold")

    translations = probe["window_length"] - len(node_words) + 1
    if witness["valid_terminal_translation_count_lower"] != translations:
        raise RuntimeError("cluster verifier: terminal translation count mismatch")
    if any(prefix_weights[length - 1] > DISTANCE_THRESHOLD for length in range(3, probe["window_length"] + 1)):
        raise RuntimeError("cluster verifier: later translation is not sub-threshold")

    label_multiplicity = math.prod(math.factorial(count) for count in expected.values())
    falling = math.prod(range(SAMPLE_SPACE - 32, SAMPLE_SPACE + 1))
    event_probability = Fraction(translations * label_multiplicity, falling)
    log2_probability = math.log2(event_probability.numerator) - math.log2(event_probability.denominator)
    aggregate_contribution = len(witness["family_ids"]) * event_probability
    log2_aggregate = (
        math.log2(aggregate_contribution.numerator)
        - math.log2(aggregate_contribution.denominator)
    )
    payload = {
        "schema": "riffle-dp-g4-g2-support33-cluster-family-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "probe_sha256": digest(args.probe),
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "verifier_source_sha256": digest(Path(__file__).resolve()),
        "window_length": probe["window_length"],
        "weight": total_weight,
        "distance_threshold": DISTANCE_THRESHOLD,
        "margin_above_distance": total_weight - DISTANCE_THRESHOLD,
        "terminal_translation_count": translations,
        "compatible_outer_family_ids": witness["family_ids"],
        "packet_value_histogram": dict(sorted(expected.items())),
        "node_words_hex": witness["node_words_hex"],
        "single_family_event_probability": {
            "numerator": str(event_probability.numerator),
            "denominator": str(event_probability.denominator),
            "log2": log2_probability,
        },
        "aggregate_first_moment_contribution": {
            "numerator": str(aggregate_contribution.numerator),
            "denominator": str(aggregate_contribution.denominator),
            "log2": log2_aggregate,
        },
        "scope_limitation": (
            "This exact lower family disproves deterministic distance. Its first-"
            "moment contribution is too small to refute the random-permutation claim."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"window_length={probe['window_length']}")
    print(f"weight={total_weight}")
    print(f"margin={total_weight-DISTANCE_THRESHOLD}")
    print(f"terminal_translations={translations}")
    print(f"compatible_outer_families={len(witness['family_ids'])}")
    print(f"single_family_log2_probability={log2_probability:.12f}")
    print(f"aggregate_log2_first_moment_contribution={log2_aggregate:.12f}")
    print(f"output={args.output}")
    print("status=VERIFIED_EXACT_SUPPORT33_CLUSTER_FAMILY")


if __name__ == "__main__":
    main()
