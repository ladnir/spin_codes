#!/usr/bin/env python3
"""Independently verify the terminal weight-three G1 ledger."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
EXPLORATIONS = REPOSITORY / "explorations"
MANIFEST = EXPLORATIONS / "riffle_dp_g4_construction_manifest.json"
LEDGER = EXPLORATIONS / "riffle_dp_g4_g1_terminal_ledger.json"
CLASSIFICATION = EXPLORATIONS / "riffle_dp_g4_support38_classification.txt"
OLD_CLASSIFICATION = EXPLORATIONS / "riffle_dp_g4_g1_low_classification.txt"
LOCAL_MESSAGES = REPOSITORY / "bch_g4_le15.bin"
ENGINE = REPOSITORY / "scripts" / "riffle_dp_g4_weighted_triples.py"
CLASSIFIER = REPOSITORY / "scripts" / "classify_riffle_dp_g4_support38.cpp"
ENUMERATOR = REPOSITORY / "scripts" / "enumerate_bch_g4_low_support.cpp"

Q = 16
CODE_SIZE = 1 << 64
TRIPLE_LENGTH = 96
DATA_SYMBOLS = 1 << 14
PACKET_COUNT = 524_352
TERMINAL_PACKETS = 47_184


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def terminal_probability(support: int) -> Fraction:
    return Fraction(
        math.comb(TERMINAL_PACKETS, support),
        math.comb(PACKET_COUNT, support),
    )


def krawtchouk(degree: int, weight: int) -> int:
    return sum(
        (-1) ** h
        * (Q - 1) ** (degree - h)
        * math.comb(weight, h)
        * math.comb(TRIPLE_LENGTH - weight, degree - h)
        for h in range(
            max(0, degree - (TRIPLE_LENGTH - weight)),
            min(degree, weight) + 1,
        )
    )


def exact_dual_bound(minimum_support: int, degree: int) -> Fraction:
    base = terminal_probability(minimum_support)
    k_zero = krawtchouk(degree, 0)
    lines = [
        (
            terminal_probability(weight) / base,
            Fraction(krawtchouk(degree, weight), k_zero),
        )
        for weight in range(minimum_support, TRIPLE_LENGTH + 1)
    ]
    candidates = {Fraction()}
    for left in range(len(lines)):
        left_c, left_r = lines[left]
        for right in range(left + 1, len(lines)):
            right_c, right_r = lines[right]
            if left_r == right_r:
                continue
            candidate = (right_c - left_c) / (left_r - right_r)
            if candidate >= 0:
                candidates.add(candidate)

    best = None
    for dual_weight in candidates:
        intercept = max(c_value + dual_weight * r_value for c_value, r_value in lines)
        objective = (CODE_SIZE - 1) * intercept + dual_weight
        if best is None or objective < best:
            best = objective
    if best is None:
        raise RuntimeError("dual certificate search produced no candidate")
    return base * best


def parse_layers(path: Path) -> tuple[dict[int, dict[str, int]], dict[str, str]]:
    layers = {}
    scalars = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        parsed = {}
        for field in fields:
            if "=" not in field:
                continue
            name, value = field.split("=", 1)
            parsed[name] = value
        if "support" in parsed:
            support = int(parsed.pop("support"))
            layers[support] = {
                name: int(value)
                for name, value in parsed.items()
                if value != "NOT_CLASSIFIED"
            }
        elif len(parsed) == 1:
            scalars.update(parsed)
    return layers, scalars


def fraction_field(value: dict) -> Fraction:
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def verify_local_file() -> None:
    records = np.fromfile(
        LOCAL_MESSAGES,
        dtype=np.dtype([("support", "u1"), ("message", "<u8")], align=False),
    )
    if records.nbytes != LOCAL_MESSAGES.stat().st_size:
        raise RuntimeError("local-message file contains a partial record")
    counts = np.bincount(records["support"], minlength=16)
    expected = {11: 20, 12: 1_526, 13: 37_014, 14: 742_031, 15: 13_430_995}
    for support in range(16):
        if int(counts[support]) != expected.get(support, 0):
            raise RuntimeError(f"unexpected local count at support {support}")


def main() -> None:
    ledger = json.loads(LEDGER.read_text())
    if ledger["decision"] != "PASS_TERMINAL_WEIGHT_THREE":
        raise RuntimeError("ledger does not claim the terminal G1 decision")
    if ledger["manifest_sha256"] != digest(MANIFEST):
        raise RuntimeError("manifest hash mismatch")
    if ledger["engine_source_sha256"] != digest(ENGINE):
        raise RuntimeError("weighted-triple engine hash mismatch")
    bound_sources = ledger["exact_classification"]
    expected_hashes = {
        "classification_sha256": digest(CLASSIFICATION),
        "classifier_source_sha256": digest(CLASSIFIER),
        "local_enumerator_source_sha256": digest(ENUMERATOR),
        "local_message_file_sha256": digest(LOCAL_MESSAGES),
    }
    for field, expected in expected_hashes.items():
        if bound_sources[field] != expected:
            raise RuntimeError(f"bound-source mismatch: {field}")
    verify_local_file()

    layers, scalars = parse_layers(CLASSIFICATION)
    if scalars.get("status") != "EXACT_SUPPORT38_CLASSIFICATION":
        raise RuntimeError("support-38 classifier status is not exact")
    old_layers, old_scalars = parse_layers(OLD_CLASSIFICATION)
    if old_scalars.get("status") != "EXACT_ALL_FOUR_CATEGORIES":
        raise RuntimeError("independent low classifier status is not exact")
    for support in range(33, 36):
        for category in ("one", "pair_p0", "pair_p1", "triple"):
            field = f"{category}_new_dirty_classes"
            if layers[support][field] != old_layers[support][field]:
                raise RuntimeError(
                    f"support-{support} classifiers disagree for {category}"
                )
    if layers[35]["additive_unordered_triples"] * 6 != int(
        old_scalars["triple_additive_ordered_candidates"]
    ):
        raise RuntimeError("ordered and unordered additive receipts disagree")

    repeated_classes = DATA_SYMBOLS + math.comb(DATA_SYMBOLS, 2)
    repeated_classified = sum(
        layers[support].get("one_new_dirty_classes", 0)
        + layers[support].get("pair_p1_new_dirty_classes", 0)
        for support in range(33, 37)
    )
    bucket_counts = {
        33: sum(
            layers[33][f"{category}_new_dirty_classes"]
            for category in ("one", "pair_p1", "pair_p0", "triple")
        ),
        35: sum(
            layers[35][f"{category}_new_dirty_classes"]
            for category in ("one", "pair_p1", "pair_p0", "triple")
        ),
        36: sum(
            layers[36][f"{category}_new_dirty_classes"]
            for category in ("one", "pair_p1", "pair_p0", "triple")
        ),
        37: repeated_classes - repeated_classified,
        38: layers[38]["pair_p0_new_dirty_classes"]
        + layers[38]["triple_new_dirty_classes"],
    }
    coefficient_classes = math.comb(DATA_SYMBOLS + 2, 3)
    bucket_counts[39] = coefficient_classes - sum(bucket_counts.values())
    expected_buckets = {
        33: 3,
        35: 11,
        36: 32_779,
        37: 134_225_892,
        38: 49_143,
        39: 733_007_667_212,
    }
    if bucket_counts != expected_buckets:
        raise RuntimeError(f"unexpected proof partition: {bucket_counts}")

    degrees = {33: 45, 35: 43, 36: 42, 37: 41, 38: 40, 39: 39}
    contributions = {
        support: count * exact_dual_bound(support, degrees[support])
        for support, count in bucket_counts.items()
    }
    upper_bound = sum(contributions.values(), Fraction())
    recorded_upper = fraction_field(ledger["proof_upper_bound"])
    if upper_bound != recorded_upper:
        raise RuntimeError("exact proof upper bound does not reproduce")
    if upper_bound > Fraction(1, 1 << 40):
        raise RuntimeError("terminal G1 upper bound exceeds 2^-40")
    for support, contribution in contributions.items():
        recorded_count = int(ledger["proof_partition"][str(support)]["class_count"])
        if recorded_count != bucket_counts[support] or contribution <= 0:
            raise RuntimeError(f"invalid proof bucket at support {support}")

    exact_lower = 26 * terminal_probability(33) + 39 * terminal_probability(35)
    if exact_lower != fraction_field(ledger["refutation_lower_family"]):
        raise RuntimeError("exact refutation family does not reproduce")
    if exact_lower > upper_bound:
        raise RuntimeError("lower family exceeds proof upper bound")

    print("candidate=Riffle DP g=4@g0-v1")
    print("decision=PASS_TERMINAL_WEIGHT_THREE")
    print("proof_upper_is_below_2^-40=true")
    print("status=EXACT_G1_LEDGER_VERIFIED")


if __name__ == "__main__":
    main()
