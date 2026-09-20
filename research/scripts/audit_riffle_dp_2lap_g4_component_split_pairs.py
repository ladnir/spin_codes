#!/usr/bin/env python3
"""Audit the exact dimension-33 component-split pair prototype."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, binary_rank, build_apply
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import observations
from prepare_riffle_dp_2lap_g4_component_split_words import limbs
from probe_riffle_dp_g4_component_mixing import transpose_columns
from solve_riffle_dp_2lap_g4_component_split_sat import (
    balanced_split,
    component_basis,
    generator_words,
)


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
INPUT = EXPLORATIONS / "riffle_dp_2lap_g4_component_split_words_s53_v01.txt"
RECEIPT = EXPLORATIONS / "riffle_dp_2lap_g4_component_split_pairs_s53_v01.json"
SAT_PROBE = EXPLORATIONS / "riffle_dp_2lap_g4_component_split_sat_s53_v01.json"
SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_split_pairs.cpp"
EXECUTABLE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_split_pairs.exe"
PREPARE = ROOT / "scripts" / "prepare_riffle_dp_2lap_g4_component_split_words.py"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_component_split_pairs_audit.json"
SUPPORT = 0x53
PACKET_VALUE = 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def linear_combination(basis: tuple[int, ...], traversal_index: int) -> int:
    information = traversal_index ^ (traversal_index >> 1)
    result = 0
    for index, character in enumerate(basis):
        if (information >> index) & 1:
            result ^= character
    return result


def main() -> None:
    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    left_mask, right_mask = balanced_split(SUPPORT)
    left_basis = component_basis(transpose_step, left_mask)
    right_basis = component_basis(transpose_step, right_mask)
    if (left_mask, right_mask) != (0x40, 0x13):
        raise RuntimeError("component-split pair audit: split changed")
    if len(left_basis) != 20 or len(right_basis) != 13:
        raise RuntimeError("component-split pair audit: dimensions changed")
    if binary_rank(left_basis + right_basis) != 33:
        raise RuntimeError("component-split pair audit: basis rank changed")

    states = observations(step, PACKET_VALUE, 24)
    left_words = generator_words(states, left_basis)
    right_words = generator_words(states, right_basis)
    expected_lines = [
        "RIFFLE_DP_2LAP_G4_COMPONENT_SPLIT_WORDS_V1",
        "53 1 40 13 20 13",
    ]
    for side, words in (("L", left_words), ("R", right_words)):
        for word in words:
            expected_lines.append(
                side + " " + " ".join(f"{limb:016x}" for limb in limbs(word))
            )
    if INPUT.read_text() != "\n".join(expected_lines) + "\n":
        raise RuntimeError("component-split pair audit: generator input changed")

    receipt = json.loads(RECEIPT.read_text())
    if (
        receipt["schema"] != "riffle-dp-2lap-g4-component-split-pairs-v1"
        or receipt["support_mask_hex"] != "0x53"
        or receipt["packet_value"] != PACKET_VALUE
        or receipt["left_support_mask_hex"] != "0x40"
        or receipt["right_support_mask_hex"] != "0x13"
        or receipt["left_dimension"] != 20
        or receipt["right_dimension"] != 13
        or receipt["pair_count"] != (1 << 33) - 1
        or receipt["minimum_two_sided_weight"] != 120
        or receipt["result"] != "PASS"
    ):
        raise RuntimeError("component-split pair audit: receipt changed")
    left_index = int(receipt["left_gray_index_hex"], 16)
    right_index = int(receipt["right_gray_index_hex"], 16)
    character = linear_combination(left_basis, left_index) ^ linear_combination(
        right_basis, right_index
    )
    weight = sum((character & state).bit_count() & 1 for state in states)
    if weight != receipt["witness_weight"]:
        raise RuntimeError("component-split pair audit: witness weight changed")
    if min(weight, 384 - weight) != receipt["minimum_two_sided_weight"]:
        raise RuntimeError("component-split pair audit: witness side changed")

    sat_probe = json.loads(SAT_PROBE.read_text())
    if sat_probe["result"] != "UNKNOWN" or any(
        row["result"] != "UNKNOWN" for row in sat_probe["rows"]
    ):
        raise RuntimeError("component-split pair audit: SAT diagnostic changed")

    full_pair_count = (1 << 64) - 1
    pair_rate = receipt["pairs_per_second"]
    projected_seconds = full_pair_count / pair_rate
    projected_years = projected_seconds / (365.25 * 24 * 60 * 60)
    payload = {
        "schema": "riffle-dp-2lap-g4-component-split-pairs-audit-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_DIMENSION33_ENUMERATION_WITH_DIAGNOSTIC_FULL64_PROJECTION",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "pair_source_sha256": digest(SOURCE),
        "pair_executable_sha256": digest(EXECUTABLE),
        "prepare_source_sha256": digest(PREPARE),
        "generator_input_sha256": digest(INPUT),
        "pair_receipt_sha256": digest(RECEIPT),
        "sat_probe_sha256": digest(SAT_PROBE),
        "exact_dimension33_result": {
            "support_mask_hex": "0x53",
            "packet_value": PACKET_VALUE,
            "split_masks_hex": ["0x40", "0x13"],
            "split_dimensions": [20, 13],
            "pair_count": receipt["pair_count"],
            "minimum_two_sided_weight": receipt["minimum_two_sided_weight"],
            "elapsed_seconds": receipt["elapsed_seconds"],
            "pairs_per_second": pair_rate,
            "witness_character_hex": hex(character),
            "witness_weight": weight,
            "result": "PASS_EXACT",
        },
        "generic_sat_result": {
            "timeout_seconds_per_side": 10,
            "weight_side": sat_probe["rows"][0]["result"],
            "complement_side": sat_probe["rows"][1]["result"],
            "result": "NO_EXACT_CONCLUSION_WITHIN_BUDGET",
        },
        "full64_naive_projection": {
            "balanced_split_dimensions": [32, 32],
            "pair_count": full_pair_count,
            "projected_seconds_at_dimension33_rate": projected_seconds,
            "projected_years_at_dimension33_rate": projected_years,
            "balanced_right_list_codeword_bytes": 48 * (1 << 32),
            "balanced_right_list_codeword_gibibytes": 192,
            "balanced_right_list_with_64bit_index_gibibytes": 224,
            "unbalanced_44_plus_20_right_list_mibibytes": 56,
            "unbalanced_pair_count": full_pair_count,
            "interpretation": (
                "An unbalanced split fixes memory but not the 2^64 pair count. "
                "The timing projection is diagnostic and optimistic because the "
                "dimension-33 right list fits in cache."
            ),
        },
        "result": "EXACT_SPLIT_VALIDATED_BUT_NAIVE_FULL64_METHOD_INFEASIBLE",
        "scope_limitation": (
            "The dimension-33 minimum and enumeration count are exact. The "
            "full-64 time estimate extrapolates measured throughput and is not "
            "a certificate or hardware-independent prediction."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("dimension33_pair_enumeration=PASS_EXACT")
    print(f"minimum_two_sided_weight={receipt['minimum_two_sided_weight']}")
    print(f"full64_projected_years={projected_years:.6f}")
    print(f"output={OUTPUT}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
