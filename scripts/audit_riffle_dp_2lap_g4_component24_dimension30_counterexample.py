#!/usr/bin/env python3
"""Audit the dimension-at-most-30 tranche and its exact counterexample."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    binary_rank,
    build_apply,
    kernel_basis,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
PRIMARY_RUN = EXPLORATIONS / "riffle_dp_2lap_g4_component24_primary_run.json"
VERIFICATION = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_component24_s49_v15_counterexample_verification.json"
)
OUTPUT = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_component24_dimension30_counterexample_ledger.json"
)
COMPONENT_DEGREES = (1, 2, 4, 9, 10, 18, 20)
COMPONENT_FACTORS = (0x3, 0x7, 0x13, 0x373, 0x519, 0x7C9C3, 0x1E1FFF)
MAXIMAL_MASKS = (0x50, 0x49, 0x31, 0x32, 0x27, 0x47, 0x2B, 0x1F)
FULLY_CERTIFIED_MASK = 0x50
PARTIALLY_RUN_MASK = 0x49
WINDOW_NODES = 24
SLOTS = 16
LENGTH = WINDOW_NODES * SLOTS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def support_degrees(mask: int) -> tuple[int, ...]:
    return tuple(
        degree
        for index, degree in enumerate(COMPONENT_DEGREES)
        if (mask >> index) & 1
    )


def support_dimension(mask: int) -> int:
    return sum(support_degrees(mask))


def polynomial_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def support_annihilator(mask: int) -> int:
    result = 1
    for index, factor in enumerate(COMPONENT_FACTORS):
        if (mask >> index) & 1:
            result = polynomial_multiply(result, factor)
    return result


def observation_columns(step, packet_value: int) -> tuple[int, ...]:
    result = []
    for slot in range(SLOTS):
        state = packet_value << (4 * slot)
        for _ in range(WINDOW_NODES):
            state = step(state)
            result.append(state)
    return tuple(result)


def code_columns(observations: tuple[int, ...], basis: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(((character & state).bit_count() & 1) << row for row, character in enumerate(basis))
        for state in observations
    )


def audit_receipt(receipt: dict, columns: tuple[int, ...], expect_pass: bool) -> None:
    dimension = receipt["code_dimension"]
    packet_value = receipt["packet_value"]
    required = 73 if packet_value == 15 else 97
    information_set_count = LENGTH // dimension
    radius = (required - 1) // information_set_count
    ball = sum(math.comb(dimension, weight) for weight in range(radius + 1))
    if receipt["schema"] != "riffle-dp-2lap-g4-component-24node-certificate-v1":
        raise RuntimeError("dimension-30 audit: receipt schema changed")
    if receipt["required_two_sided_distance"] != required:
        raise RuntimeError("dimension-30 audit: threshold changed")
    if receipt["information_set_count"] != information_set_count:
        raise RuntimeError("dimension-30 audit: information-set count changed")
    if receipt["information_radius"] != radius or receipt["hamming_ball_size"] != ball:
        raise RuntimeError("dimension-30 audit: Hamming ball changed")
    if expect_pass:
        if receipt["result"] != "PASS":
            raise RuntimeError("dimension-30 audit: expected passing receipt")
        if receipt["total_information_vectors_checked"] != information_set_count * ball:
            raise RuntimeError("dimension-30 audit: exhaustive count changed")
    elif (
        receipt["result"] != "COUNTEREXAMPLE"
        or receipt["violating_distance"] != 72
        or int(receipt["violating_character_hex"], 16) != 0xA685AAC60E5ACFB3
    ):
        raise RuntimeError("dimension-30 audit: counterexample changed")

    information_sets = receipt["information_sets"]
    flat = [coordinate for information_set in information_sets for coordinate in information_set]
    if len(information_sets) != information_set_count:
        raise RuntimeError("dimension-30 audit: partition count changed")
    if len(flat) != information_set_count * dimension or len(flat) != len(set(flat)):
        raise RuntimeError("dimension-30 audit: sets are not disjoint")
    if set(flat) | set(receipt["unused_coordinates"]) != set(range(LENGTH)):
        raise RuntimeError("dimension-30 audit: partition does not cover coordinates")
    for information_set in information_sets:
        if binary_rank(tuple(columns[index] for index in information_set)) != dimension:
            raise RuntimeError("dimension-30 audit: singular information set")


def main() -> None:
    target_masks = tuple(
        mask
        for mask in range(1, 1 << len(COMPONENT_DEGREES))
        if support_dimension(mask) <= 30
    )
    if len(target_masks) != 57:
        raise RuntimeError("dimension-30 audit: target support count changed")
    if any(not any(mask & ~maximal == 0 for maximal in MAXIMAL_MASKS) for mask in target_masks):
        raise RuntimeError("dimension-30 audit: maximal-support coverage gap")

    run = json.loads(PRIMARY_RUN.read_text())
    if run["result"] != "COUNTEREXAMPLE" or len(run["rows"]) != 30:
        raise RuntimeError("dimension-30 audit: primary stopping point changed")
    row_by_identity = {
        (int(row["component_support_mask_hex"], 16), row["packet_value"]): row
        for row in run["rows"]
    }

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(step_columns))
    observations = {
        value: observation_columns(step, value) for value in range(1, 16)
    }
    audited_receipts = []
    for mask in (FULLY_CERTIFIED_MASK, PARTIALLY_RUN_MASK):
        annihilator = support_annihilator(mask)
        annihilator_columns = tuple(
            apply_polynomial(transpose_apply, annihilator, 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(annihilator_columns)
        if len(basis) != support_dimension(mask):
            raise RuntimeError("dimension-30 audit: component kernel changed")
        for packet_value in range(1, 16):
            row = row_by_identity[(mask, packet_value)]
            path = ROOT / row["receipt"]
            if digest(path) != row["receipt_sha256"]:
                raise RuntimeError("dimension-30 audit: receipt hash mismatch")
            receipt = json.loads(path.read_text())
            columns = code_columns(observations[packet_value], basis)
            audit_receipt(
                receipt,
                columns,
                expect_pass=not (mask == PARTIALLY_RUN_MASK and packet_value == 15),
            )
            audited_receipts.append(row["receipt_sha256"])

    verification = json.loads(VERIFICATION.read_text())
    counterexample_receipt = (
        EXPLORATIONS / "riffle_dp_2lap_g4_component24_primary_s49_v15.json"
    )
    if (
        verification["status"] != "PASS_EXACT_COUNTEREXAMPLE_REPLAY"
        or verification["certificate_receipt_sha256"] != digest(counterexample_receipt)
        or verification["full_orbit"]["absolute_bias"] != 0.625
    ):
        raise RuntimeError("dimension-30 audit: independent replay changed")

    ledger = []
    for mask in target_masks:
        degrees = support_degrees(mask)
        if mask & 1:
            value15_status = "REFUTED_DISTANCE_73_BY_DEGREE_ONE_WITNESS"
        elif mask & ~FULLY_CERTIFIED_MASK == 0:
            value15_status = "CERTIFIED_DISTANCE_AT_LEAST_73"
        else:
            value15_status = "NOT_RUN_AFTER_COUNTEREXAMPLE"

        if mask & ~FULLY_CERTIFIED_MASK == 0 or mask & ~PARTIALLY_RUN_MASK == 0:
            non15_status = "CERTIFIED_DISTANCE_AT_LEAST_97"
        else:
            non15_status = "NOT_RUN_AFTER_COUNTEREXAMPLE"
        ledger.append(
            {
                "component_support_mask_hex": hex(mask),
                "component_degrees": list(degrees),
                "dimension": sum(degrees),
                "values_1_through_14_status": non15_status,
                "value_15_distance_73_status": value15_status,
            }
        )

    payload = {
        "schema": "riffle-dp-2lap-g4-component24-dimension30-ledger-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_COUNTEREXAMPLE_WITH_AUDITED_PARTIAL_LEDGER",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "primary_run_sha256": digest(PRIMARY_RUN),
        "independent_verification_sha256": digest(VERIFICATION),
        "target_maximum_dimension": 30,
        "target_support_count": len(target_masks),
        "planned_maximal_support_masks_hex": [hex(mask) for mask in MAXIMAL_MASKS],
        "fully_certified_strong_threshold_support_count": sum(
            mask & ~FULLY_CERTIFIED_MASK == 0 for mask in target_masks
        ),
        "value15_refuted_support_count": sum(mask & 1 != 0 for mask in target_masks),
        "value15_unrun_support_count": sum(
            mask & 1 == 0 and mask & ~FULLY_CERTIFIED_MASK != 0
            for mask in target_masks
        ),
        "non15_certified_support_count": sum(
            mask & ~FULLY_CERTIFIED_MASK == 0 or mask & ~PARTIALLY_RUN_MASK == 0
            for mask in target_masks
        ),
        "audited_receipt_count": len(audited_receipts),
        "support_ledger": ledger,
        "counterexample": {
            "component_support_mask_hex": "0x49",
            "component_degrees": [1, 9, 20],
            "packet_value": 15,
            "character_hex": "0xa685aac60e5acfb3",
            "observed_24_node_weight": 72,
            "required_24_node_distance": 73,
            "pure_degree_one_character": True,
            "full_orbit_bias": 0.625,
            "interpretation": (
                "The stronger distance-73 lemma is false. The witness attains, "
                "rather than violates, the desired global value-15 cap."
            ),
        },
        "result": "COUNTEREXAMPLE",
        "revised_proof_target": (
            "Use 24-node distance 72 for value 15 together with a distance-36 "
            "bound for the final 12-node block. Retain 24-node distance 97 for "
            "values 1 through 14."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"target_supports={len(target_masks)}")
    print(f"audited_receipts={len(audited_receipts)}")
    print(f"output={OUTPUT}")
    print("status=PASS_EXACT_COUNTEREXAMPLE_WITH_AUDITED_PARTIAL_LEDGER")


if __name__ == "__main__":
    main()
