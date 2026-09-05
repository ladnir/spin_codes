#!/usr/bin/env python3
"""Audit both exhaustive C18+C20 24-node certificate families."""

from __future__ import annotations

import hashlib
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
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_c18_c20_24node_audit.json"
PRIMARY_SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_c18_c20_24node.cpp"
PRIMARY_EXE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_c18_c20_24node.exe"
INDEPENDENT_SOURCE = ROOT / "scripts" / "verify_riffle_dp_2lap_g4_c18_c20_24node.cpp"
INDEPENDENT_EXE = ROOT / "scripts" / "verify_riffle_dp_2lap_g4_c18_c20_24node.exe"
FACTOR18 = 0x7C9C3
FACTOR20 = 0x1E1FFF
DIMENSION = 38
WINDOW_NODES = 24
SLOTS = 16
LENGTH = WINDOW_NODES * SLOTS
INFORMATION_SET_COUNT = 10


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def observation_columns(step, packet_value: int) -> tuple[int, ...]:
    columns = []
    for slot in range(SLOTS):
        state = packet_value << (4 * slot)
        for _ in range(WINDOW_NODES):
            state = step(state)
            columns.append(state)
    return tuple(columns)


def code_columns(observations: tuple[int, ...], basis: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(((character & state).bit_count() & 1) << row for row, character in enumerate(basis))
        for state in observations
    )


def load_receipt(packet_value: int, independent: bool) -> tuple[Path, dict]:
    suffix = "_independent" if independent else ""
    path = EXPLORATIONS / (
        f"riffle_dp_2lap_g4_c18_c20_24node_v{packet_value:02d}{suffix}.json"
    )
    return path, json.loads(path.read_text())


def audit_receipt(receipt: dict, columns: tuple[int, ...], packet_value: int, independent: bool) -> None:
    required = 73 if packet_value == 15 else 97
    radius = (required - 1) // INFORMATION_SET_COUNT
    expected_label = (
        "EXACT_INDEPENDENT_EXHAUSTIVE_VERIFICATION"
        if independent
        else "EXACT_EXHAUSTIVE_CERTIFICATE"
    )
    expected_schema = (
        "riffle-dp-2lap-g4-c18-c20-24node-independent-verification-v1"
        if independent
        else "riffle-dp-2lap-g4-c18-c20-24node-certificate-v1"
    )
    expected_ball = sum(math.comb(DIMENSION, weight) for weight in range(radius + 1))
    checks = {
        receipt["schema"] == expected_schema,
        receipt["evidence_label"] == expected_label,
        receipt["packet_value"] == packet_value,
        receipt["code_length"] == LENGTH,
        receipt["code_dimension"] == DIMENSION,
        receipt["component_degrees"] == [18, 20],
        receipt["required_two_sided_distance"] == required,
        receipt["information_radius"] == radius,
        receipt["hamming_ball_size"] == expected_ball,
        receipt["information_set_count"] == INFORMATION_SET_COUNT,
        receipt["total_information_vectors_checked"] == INFORMATION_SET_COUNT * expected_ball,
        receipt["passed"] is True,
        receipt["result"] == "PASS",
    }
    if checks != {True}:
        raise RuntimeError(f"24-node audit: receipt fields failed for value {packet_value}")

    information_sets = receipt["information_sets"]
    if len(information_sets) != INFORMATION_SET_COUNT:
        raise RuntimeError("24-node audit: information-set count changed")
    flat = [coordinate for information_set in information_sets for coordinate in information_set]
    if len(flat) != 380 or len(set(flat)) != 380:
        raise RuntimeError("24-node audit: information sets are not disjoint")
    if set(flat) & set(receipt["unused_coordinates"]):
        raise RuntimeError("24-node audit: unused coordinates overlap a set")
    if set(flat) | set(receipt["unused_coordinates"]) != set(range(LENGTH)):
        raise RuntimeError("24-node audit: coordinate partition is incomplete")
    for information_set in information_sets:
        if len(information_set) != DIMENSION:
            raise RuntimeError("24-node audit: information-set size changed")
        if binary_rank(tuple(columns[index] for index in information_set)) != DIMENSION:
            raise RuntimeError("24-node audit: singular claimed information set")


def main() -> None:
    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(step_columns))
    annihilator = polynomial_multiply(FACTOR18, FACTOR20)
    annihilator_columns = tuple(
        apply_polynomial(transpose_apply, annihilator, 1 << bit)
        for bit in range(64)
    )
    basis = kernel_basis(annihilator_columns)
    if len(basis) != DIMENSION:
        raise RuntimeError("24-node audit: component dimension changed")

    rows = []
    total_primary_candidates = 0
    total_independent_candidates = 0
    for packet_value in range(1, 16):
        columns = code_columns(observation_columns(step, packet_value), basis)
        if binary_rank(columns) != DIMENSION:
            raise RuntimeError("24-node audit: observation code lost dimension")
        primary_path, primary = load_receipt(packet_value, False)
        independent_path, independent = load_receipt(packet_value, True)
        audit_receipt(primary, columns, packet_value, False)
        audit_receipt(independent, columns, packet_value, True)
        if primary["information_sets"] == independent["information_sets"]:
            raise RuntimeError("24-node audit: verifier reused the primary partition")
        total_primary_candidates += primary["total_information_vectors_checked"]
        total_independent_candidates += independent["total_information_vectors_checked"]
        rows.append(
            {
                "packet_value": packet_value,
                "required_two_sided_distance": primary["required_two_sided_distance"],
                "information_radius": primary["information_radius"],
                "primary_receipt_sha256": digest(primary_path),
                "independent_receipt_sha256": digest(independent_path),
                "partitions_are_distinct": True,
                "primary_minimum_enumerated_side": min(
                    primary["minimum_enumerated_weight"],
                    primary["minimum_enumerated_complement_weight"],
                ),
                "independent_minimum_enumerated_side": min(
                    independent["minimum_enumerated_weight"],
                    independent["minimum_enumerated_complement_weight"],
                ),
                "passed": True,
            }
        )

    payload = {
        "schema": "riffle-dp-2lap-g4-c18-c20-24node-audit-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_EXHAUSTIVE_CERTIFICATE_WITH_INDEPENDENT_REPLAY",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "primary_executable_sha256": digest(PRIMARY_EXE),
        "independent_source_sha256": digest(INDEPENDENT_SOURCE),
        "independent_executable_sha256": digest(INDEPENDENT_EXE),
        "component_degrees": [18, 20],
        "code_dimension": DIMENSION,
        "window_nodes": WINDOW_NODES,
        "code_length": LENGTH,
        "information_set_count": INFORMATION_SET_COUNT,
        "packet_value_rows": rows,
        "total_primary_information_vectors_checked": total_primary_candidates,
        "total_independent_information_vectors_checked": total_independent_candidates,
        "result": "PASS",
        "proved_statement": (
            "For every nonzero packet value v and every nonzero character in "
            "C18+C20, the 24-node observation has two-sided weight at least 73 "
            "when v=15 and at least 97 otherwise."
        ),
        "scope_limitation": (
            "The certificate covers only the C18+C20 component subcode. It does "
            "not cover the full 64-dimensional character code or other uncovered "
            "component supports."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"packet_values={len(rows)}")
    print(f"primary_candidates={total_primary_candidates}")
    print(f"independent_candidates={total_independent_candidates}")
    print(f"output={OUTPUT}")
    print("status=PASS_EXACT_EXHAUSTIVE_CERTIFICATE_WITH_INDEPENDENT_REPLAY")


if __name__ == "__main__":
    main()
