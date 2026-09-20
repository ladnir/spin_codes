#!/usr/bin/env python3
"""Audit the corrected dimension-at-most-30 endpoint certificate family."""

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
PRIMARY_RUN = EXPLORATIONS / "riffle_dp_2lap_g4_component_endpoint_primary_run.json"
INDEPENDENT_RUN = EXPLORATIONS / "riffle_dp_2lap_g4_component_endpoint_independent_run.json"
OLD_C18_C20_AUDIT = EXPLORATIONS / "riffle_dp_2lap_g4_c18_c20_24node_audit.json"
TIGHT_WITNESS_AUDIT = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_component24_s49_v15_counterexample_verification.json"
)
OUTPUT = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_component_endpoint_dimension30_audit.json"
)
SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_24node.cpp"
RUNNER = ROOT / "scripts" / "run_riffle_dp_2lap_g4_component_24node.py"

COMPONENT_DEGREES = (1, 2, 4, 9, 10, 18, 20)
COMPONENT_FACTORS = (0x3, 0x7, 0x13, 0x373, 0x519, 0x7C9C3, 0x1E1FFF)
MAXIMAL_MASKS = (0x50, 0x49, 0x31, 0x32, 0x27, 0x47, 0x2B, 0x1F)
SLOTS = 16
TIGHT_CHARACTER = 0xA685AAC60E5ACFB3


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


def observations(step, packet_value: int, window_nodes: int) -> tuple[int, ...]:
    result = []
    for slot in range(SLOTS):
        state = packet_value << (4 * slot)
        for _ in range(window_nodes):
            state = step(state)
            result.append(state)
    return tuple(result)


def code_columns(states: tuple[int, ...], basis: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(
            (((character & state).bit_count() & 1) << row)
            for row, character in enumerate(basis)
        )
        for state in states
    )


def expected_threshold(window_nodes: int, packet_value: int) -> int:
    if window_nodes == 24:
        return 72 if packet_value == 15 else 97
    if window_nodes == 12 and packet_value == 15:
        return 36
    raise RuntimeError("endpoint audit: unexpected certificate case")


def audit_receipt(receipt: dict, columns: tuple[int, ...], mode: str) -> None:
    dimension = receipt["code_dimension"]
    window_nodes = receipt["window_nodes"]
    length = SLOTS * window_nodes
    required = expected_threshold(window_nodes, receipt["packet_value"])
    information_set_count = length // dimension
    radius = (required - 1) // information_set_count
    ball = sum(math.comb(dimension, weight) for weight in range(radius + 1))
    expected_mode = (
        "primary_split_table" if mode == "primary" else "independent_recursive"
    )
    if receipt["schema"] != "riffle-dp-2lap-g4-component-endpoint-certificate-v2":
        raise RuntimeError("endpoint audit: receipt schema changed")
    if receipt["verification_mode"] != expected_mode:
        raise RuntimeError("endpoint audit: verification mode changed")
    if receipt["code_length"] != length:
        raise RuntimeError("endpoint audit: code length changed")
    if receipt["required_two_sided_distance"] != required:
        raise RuntimeError("endpoint audit: threshold changed")
    if receipt["information_set_count"] != information_set_count:
        raise RuntimeError("endpoint audit: information-set count changed")
    if receipt["information_radius"] != radius or receipt["hamming_ball_size"] != ball:
        raise RuntimeError("endpoint audit: Hamming ball changed")
    if receipt["result"] != "PASS" or not receipt["passed"]:
        raise RuntimeError("endpoint audit: non-passing receipt")
    if receipt["total_information_vectors_checked"] != information_set_count * ball:
        raise RuntimeError("endpoint audit: exhaustive candidate count changed")

    information_sets = receipt["information_sets"]
    flat = [coordinate for information_set in information_sets for coordinate in information_set]
    unused = receipt["unused_coordinates"]
    if len(information_sets) != information_set_count:
        raise RuntimeError("endpoint audit: partition count changed")
    if len(flat) != information_set_count * dimension or len(flat) != len(set(flat)):
        raise RuntimeError("endpoint audit: information sets are not disjoint")
    if set(flat) & set(unused) or set(flat) | set(unused) != set(range(length)):
        raise RuntimeError("endpoint audit: partition does not cover the code coordinates")
    for information_set in information_sets:
        if binary_rank(tuple(columns[index] for index in information_set)) != dimension:
            raise RuntimeError("endpoint audit: singular listed information set")


def main() -> None:
    target_masks = tuple(
        mask
        for mask in range(1, 1 << len(COMPONENT_DEGREES))
        if support_dimension(mask) <= 30
    )
    if len(target_masks) != 57:
        raise RuntimeError("endpoint audit: target support count changed")
    if any(not any(mask & ~maximal == 0 for maximal in MAXIMAL_MASKS) for mask in target_masks):
        raise RuntimeError("endpoint audit: maximal-support coverage gap")
    maximal_target_masks = tuple(
        mask
        for mask in target_masks
        if not any(mask != other and mask & ~other == 0 for other in target_masks)
    )
    if set(maximal_target_masks) != set(MAXIMAL_MASKS):
        raise RuntimeError("endpoint audit: maximal-support list changed")

    runs = {
        "primary": json.loads(PRIMARY_RUN.read_text()),
        "independent": json.loads(INDEPENDENT_RUN.read_text()),
    }
    row_maps = {}
    for mode, run in runs.items():
        if (
            run["schema"] != "riffle-dp-2lap-g4-component-endpoint-run-v2"
            or run["mode"] != mode
            or run["result"] != "PASS"
            or len(run["rows"]) != 128
            or run["covered_support_count_by_completed_certificates"] != 57
        ):
            raise RuntimeError("endpoint audit: aggregate run changed")
        if run["runner_source_sha256"] != digest(RUNNER):
            raise RuntimeError("endpoint audit: runner hash mismatch")
        if run["certificate_source_sha256"] != digest(SOURCE):
            raise RuntimeError("endpoint audit: certificate source hash mismatch")
        for key, expected_hash in run["executable_sha256_by_dimension"].items():
            window, dimension = key.split("_d")
            stem = "certify" if mode == "primary" else "verify"
            binary = ROOT / "scripts" / (
                f"{stem}_riffle_dp_2lap_g4_component_{window}_d{dimension}.exe"
            )
            if digest(binary) != expected_hash:
                raise RuntimeError("endpoint audit: executable hash mismatch")
        row_maps[mode] = {
            (
                row["window_nodes"],
                int(row["component_support_mask_hex"], 16),
                row["packet_value"],
            ): row
            for row in run["rows"]
        }

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(step_columns))
    state_cache = {
        (window, value): observations(step, value, window)
        for window, values in ((24, range(1, 16)), (12, (15,)))
        for value in values
    }

    rows = []
    total_audited_receipts = 0
    for mask in MAXIMAL_MASKS:
        annihilator = support_annihilator(mask)
        basis = kernel_basis(
            tuple(
                apply_polynomial(transpose_apply, annihilator, 1 << bit)
                for bit in range(64)
            )
        )
        if len(basis) != support_dimension(mask):
            raise RuntimeError("endpoint audit: component kernel dimension changed")
        for window, values in ((24, range(1, 16)), (12, (15,))):
            for value in values:
                identity = (window, mask, value)
                columns = code_columns(state_cache[(window, value)], basis)
                receipts = {}
                for mode in ("primary", "independent"):
                    row = row_maps[mode][identity]
                    path = ROOT / row["receipt"]
                    if digest(path) != row["receipt_sha256"]:
                        raise RuntimeError("endpoint audit: receipt hash mismatch")
                    receipt = json.loads(path.read_text())
                    if int(receipt["component_support_mask_hex"], 16) != mask:
                        raise RuntimeError("endpoint audit: receipt support mismatch")
                    if tuple(receipt["component_degrees"]) != support_degrees(mask):
                        raise RuntimeError("endpoint audit: receipt degree list mismatch")
                    if int(receipt["component_annihilator_hex"], 16) != annihilator:
                        raise RuntimeError("endpoint audit: receipt annihilator mismatch")
                    audit_receipt(receipt, columns, mode)
                    receipts[mode] = receipt
                    total_audited_receipts += 1
                distinct = (
                    receipts["primary"]["information_sets"]
                    != receipts["independent"]["information_sets"]
                )
                if not distinct:
                    raise RuntimeError("endpoint audit: primary and independent partitions coincide")
                rows.append(
                    {
                        "window_nodes": window,
                        "component_support_mask_hex": hex(mask),
                        "component_degrees": list(support_degrees(mask)),
                        "packet_value": value,
                        "required_two_sided_distance": expected_threshold(window, value),
                        "primary_receipt_sha256": row_maps["primary"][identity]["receipt_sha256"],
                        "independent_receipt_sha256": row_maps["independent"][identity]["receipt_sha256"],
                        "partitions_are_distinct": distinct,
                        "passed": True,
                    }
                )

    factor_one_residual = apply_polynomial(transpose_apply, COMPONENT_FACTORS[0], TIGHT_CHARACTER)
    tight_weights = {}
    for window in (12, 24):
        tight_weights[str(window)] = sum(
            (TIGHT_CHARACTER & state).bit_count() & 1
            for state in state_cache[(window, 15)]
        )
    if factor_one_residual != 0 or tight_weights != {"12": 36, "24": 72}:
        raise RuntimeError("endpoint audit: degree-one tight witness changed")

    old_audit = json.loads(OLD_C18_C20_AUDIT.read_text())
    if old_audit["result"] != "PASS" or old_audit["code_dimension"] != 38:
        raise RuntimeError("endpoint audit: C18+C20 certificate changed")
    witness_audit = json.loads(TIGHT_WITNESS_AUDIT.read_text())
    if (
        witness_audit["status"] != "PASS_EXACT_COUNTEREXAMPLE_REPLAY"
        or not witness_audit["full_orbit"]["equals_required_value15_cap"]
    ):
        raise RuntimeError("endpoint audit: full-orbit tight-witness replay changed")

    support_ledger = [
        {
            "component_support_mask_hex": hex(mask),
            "component_degrees": list(support_degrees(mask)),
            "dimension": support_dimension(mask),
            "values_1_through_14_24node_distance": 97,
            "value_15_24node_distance": 72,
            "value_15_12node_endpoint_distance": 36,
            "status": "CERTIFIED_BY_MAXIMAL_SUBCODE_INCLUSION",
        }
        for mask in target_masks
    ]
    combined_masks = set(target_masks) | {0x60}
    payload = {
        "schema": "riffle-dp-2lap-g4-component-endpoint-dimension30-audit-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_EXHAUSTIVE_CERTIFICATE_WITH_INDEPENDENT_VERIFICATION",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "certificate_source_sha256": digest(SOURCE),
        "primary_run_sha256": digest(PRIMARY_RUN),
        "independent_run_sha256": digest(INDEPENDENT_RUN),
        "old_c18_c20_audit_sha256": digest(OLD_C18_C20_AUDIT),
        "tight_witness_audit_sha256": digest(TIGHT_WITNESS_AUDIT),
        "target_maximum_dimension": 30,
        "target_support_count": len(target_masks),
        "maximal_support_masks_hex": [hex(mask) for mask in MAXIMAL_MASKS],
        "certificate_case_count_per_mode": len(rows),
        "audited_receipt_count": total_audited_receipts,
        "information_vectors_checked_per_mode": runs["primary"]["total_information_vectors_checked"],
        "information_vectors_checked_both_modes": sum(
            run["total_information_vectors_checked"] for run in runs.values()
        ),
        "certificate_rows": rows,
        "support_ledger": support_ledger,
        "tight_degree_one_witness": {
            "character_hex": hex(TIGHT_CHARACTER),
            "factor_one_residual_hex": hex(factor_one_residual),
            "packet_value": 15,
            "weight_24_nodes": tight_weights["24"],
            "weight_12_nodes": tight_weights["12"],
            "full_orbit_absolute_bias": witness_audit["full_orbit"]["absolute_bias"],
            "interpretation": "Both corrected value-15 local bounds and the global 5/8 cap are tight.",
        },
        "global_implications": {
            "full_nodes": 32772,
            "complete_24_node_blocks": 1365,
            "remaining_nodes": 12,
            "values_1_through_14_minimum_weight": 1365 * 97,
            "value_15_minimum_weight": 1365 * 72 + 36,
            "required_three_sixteenths_weight": 98316,
            "value_15_absolute_bias_cap": 0.625,
        },
        "merged_component_coverage": {
            "all_nonempty_component_supports": 127,
            "dimension_at_most_30_supports": len(target_masks),
            "additional_c18_c20_support_mask_hex": "0x60",
            "distinct_certified_supports": len(combined_masks),
            "remaining_supports": 127 - len(combined_masks),
        },
        "result": "PASS",
        "scope_limitation": (
            "The exact result covers the 57 supports of total dimension at most 30. "
            "The earlier C18+C20 certificate adds support 0x60. The other 69 "
            "nonempty component supports remain outside these finite certificates."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"target_supports={len(target_masks)}")
    print(f"certificate_cases_per_mode={len(rows)}")
    print(f"audited_receipts={total_audited_receipts}")
    print(f"combined_certified_supports={len(combined_masks)}")
    print(f"output={OUTPUT}")
    print("status=PASS")


if __name__ == "__main__":
    main()
