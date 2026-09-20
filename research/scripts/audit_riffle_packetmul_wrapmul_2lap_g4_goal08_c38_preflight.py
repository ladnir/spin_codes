#!/usr/bin/env python3
"""Audit the representative C38 zero-anchor information-set preflight."""

from __future__ import annotations

import hashlib
import json
from math import comb
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns
from probe_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor import (
    MASK64,
    delete_node,
    kernel_basis,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PACKER = ROOT / "scripts" / "prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_information_sets.py"
UNIFORM_SOURCE = ROOT / "scripts" / "prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_uniform35.py"
REFUTATION_SOURCE = ROOT / "scripts" / "probe_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_refutation.py"
REFUTATION_RECEIPT = CANDIDATE / "receipts" / "goal08_c38_refutation_search.json"
OBSERVATION_SOURCE = ROOT / "scripts" / "analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows.py"
STATE_SOURCE = ROOT / "scripts" / "analyze_systematic_group_kernel.py"
GOAL06_AUDIT = CANDIDATE / "receipts" / "goal06_amortized_route_audit.json"
ANCHORS = (0, 19, 37)
OLD_ANCHORS = (0, 12, 23)
NODES = 38
DIMENSION = 64
LENGTH = 37 * 64
BOUND = 227
OUTPUT = CANDIDATE / "receipts" / "goal08_c38_preflight_audit.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def rank(columns: list[int]) -> int:
    pivots: dict[int, int] = {}
    result = 0
    for column in columns:
        value = column
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                result += 1
                break
    return result


def certificate_profile(information_sets: int, bound: int) -> dict:
    base = bound // information_sets
    remainder = bound - base * information_sets
    fixed = remainder + 1
    candidates = information_sets * sum(comb(64, weight) for weight in range(base))
    candidates += fixed * comb(64, base)
    return {
        "information_set_count": information_sets,
        "bound": bound,
        "base_restriction_weight": base,
        "fixed_base_weight_set_count": fixed,
        "candidate_count": candidates,
    }


def reconstruct_columns(anchor: int) -> list[int]:
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    basis = kernel_basis(
        tuple((row >> (64 * anchor)) & MASK64 for row in full_rows),
        64,
    )
    if len(basis) != DIMENSION:
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} has wrong kernel dimension")
    generators = tuple(
        delete_node(xor_selected(full_rows, state), anchor)
        for state in basis
    )
    return [
        sum(
            ((generator >> coordinate) & 1) << bit
            for bit, generator in enumerate(generators)
        )
        for coordinate in range(LENGTH)
    ]


def audit_new_anchor(anchor: int) -> dict:
    path = CANDIDATE / "receipts" / f"goal08_c38_anchor_{anchor:02d}_information_sets.json"
    receipt = load(path)
    if receipt["source_sha256"] != digest(PACKER):
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} source hash changed")
    if receipt["anchor_node"] != anchor or receipt["window_nodes"] != NODES:
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} metadata mismatch")
    if receipt["shortened_code"] != {"length": LENGTH, "dimension": DIMENSION}:
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} code parameters mismatch")
    maximum = receipt["maximum_disjoint_information_sets"]
    sets = receipt["information_sets"]
    if len(sets) != maximum or any(len(row) != DIMENSION for row in sets):
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} base count mismatch")
    flat = [coordinate for row in sets for coordinate in row]
    if len(flat) != len(set(flat)):
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} sets overlap")
    if any(not 0 <= coordinate < LENGTH for coordinate in flat):
        raise RuntimeError(f"Goal 08 audit: anchor {anchor} coordinate out of range")

    columns = reconstruct_columns(anchor)
    for set_index, row in enumerate(sets):
        if rank([columns[coordinate] for coordinate in row]) != DIMENSION:
            raise RuntimeError(
                f"Goal 08 audit: anchor {anchor} set {set_index} is not an information set"
            )

    theoretical = LENGTH // DIMENSION
    if maximum < theoretical:
        successor = receipt["infeasible_successor"]
        if successor is None:
            raise RuntimeError(f"Goal 08 audit: anchor {anchor} lacks successor obstruction")
        if successor["requested_set_count"] != maximum + 1 or successor["success"]:
            raise RuntimeError(f"Goal 08 audit: anchor {anchor} successor metadata mismatch")
        if successor["maximum_union_size"] >= DIMENSION * (maximum + 1):
            raise RuntimeError(f"Goal 08 audit: anchor {anchor} successor is not obstructed")

    profile = certificate_profile(maximum, BOUND)
    radii = {
        str(anchor_weight): (BOUND - anchor_weight) // maximum
        for anchor_weight in range(6)
    }
    return {
        "anchor_node": anchor,
        "receipt_sha256": digest(path),
        "maximum_disjoint_information_sets": maximum,
        "successor_union_size": (
            None
            if receipt["infeasible_successor"] is None
            else receipt["infeasible_successor"]["maximum_union_size"]
        ),
        "systematic_certificate_profile": profile,
        "restriction_radius_by_anchor_weight_zero_through_five": radii,
    }


def old_profile(anchor: int) -> dict:
    path = CANDIDATE / "receipts" / f"goal06_zero_anchor_node_{anchor:02d}_c24.json"
    receipt = load(path)
    return {
        "anchor_node": anchor,
        "receipt_sha256": digest(path),
        "information_set_count": receipt["information_set_count"],
        "base_restriction_weight": receipt["base_restriction_weight"],
        "fixed_base_weight_set_count": receipt["fixed_base_weight_set_count"],
        "candidate_count": receipt["expected_candidate_count"],
    }


def audit_uniform_anchor(anchor: int) -> dict:
    path = CANDIDATE / "receipts" / f"goal08_c38_anchor_{anchor:02d}_uniform35.json"
    receipt = load(path)
    if receipt["source_sha256"] != digest(UNIFORM_SOURCE):
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} source hash changed")
    if receipt["packer_source_sha256"] != digest(PACKER):
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} packer hash changed")
    if receipt["anchor_node"] != anchor or receipt["information_set_count"] != 35:
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} metadata mismatch")
    if receipt["shortened_code"] != {"length": LENGTH, "dimension": DIMENSION}:
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} parameters mismatch")
    sets = receipt["information_sets"]
    if len(sets) != 35 or any(len(row) != DIMENSION for row in sets):
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} set sizes mismatch")
    flat = [coordinate for row in sets for coordinate in row]
    if len(flat) != len(set(flat)):
        raise RuntimeError(f"Goal 08 audit: uniform anchor {anchor} sets overlap")
    columns = reconstruct_columns(anchor)
    for set_index, row in enumerate(sets):
        if rank([columns[coordinate] for coordinate in row]) != DIMENSION:
            raise RuntimeError(
                f"Goal 08 audit: uniform anchor {anchor} set {set_index} is dependent"
            )
    return {
        "anchor_node": anchor,
        "receipt_sha256": digest(path),
        "elapsed_seconds": receipt["elapsed_seconds"],
        "origin_method": receipt["origin"]["method"],
    }


def main() -> None:
    new = [audit_new_anchor(anchor) for anchor in ANCHORS]
    uniform = [audit_uniform_anchor(anchor) for anchor in range(NODES)]
    old = [old_profile(anchor) for anchor in OLD_ANCHORS]
    comparisons = []
    for position, (old_row, new_row) in enumerate(zip(old, new, strict=True)):
        new_profile = new_row["systematic_certificate_profile"]
        comparisons.append({
            "position": ("initial", "center", "terminal")[position],
            "old_anchor": old_row["anchor_node"],
            "new_anchor": new_row["anchor_node"],
            "old_information_sets": old_row["information_set_count"],
            "new_information_sets": new_row["maximum_disjoint_information_sets"],
            "old_base_restriction_weight": old_row["base_restriction_weight"],
            "new_base_restriction_weight": new_profile["base_restriction_weight"],
            "old_candidate_count": old_row["candidate_count"],
            "new_candidate_count": new_profile["candidate_count"],
            "new_to_old_candidate_ratio": (
                new_profile["candidate_count"] / old_row["candidate_count"]
            ),
        })

    total_old = sum(row["candidate_count"] for row in old)
    total_new = sum(row["systematic_certificate_profile"]["candidate_count"] for row in new)
    refutation = load(REFUTATION_RECEIPT)
    if refutation["source_sha256"] != digest(REFUTATION_SOURCE):
        raise RuntimeError("Goal 08 audit: refutation source hash changed")
    parity = build_apply(systematic_state_columns())
    refutation_rows = observation_rows(NODES, parity)
    refutation_state = int(refutation["initial_a_hex"], 16) | (
        int(refutation["initial_b_hex"], 16) << 64
    )
    refutation_word = xor_selected(refutation_rows, refutation_state)
    refutation_weights = [
        ((refutation_word >> (64 * node)) & MASK64).bit_count()
        for node in range(NODES)
    ]
    if refutation_weights != refutation["node_weights"]:
        raise RuntimeError("Goal 08 audit: refutation witness node weights disagree")
    if sum(refutation_weights) != refutation["best_weight"]:
        raise RuntimeError("Goal 08 audit: refutation witness total disagrees")
    if refutation["counterexample_found"] != (refutation["best_weight"] < 228):
        raise RuntimeError("Goal 08 audit: refutation result flag disagrees")
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal08-c38-preflight-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "AUDITED_EXACT_INFORMATION_SET_PREFLIGHT",
        "source_sha256": digest(Path(__file__).resolve()),
        "packer_source_sha256": digest(PACKER),
        "observation_source_sha256": digest(OBSERVATION_SOURCE),
        "state_source_sha256": digest(STATE_SOURCE),
        "goal06_audit_sha256": digest(GOAL06_AUDIT),
        "new_c38_anchors": new,
        "uniform_c38_all_anchors": {
            "anchor_count": len(uniform),
            "information_sets_per_anchor": 35,
            "audited_information_set_count": 35 * len(uniform),
            "anchors": uniform,
            "systematic_certificate_profile": certificate_profile(35, BOUND),
            "restriction_radius_by_anchor_weight_zero_through_five": {
                str(anchor_weight): (BOUND - anchor_weight) // 35
                for anchor_weight in range(6)
            },
            "total_packing_seconds": sum(row["elapsed_seconds"] for row in uniform),
            "maximum_single_anchor_seconds": max(row["elapsed_seconds"] for row in uniform),
        },
        "old_c24_anchors": old,
        "position_comparison": comparisons,
        "representative_total_old_candidates": total_old,
        "representative_total_new_candidates": total_new,
        "representative_new_to_old_candidate_ratio": total_new / total_old,
        "restriction_list_sizes": {
            "radius_6": sum(comb(64, weight) for weight in range(7)),
            "radius_7": sum(comb(64, weight) for weight in range(8)),
            "radius_7_to_radius_6_ratio": (
                sum(comb(64, weight) for weight in range(8))
                / sum(comb(64, weight) for weight in range(7))
            ),
        },
        "diagnostic_refutation_search": {
            "receipt_sha256": digest(REFUTATION_RECEIPT),
            "random_restarts": refutation["random_restarts"],
            "best_weight": refutation["best_weight"],
            "target": refutation["distance_target"],
            "counterexample_found": refutation["counterexample_found"],
            "known_transient_c38_weight": refutation["known_transient_c38_weight"],
        },
        "result": "C38_HAS_THIRTY_FIVE_DISJOINT_INFORMATION_SETS_AT_EVERY_ANCHOR",
        "scope_limitation": (
            "The audit reconstructs and verifies every reported information set and "
            "derives the systematic certificate counts. The exact maximum claim also "
            "depends on the audited matroid-union implementation. No C38 distance bound "
            "or nonzero-anchor certificate is established."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "status": "PASS",
        "information_sets": {
            str(row["anchor_node"]): row["maximum_disjoint_information_sets"]
            for row in new
        },
        "representative_candidate_ratio": payload["representative_new_to_old_candidate_ratio"],
        "output": str(OUTPUT.relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()
