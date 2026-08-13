#!/usr/bin/env python3
"""Produce dominance-informed incomplete BSP shards for high-dimensional g=8.

The producer combines the frozen base witness atlas with repeatable
supplementary discovery atlases.  Competing vertex-leading affine witnesses
suggest a dominance cut.  The cut is quantized to a primitive integer form;
only that exact integer form defines ownership.  No soundness claim depends
on its proximity to the floating dominance plane.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import produce_packet_group_g8_highdim_adaptive_shards as coordinate
from certify_packet_group_g8_support_shard import (
    enumerate_vertices,
    root_constraints,
    split_constraint,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTARY_SCHEMA = "permute-conv.packet-group-g8-supplementary-witness-atlas.v1"
BATCH_SCHEMA = "permute-conv.packet-group-g8-highdim-dominance-batch.v1"
SPLIT_POLICY = "vertex-leading-quantized-dominance-v1"
CLASSES = 9


def parse_int_list(text: str) -> tuple[int, ...]:
    return coordinate.parse_int_list(text)


def parse_mask_list(text: str) -> tuple[int, ...]:
    return coordinate.parse_mask_list(text)


def load_supplementary_atlas(
    path: Path,
    manifest_digest: str,
    base_digest: str,
) -> tuple[np.ndarray, np.ndarray, tuple[dict[str, Any], ...], dict[str, Any]]:
    digest = coordinate.sha256_path(path)
    artifact = json.loads(path.read_bytes())
    if artifact.get("schema") != SUPPLEMENTARY_SCHEMA:
        raise ValueError(f"unexpected supplementary atlas schema: {path}")
    if artifact.get("status") != "DIAGNOSTIC_BINARY64_SUPPLEMENTARY_WITNESS_ATLAS":
        raise ValueError(f"supplementary atlas has an unexpected status: {path}")
    if artifact.get("manifest_sha256") != manifest_digest:
        raise ValueError(f"supplementary atlas binds another manifest: {path}")
    base_binding = artifact.get("base_atlas_binding", {})
    if base_binding.get("sha256") != base_digest:
        raise ValueError(f"supplementary atlas binds another base atlas: {path}")
    rows = artifact.get("rows")
    bindings = artifact.get("row_bindings")
    if not isinstance(rows, list) or not isinstance(bindings, list) or len(rows) != len(bindings):
        raise ValueError(f"supplementary row catalogue is malformed: {path}")
    source_id = f"supplementary-{digest[:16]}"
    constants = []
    charges = []
    references = []
    declared = {int(row["row"]): row for row in bindings}
    if len(declared) != len(bindings) or set(declared) != set(range(len(rows))):
        raise ValueError(f"supplementary row indices are not canonical: {path}")
    for index, row in enumerate(rows):
        binding = declared.get(index)
        row_digest = coordinate.sha256_bytes(coordinate.canonical_bytes(row))
        if binding is None or binding.get("row_sha256") != row_digest:
            raise ValueError(f"supplementary row digest mismatch at {path}:{index}")
        affine = row.get("affine", {})
        constant = float(affine.get("constant_log2"))
        charge = np.asarray(affine.get("charge_log2"), dtype=np.float64)
        fugacities = np.asarray(row.get("inner", {}).get("fugacities"), dtype=np.float64)
        if (
            charge.shape != (CLASSES,)
            or fugacities.shape != (CLASSES,)
            or not math.isfinite(constant)
            or not np.all(np.isfinite(charge))
            or np.any(fugacities <= 0.0)
        ):
            raise ValueError(f"supplementary affine row is malformed at {path}:{index}")
        constants.append(constant)
        charges.append(charge)
        references.append(
            {
                "source_id": source_id,
                "source_sha256": digest,
                "row": index,
                "row_sha256": row_digest,
            }
        )
    source = {
        "source_id": source_id,
        "path": str(path),
        "sha256": digest,
        "schema": SUPPLEMENTARY_SCHEMA,
        "row_count": len(rows),
        "rows": [
            {
                "row": index,
                "row_sha256": declared[index]["row_sha256"],
            }
            for index in range(len(rows))
        ],
        "row_bindings_sha256": coordinate.sha256_bytes(
            coordinate.canonical_bytes(bindings)
        ),
    }
    return (
        np.asarray(constants, dtype=np.float64),
        np.asarray(charges, dtype=np.float64),
        tuple(references),
        source,
    )


def combined_bank(
    base: coordinate.AtlasBank,
    supplementary: list[tuple[np.ndarray, np.ndarray, tuple[dict[str, Any], ...], dict]],
) -> coordinate.AtlasBank:
    constants = [base.constants]
    charges = [base.charges]
    references = list(base.references)
    seen_refs = {(row["source_sha256"], int(row["row"])) for row in references}
    for local_constants, local_charges, local_references, _source in supplementary:
        constants.append(local_constants)
        charges.append(local_charges)
        for reference in local_references:
            key = (reference["source_sha256"], int(reference["row"]))
            if key in seen_refs:
                raise ValueError("duplicate supplementary witness row")
            seen_refs.add(key)
            references.append(reference)
    return coordinate.AtlasBank(
        constants=np.concatenate(constants),
        charges=np.vstack(charges),
        references=tuple(references),
    )


def profile_as_float(vertex: tuple[Fraction, ...]) -> list[float]:
    return [float(value) for value in vertex]


def serialize_vertex(vertex: tuple[Fraction, ...]) -> list[str]:
    return [f"{value.numerator}/{value.denominator}" for value in vertex]


def integral_vertex(vertex: tuple[Fraction, ...]) -> list[int] | None:
    if any(value.denominator != 1 for value in vertex):
        return None
    return [value.numerator for value in vertex]


def score_vertices(
    bank: coordinate.AtlasBank,
    vertices: tuple[tuple[Fraction, ...], ...],
    total: int,
    uniform_target: float,
    hardening_reserve: float,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    matrix = np.asarray([profile_as_float(vertex) for vertex in vertices], dtype=np.float64)
    normalizations = coordinate.diagnostic_normalization(matrix, total)
    scores = bank.constants[:, None] - bank.charges @ matrix.T - normalizations[None, :]
    maxima = np.max(scores, axis=1)
    witness = int(np.argmin(maxima))
    worst = int(np.argmax(scores[witness]))
    upper = float(scores[witness, worst])
    target_with_reserve = uniform_target - hardening_reserve
    leaders = np.argmin(scores, axis=0)
    worst_vertex = vertices[worst]
    diagnostic = {
        "candidate_status": (
            "PASS_BINARY64_NEEDS_OUTWARD_REPLAY"
            if upper <= target_with_reserve
            else "FAIL_BINARY64_REQUIRES_SPLIT"
        ),
        "candidate_selector": {
            "kind": "witness",
            "witness": bank.references[witness],
        },
        "candidate_witness_row_in_combined_bank": witness,
        "candidate_upper_log2": upper,
        "uniform_profile_target_log2": uniform_target,
        "hardening_reserve_bits": hardening_reserve,
        "target_with_reserve_log2": target_with_reserve,
        "uniform_target_gap_bits": upper - uniform_target,
        "positive_residual_bits": max(0.0, upper - target_with_reserve),
        "worst_vertex_rational": serialize_vertex(worst_vertex),
        "exact_vertex_count": len(vertices),
        "vertex_leading_witness_count": int(len(set(map(int, leaders)))),
        "vertex_sha256": coordinate.sha256_bytes(
            coordinate.canonical_bytes([serialize_vertex(row) for row in vertices])
        ),
        "arithmetic": "binary64-discovery-over-exact-rational-vertices",
    }
    integer_worst = integral_vertex(worst_vertex)
    if integer_worst is not None:
        diagnostic["worst_vertex"] = integer_worst
    return diagnostic, scores, leaders


def canonicalize_integer_cut(
    coefficients: list[int],
    threshold: int,
    active: tuple[int, ...],
    total: int,
) -> tuple[tuple[int, ...], int, bool] | None:
    """Apply the verifier's support-relative canonical integer-cut rule."""

    active_set = set(active)
    coefficients = [
        int(value) if index in active_set else 0
        for index, value in enumerate(coefficients)
    ]
    pivot = max(active)
    pivot_value = coefficients[pivot]
    for index in active:
        coefficients[index] -= pivot_value
    threshold -= pivot_value * total
    if not any(coefficients[index] for index in active):
        return None
    divisor = math.gcd(*(abs(coefficients[index]) for index in active))
    coefficients = [value // divisor for value in coefficients]
    threshold = threshold // divisor
    first = next(coefficients[index] for index in active if coefficients[index])
    swapped = False
    if first < 0:
        coefficients = [-value for value in coefficients]
        threshold = -threshold - 1
        swapped = True
    if coefficients[pivot] != 0 or any(
        coefficients[index] for index in range(CLASSES) if index not in active_set
    ):
        raise AssertionError("support-relative integer-cut gauge failed")
    if max(abs(value).bit_length() for value in coefficients) > 256:
        return None
    if abs(threshold).bit_length() > 256:
        return None
    return tuple(coefficients), threshold, swapped


def quantized_dominance_cut(
    bank: coordinate.AtlasBank,
    left_witness: int,
    right_witness: int,
    active: tuple[int, ...],
    total: int,
    scale: int,
) -> tuple[tuple[int, ...], int, dict[str, Any]] | None:
    difference = bank.charges[left_witness] - bank.charges[right_witness]
    constant = float(bank.constants[left_witness] - bank.constants[right_witness])
    if not np.all(np.isfinite(difference)) or not math.isfinite(constant):
        return None

    def scaled_round(value: float) -> int:
        # Fraction.from_float binds the exact binary64 proposal. Fraction's
        # round uses ties-to-even, so discovery is reproducible across hosts.
        return round(Fraction.from_float(float(value)) * scale)

    def scaled_floor(value: float) -> int:
        exact = Fraction.from_float(float(value)) * scale
        return exact.numerator // exact.denominator

    raw_coefficients = [
        scaled_round(float(value)) if index in active else 0
        for index, value in enumerate(difference)
    ]
    raw_threshold = scaled_floor(constant)
    canonical = canonicalize_integer_cut(
        raw_coefficients, raw_threshold, active, total
    )
    if canonical is None:
        return None
    coefficients, threshold, children_swapped = canonical
    diagnostic = {
        "source_witnesses": [
            bank.references[left_witness],
            bank.references[right_witness],
        ],
        "quantization_scale": str(scale),
        "support_relative_pivot": max(active),
        "common_normalization_cancelled": True,
        "children_swapped_during_canonicalization": children_swapped,
        "raw_coefficients": [str(value) for value in raw_coefficients],
        "raw_threshold": str(raw_threshold),
        "soundness_scope": (
            "the floating dominance plane guides discovery only; exact ownership and "
            "verification use solely the serialized support-relative canonical "
            "primitive integer split"
        ),
    }
    return coefficients, threshold, diagnostic


def choose_dominance_split(
    bank: coordinate.AtlasBank,
    vertices: tuple[tuple[Fraction, ...], ...],
    leaders: np.ndarray,
    active: tuple[int, ...],
    total: int,
    scale: int,
) -> tuple[tuple[int, ...], int, dict[str, Any]] | None:
    unique = sorted(set(map(int, leaders)))
    if len(unique) < 2:
        return None
    best = None
    for left_witness, right_witness in itertools.combinations(unique, 2):
        cut = quantized_dominance_cut(
            bank, left_witness, right_witness, active, total, scale
        )
        if cut is None:
            continue
        coefficients, threshold, diagnostic = cut
        values = [
            sum(Fraction(coefficient) * value for coefficient, value in zip(coefficients, vertex))
            for vertex in vertices
        ]
        left_indices = [index for index, value in enumerate(values) if value <= threshold]
        right_indices = [index for index, value in enumerate(values) if value >= threshold + 1]
        if not left_indices or not right_indices:
            continue
        separated = sum(
            int(leaders[index]) == left_witness for index in left_indices
        ) + sum(int(leaders[index]) == right_witness for index in right_indices)
        reverse = sum(
            int(leaders[index]) == right_witness for index in left_indices
        ) + sum(int(leaders[index]) == left_witness for index in right_indices)
        alignment = max(separated, reverse)
        quality = (min(len(left_indices), len(right_indices)), alignment, -max(len(left_indices), len(right_indices)))
        if best is None or quality > best[0]:
            best = (
                quality,
                coefficients,
                threshold,
                {
                    **diagnostic,
                    "left_vertex_count": len(left_indices),
                    "right_vertex_count": len(right_indices),
                    "leader_alignment_count": alignment,
                },
            )
    if best is None:
        return None
    return best[1], best[2], best[3]


def coordinate_fallback(
    vertices: tuple[tuple[Fraction, ...], ...],
    active: tuple[int, ...],
    total: int,
) -> tuple[tuple[int, ...], int, dict[str, Any]] | None:
    best = None
    for index in active:
        values = sorted(set(vertex[index] for vertex in vertices))
        if len(values) < 2:
            continue
        middle = values[len(values) // 2 - 1]
        threshold = math.floor(middle)
        left = sum(vertex[index] <= threshold for vertex in vertices)
        right = sum(vertex[index] >= threshold + 1 for vertex in vertices)
        if not left or not right:
            continue
        quality = min(left, right)
        if best is None or quality > best[0]:
            coefficients = [0] * CLASSES
            coefficients[index] = 1
            canonical = canonicalize_integer_cut(
                coefficients, threshold, active, total
            )
            if canonical is None:
                continue
            canonical_coefficients, canonical_threshold, swapped = canonical
            best = (
                quality,
                canonical_coefficients,
                canonical_threshold,
                {
                    "fallback": "vertex-balanced-coordinate-v1",
                    "coordinate": index,
                    "children_swapped_during_canonicalization": swapped,
                    "left_vertex_count": left,
                    "right_vertex_count": right,
                },
            )
    return None if best is None else (best[1], best[2], best[3])


def box_upper_record(
    active: tuple[int, ...], lower: tuple[int, ...], upper: tuple[int, ...], total: int
) -> dict[str, Any]:
    record = coordinate.count_record(active, lower, upper, total)
    record["diagnostic_exact_for_coordinate_path"] = False
    record["diagnostic_count_scope"] = (
        "recomputed coordinate-box upper containing the affine BSP node; no exact "
        "leaf-count conservation claim"
    )
    return record


@dataclass
class WorkNode:
    node_id: str
    depth: int
    constraints: tuple
    lower: tuple[int, ...]
    upper: tuple[int, ...]
    count_upper: int
    vertices: tuple[tuple[Fraction, ...], ...]
    diagnostic: dict[str, Any]
    scores: np.ndarray
    leaders: np.ndarray


def unresolved_node(work: WorkNode, reason: str, active: tuple[int, ...], total: int):
    diagnostic = dict(work.diagnostic)
    diagnostic["unresolved_reason"] = reason
    diagnostic["collapsed_candidate_contribution_log2"] = (
        math.log2(work.count_upper) + float(diagnostic["candidate_upper_log2"])
    )
    return {
        "node_id": work.node_id,
        "state": "UNRESOLVED",
        "count": box_upper_record(active, work.lower, work.upper, total),
        "diagnostic": diagnostic,
    }


def make_child(
    node_id: str,
    depth: int,
    constraints: tuple,
    lower: tuple[int, ...],
    upper: tuple[int, ...],
    active: tuple[int, ...],
    total: int,
    bank: coordinate.AtlasBank,
    uniform_target: float,
    reserve: float,
    maximum_systems: int,
) -> WorkNode | None:
    vertices = enumerate_vertices(active, constraints, maximum_systems)
    if not vertices:
        return None
    count_upper = coordinate.exact_box_count(active, lower, upper, total)
    diagnostic, scores, leaders = score_vertices(
        bank, vertices, total, uniform_target, reserve
    )
    return WorkNode(
        node_id,
        depth,
        constraints,
        lower,
        upper,
        count_upper,
        vertices,
        diagnostic,
        scores,
        leaders,
    )


def produce_shard(
    manifest_digest: str,
    census_row: dict[str, Any],
    bank: coordinate.AtlasBank,
    supplementary_sources: list[dict[str, Any]],
    total: int,
    global_count: int,
    probability_target: float,
    max_nodes: int,
    max_depth: int,
    reserve: float,
    dominance_scale: int,
    maximum_systems: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    active = tuple(map(int, census_row["active_classes"]))
    dimension = int(census_row["dimension"])
    if dimension not in (6, 7, 8) or int(census_row["weight_cut_exclusions"]):
        raise ValueError("dominance producer requires an exclusion-free dimension 6--8 root")
    lower = tuple(1 if index in active else 0 for index in range(CLASSES))
    maximum = total - len(active) + 1
    upper = tuple(maximum if index in active else 0 for index in range(CLASSES))
    root_count = coordinate.exact_box_count(active, lower, upper, total)
    if root_count != int(census_row["count"]):
        raise ValueError("root count changed")
    uniform_target = probability_target - math.log2(global_count)
    constraints = root_constraints(active, total, 21)
    root = make_child(
        "r", 0, constraints, lower, upper, active, total, bank, uniform_target,
        reserve, maximum_systems
    )
    if root is None:
        raise ValueError("high-dimensional root is empty")
    nodes: dict[str, dict[str, Any]] = {}
    pending = []
    serial = 0

    def push(work: WorkNode):
        nonlocal serial
        priority = -(
            math.log2(work.count_upper) + float(work.diagnostic["candidate_upper_log2"])
        )
        heapq.heappush(pending, (priority, serial, work))
        serial += 1

    push(root)
    allocated = 1
    split_counts = {"dominance": 0, "coordinate": 0}
    passing = 0
    residual = 0
    while pending:
        _priority, _serial, work = heapq.heappop(pending)
        if work.diagnostic["candidate_status"] == "PASS_BINARY64_NEEDS_OUTWARD_REPLAY":
            nodes[work.node_id] = unresolved_node(
                work, "AWAITING_INDEPENDENT_OUTWARD_REPLAY", active, total
            )
            passing += 1
            continue
        if work.depth >= max_depth or allocated + 2 > max_nodes:
            reason = "MAX_DEPTH_REACHED" if work.depth >= max_depth else "NODE_BUDGET_EXHAUSTED"
            nodes[work.node_id] = unresolved_node(work, reason, active, total)
            residual += 1
            continue
        split = choose_dominance_split(
            bank, work.vertices, work.leaders, active, total, dominance_scale
        )
        split_kind = "dominance"
        if split is None:
            split = coordinate_fallback(work.vertices, active, total)
            split_kind = "coordinate"
        if split is None:
            nodes[work.node_id] = unresolved_node(
                work, "NO_VERTEX_SEPARATING_INTEGER_SPLIT", active, total
            )
            residual += 1
            continue
        coefficients, threshold, split_diagnostic = split
        child_specs = []
        for side, suffix in ((True, "0"), (False, "1")):
            child_constraints = work.constraints + (
                split_constraint(coefficients, threshold, active, side),
            )
            child_lower = list(work.lower)
            child_upper = list(work.upper)
            nonzero = [index for index, value in enumerate(coefficients) if value]
            if len(nonzero) == 1 and coefficients[nonzero[0]] == 1:
                index = nonzero[0]
                if side:
                    child_upper[index] = min(child_upper[index], threshold)
                else:
                    child_lower[index] = max(child_lower[index], threshold + 1)
            child = make_child(
                work.node_id + suffix,
                work.depth + 1,
                child_constraints,
                tuple(child_lower),
                tuple(child_upper),
                active,
                total,
                bank,
                uniform_target,
                reserve,
                maximum_systems,
            )
            if child is None:
                child_specs = []
                break
            child_specs.append(child)
        if len(child_specs) != 2:
            nodes[work.node_id] = unresolved_node(
                work, "QUANTIZED_SPLIT_HAS_EMPTY_REAL_CHILD", active, total
            )
            residual += 1
            continue
        for child in child_specs:
            push(child)
        nodes[work.node_id] = {
            "node_id": work.node_id,
            "state": "SPLIT",
            "left": child_specs[0].node_id,
            "right": child_specs[1].node_id,
            "split": {
                "coefficients": list(coefficients),
                "threshold": str(threshold),
            },
            "count": box_upper_record(active, work.lower, work.upper, total),
            "diagnostic": {
                **work.diagnostic,
                **split_diagnostic,
                "split_kind": split_kind,
                "split_policy": SPLIT_POLICY,
                "left_count_upper": str(child_specs[0].count_upper),
                "right_count_upper": str(child_specs[1].count_upper),
                "count_note": (
                    "child coordinate-box counts are independently verified uppers; "
                    "arbitrary affine children do not claim exact conservation"
                ),
            },
        }
        split_counts[split_kind] += 1
        allocated += 2
    ordered = [nodes[key] for key in sorted(nodes, key=lambda value: (len(value), value))]
    shard = {
        "schema": coordinate.SHARD_SCHEMA,
        "manifest_sha256": manifest_digest,
        "support_mask": census_row["mask"],
        "active_classes": list(active),
        "root_count": str(root_count),
        "root_node": "r",
        "nodes": ordered,
        "aggregation_requested": "collapsed",
        "state": "INCOMPLETE",
        "discovery": {
            "producer": Path(__file__).name,
            "producer_mode": SPLIT_POLICY,
            "supplementary_sources": supplementary_sources,
            "combined_witness_count": len(bank.references),
            "dominance_quantization_scale": str(dominance_scale),
            "maximum_vertex_systems": maximum_systems,
            "node_counts": {
                "split_dominance": split_counts["dominance"],
                "split_coordinate": split_counts["coordinate"],
                "unresolved": passing + residual,
                "diagnostic_pass_awaiting_outward": passing,
                "diagnostic_residual": residual,
            },
            "completion_blocker": (
                "supplementary sources are discovery-only and every terminal remains "
                "UNRESOLVED pending manifest binding and independent outward replay"
            ),
        },
    }
    summary = {
        "support_mask": census_row["mask"],
        "dimension": dimension,
        "root_count": str(root_count),
        "nodes": len(ordered),
        "dominance_splits": split_counts["dominance"],
        "coordinate_splits": split_counts["coordinate"],
        "diagnostic_pass_awaiting_outward": passing,
        "diagnostic_residual": residual,
        "state": "INCOMPLETE",
    }
    return shard, summary


def synthetic_self_test() -> None:
    # Two incompatible leaders have the same pivot coordinate a2, so a split
    # on that stalled corner coordinate cannot separate them.  The quantized
    # dominance difference a0-a1 does separate them.
    vertices = (
        (Fraction(1), Fraction(2), Fraction(1)) + (Fraction(0),) * 6,
        (Fraction(2), Fraction(1), Fraction(1)) + (Fraction(0),) * 6,
        (Fraction(1), Fraction(1), Fraction(2)) + (Fraction(0),) * 6,
    )
    references = (
        {"source_id": "left", "source_sha256": "a" * 64, "row": 0, "row_sha256": "b" * 64},
        {"source_id": "right", "source_sha256": "c" * 64, "row": 0, "row_sha256": "d" * 64},
    )
    charges = np.zeros((2, CLASSES))
    charges[0, 0] = 1.0
    charges[0, 1] = -1.0
    charges[1, 0] = -1.0
    charges[1, 1] = 1.0
    bank = coordinate.AtlasBank(np.zeros(2), charges, references)
    leaders = np.asarray([1, 0, 0], dtype=np.int64)
    cut = choose_dominance_split(bank, vertices, leaders, (0, 1, 2), 4, 1024)
    if cut is None:
        raise AssertionError("synthetic dominance split was not found")
    coefficients, threshold, diagnostic = cut
    sides = [
        sum(Fraction(value) * point for value, point in zip(coefficients, vertex)) <= threshold
        for vertex in vertices[:2]
    ]
    if sides[0] == sides[1] or math.gcd(*map(abs, coefficients)) != 1:
        raise AssertionError("synthetic primitive dominance split did not separate leaders")
    if coefficients[2] != 0 or any(coefficients[index] for index in range(3, CLASSES)):
        raise AssertionError("synthetic support-relative pivot gauge failed")
    if next(coefficients[index] for index in (0, 1, 2) if coefficients[index]) < 0:
        raise AssertionError("synthetic canonical sign rule failed")
    if vertices[0][2] != vertices[1][2]:
        raise AssertionError("synthetic coordinate-stall premise changed")
    print(f"synthetic_dominance_coefficients={list(coefficients)}")
    print(f"synthetic_dominance_threshold={threshold}")
    print(f"synthetic_leader_sides={sides}")
    print(f"synthetic_soundness_scope={diagnostic['soundness_scope']}")
    print("status=SYNTHETIC_DOMINANCE_SPLIT_SMOKE_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=coordinate.DEFAULT_MANIFEST)
    parser.add_argument("--base-atlas", type=Path, default=coordinate.DEFAULT_ATLAS)
    parser.add_argument("--supplementary-atlas", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dimensions", type=parse_int_list, default=(8, 7, 6))
    parser.add_argument("--masks", type=parse_mask_list)
    parser.add_argument("--shard-start", type=int, default=0)
    parser.add_argument("--max-shards", type=int, default=0)
    parser.add_argument("--max-nodes-per-shard", type=int, default=63)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--hardening-reserve-bits", type=float, default=16.0)
    parser.add_argument("--dominance-quantization-scale", type=int, default=1 << 20)
    parser.add_argument("--max-vertex-systems", type=int, default=250000)
    parser.add_argument("--synthetic-self-test", action="store_true")
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if args.output_dir is None or not args.supplementary_atlas:
        parser.error("--output-dir and at least one --supplementary-atlas are required")
    if (
        args.shard_start < 0
        or args.max_shards < 0
        or args.max_nodes_per_shard < 1
        or args.max_nodes_per_shard % 2 == 0
        or args.max_depth < 0
        or args.dominance_quantization_scale <= 0
        or args.dominance_quantization_scale > (1 << 52)
        or args.dominance_quantization_scale
        & (args.dominance_quantization_scale - 1)
        or args.max_vertex_systems <= 0
        or args.hardening_reserve_bits < 0
    ):
        parser.error("invalid batch, node, quantization, or vertex-system budget")
    dimensions = set(args.dimensions)
    if not dimensions or not dimensions <= {6, 7, 8}:
        parser.error("this producer accepts only dimensions 6, 7, and 8")
    try:
        manifest, atlas, source, bindings, manifest_digest, base_digest = (
            coordinate.load_frozen_inputs(args.manifest, args.base_atlas)
        )
        base_bank = coordinate.build_atlas_bank(atlas, source, bindings)
        supplements = [
            load_supplementary_atlas(path, manifest_digest, base_digest)
            for path in args.supplementary_atlas
        ]
        bank = combined_bank(base_bank, supplements)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"g=8 dominance producer: {error}") from error
    supplementary_sources = [row[3] for row in supplements]
    parameters = manifest["parameters"]
    total = int(parameters["M"])
    global_count = int(manifest["census"]["feasible_profile_count"])
    numerator, denominator = map(int, parameters["probability_target_log2"].split("/"))
    probability_target = numerator / denominator
    masks = set(args.masks or ())
    census_rows = [
        row for row in manifest["census"]["feasible_masks"]
        if int(row["dimension"]) in dimensions
        and (not masks or int(row["mask"], 0) in masks)
    ]
    census_rows.sort(key=lambda row: (-int(row["dimension"]), int(row["mask"], 0)))
    census_rows = census_rows[args.shard_start:]
    if args.max_shards:
        census_rows = census_rows[:args.max_shards]
    if not census_rows:
        parser.error("batch selection contains no support shards")
    records = []
    for ordinal, census_row in enumerate(census_rows, 1):
        shard, summary = produce_shard(
            manifest_digest, census_row, bank, supplementary_sources, total,
            global_count, probability_target, args.max_nodes_per_shard,
            args.max_depth, args.hardening_reserve_bits,
            args.dominance_quantization_scale, args.max_vertex_systems,
        )
        filename = f"support_{int(census_row['mask'], 0):03x}.json"
        digest = coordinate.write_canonical(args.output_dir / filename, shard)
        records.append({**summary, "path": filename, "sha256": digest})
        print(
            f"shard={ordinal}/{len(census_rows)} mask={summary['support_mask']} "
            f"nodes={summary['nodes']} dominance={summary['dominance_splits']} "
            f"coordinate={summary['coordinate_splits']} "
            f"pass={summary['diagnostic_pass_awaiting_outward']}",
            flush=True,
        )
    batch = {
        "schema": BATCH_SCHEMA,
        "status": "INCOMPLETE_BINARY64_DOMINANCE_DISCOVERY_BATCH",
        "manifest_sha256": manifest_digest,
        "base_atlas_sha256": base_digest,
        "supplementary_sources": supplementary_sources,
        "configuration": {
            "dimensions": sorted(dimensions, reverse=True),
            "requested_masks": [] if args.masks is None else [f"0x{value:03x}" for value in args.masks],
            "shard_start": args.shard_start,
            "max_shards": args.max_shards,
            "max_nodes_per_shard": args.max_nodes_per_shard,
            "max_depth": args.max_depth,
            "hardening_reserve_bits": args.hardening_reserve_bits,
            "dominance_quantization_scale": str(args.dominance_quantization_scale),
            "max_vertex_systems": args.max_vertex_systems,
            "split_policy": SPLIT_POLICY,
        },
        "shards": records,
        "all_shards_incomplete": True,
        "completion_blocker": (
            "supplementary witness sources are not in the frozen proof manifest and "
            "no independent outward replay has run"
        ),
    }
    digest = coordinate.write_canonical(args.output_dir / "batch_index.json", batch)
    print(f"batch_index={args.output_dir / 'batch_index.json'}")
    print(f"batch_sha256={digest}")
    print("status=INCOMPLETE_BINARY64_DOMINANCE_DISCOVERY_BATCH")


if __name__ == "__main__":
    main()
