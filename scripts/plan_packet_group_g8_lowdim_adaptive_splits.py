#!/usr/bin/env python3
"""Propose exact adaptive split trees for failed low-dimensional g=8 roots.

The producer uses binary64 witness scores only to choose integer cuts and
candidate selectors.  Every emitted shard is deliberately ``INCOMPLETE`` and
every terminal is ``UNRESOLVED``.  An outward replay or a later hardening pass
must replace terminals by ``CERTIFIED_LEAF`` before any completion claim.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import certify_packet_group_g8_support_shard as verifier
import plan_packet_group_g8_lowdim_root_cover as rootplan


PLAN_SCHEMA = "permute-conv.packet-group-g8-lowdim-root-cover-plan.v1"
SHARD_SCHEMA = "packet-group-g8-support-shard-v1"
REPORT_SCHEMA = "permute-conv.packet-group-g8-lowdim-adaptive-split-plan.v1"
MASS = 262144
CUTOFF = 21

Profile = tuple[Fraction, ...]
Constraint = tuple[tuple[Fraction, ...], Fraction]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def serialize_vertices(vertices: tuple[Profile, ...]) -> list[list[str]]:
    return [[fraction_text(value) for value in profile] for profile in vertices]


def load_source_plan(path: Path) -> tuple[dict[str, Any], dict[int, dict[str, Any]]]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema") != PLAN_SCHEMA:
        raise ValueError("adaptive split planner: unexpected source-plan schema")
    atlas = artifact.get("atlas", {})
    if atlas.get("sha256") != rootplan.ATLAS_SHA256:
        raise ValueError("adaptive split planner: source plan used another atlas")
    plans = artifact.get("plans")
    if not isinstance(plans, list):
        raise ValueError("adaptive split planner: source plan has no roots")
    by_mask = {int(row["support_mask"], 16): row for row in plans}
    if len(by_mask) != len(plans):
        raise ValueError("adaptive split planner: duplicate source-plan mask")
    return artifact, by_mask


def canonical_support_split(
    coefficients: list[int], threshold: int, support: tuple[int, ...]
) -> tuple[tuple[int, ...], int, bool] | None:
    """Apply the verifier's support-relative primitive-affine-gap form."""

    if len(coefficients) != rootplan.CLASSES or not support:
        raise ValueError("adaptive split planner: malformed support split")
    active = set(support)
    coefficients = [value if index in active else 0 for index, value in enumerate(coefficients)]
    pivot = max(support)
    pivot_value = coefficients[pivot]
    for index in support:
        coefficients[index] -= pivot_value
    threshold -= pivot_value * MASS
    divisor = 0
    for index in support:
        value = coefficients[index]
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        return None
    coefficients = [value // divisor for value in coefficients]
    threshold = threshold // divisor
    first = next(coefficients[index] for index in support if coefficients[index])
    swapped = False
    if first < 0:
        coefficients = [-value for value in coefficients]
        threshold = -threshold - 1
        swapped = True
    if coefficients[pivot] != 0 or any(
        coefficients[index] for index in range(rootplan.CLASSES) if index not in active
    ):
        raise RuntimeError("adaptive split planner: support gauge failed")
    if any(abs(value).bit_length() > 256 for value in coefficients) or abs(threshold).bit_length() > 256:
        return None
    return tuple(coefficients), threshold, swapped


def dominance_split(
    left: dict[str, Any], right: dict[str, Any], support: tuple[int, ...]
) -> tuple[tuple[int, ...], int, bool] | None:
    """Integerize the exact dyadic dominance plane of two discovery rows."""

    left_affine = left["affine"]
    right_affine = right["affine"]
    left_charges = [float(value) for value in left_affine["charge_log2"]]
    right_charges = [float(value) for value in right_affine["charge_log2"]]
    left_constant = float(left_affine["constant_log2"])
    right_constant = float(right_affine["constant_log2"])
    if not all(
        math.isfinite(value)
        for value in [left_constant, right_constant, *left_charges, *right_charges]
    ):
        return None
    coefficients = [
        Fraction.from_float(right_charge) - Fraction.from_float(left_charge)
        for left_charge, right_charge in zip(left_charges, right_charges)
    ]
    right_hand = Fraction.from_float(right_constant) - Fraction.from_float(left_constant)
    if not any(coefficients):
        return None
    denominator = right_hand.denominator
    for value in coefficients:
        denominator = math.lcm(denominator, value.denominator)
    integer_coefficients = [
        value.numerator * (denominator // value.denominator) for value in coefficients
    ]
    integer_right = right_hand.numerator * (denominator // right_hand.denominator)
    return canonical_support_split(integer_coefficients, integer_right, support)


def split_vertices(
    support: tuple[int, ...],
    constraints: tuple[Constraint, ...],
    split: tuple[tuple[int, ...], int, bool],
    maximum_systems: int,
) -> tuple[
    tuple[Constraint, ...],
    tuple[Constraint, ...],
    tuple[Profile, ...],
    tuple[Profile, ...],
] | None:
    coefficients, threshold, _swapped = split
    left_constraints = constraints + (
        verifier.split_constraint(coefficients, threshold, support, True),
    )
    right_constraints = constraints + (
        verifier.split_constraint(coefficients, threshold, support, False),
    )
    left_vertices = verifier.enumerate_vertices(support, left_constraints, maximum_systems)
    right_vertices = verifier.enumerate_vertices(support, right_constraints, maximum_systems)
    if not left_vertices or not right_vertices:
        return None
    parent_vertices = verifier.enumerate_vertices(support, constraints, maximum_systems)
    if set(left_vertices) == set(parent_vertices) or set(right_vertices) == set(parent_vertices):
        return None
    return left_constraints, right_constraints, left_vertices, right_vertices


def coordinate_split(
    support: tuple[int, ...],
    vertices: tuple[Profile, ...],
    values: np.ndarray,
) -> tuple[tuple[int, ...], int, bool] | None:
    spans = []
    for coordinate in support:
        low = ceil_fraction(min(profile[coordinate] for profile in vertices))
        high = max(value.numerator // value.denominator for value in (profile[coordinate] for profile in vertices))
        if low < high:
            spans.append((high - low, coordinate, low, high))
    if not spans:
        return None
    _width, coordinate, low, high = max(spans, key=lambda row: (row[0], -row[1]))
    threshold = (low + high) // 2

    # On a line, place the integer threshold near the equality point of the
    # two endpoint-winning witnesses when those witnesses cross.
    if len(support) == 2 and len(vertices) >= 2:
        ordered = sorted(range(len(vertices)), key=lambda index: vertices[index][coordinate])
        left_vertex, right_vertex = ordered[0], ordered[-1]
        left_winner = int(np.argmin(values[left_vertex]))
        right_winner = int(np.argmin(values[right_vertex]))
        if left_winner != right_winner:
            left_gap = float(values[left_vertex, left_winner] - values[left_vertex, right_winner])
            right_gap = float(values[right_vertex, left_winner] - values[right_vertex, right_winner])
            if left_gap * right_gap < 0.0 and right_gap != left_gap:
                x0 = float(vertices[left_vertex][coordinate])
                x1 = float(vertices[right_vertex][coordinate])
                crossing = x0 - left_gap * (x1 - x0) / (right_gap - left_gap)
                threshold = min(high - 1, max(low, math.floor(crossing)))
    coefficients = [0] * rootplan.CLASSES
    coefficients[coordinate] = 1
    return canonical_support_split(coefficients, threshold, support)


def select_cell_candidate(
    vertices: tuple[Profile, ...],
    candidates: list[dict[str, Any]],
    row_hashes: dict[int, str],
    denominator: int,
    reduced_cost_tolerance: float,
) -> tuple[dict[str, Any], np.ndarray, dict[str, Any], dict[str, Any]]:
    values = rootplan.witness_value_matrix(vertices, candidates)
    singleton = rootplan.best_singleton(values)
    mixture = rootplan.optimize_minimax_mixture(
        values,
        rootplan.UNIFORM_TARGET_LOG2,
        seed_per_vertex=1,
        reduced_cost_tolerance=reduced_cost_tolerance,
        weight_tolerance=1e-12,
        denominator=denominator,
    )
    singleton_selector = {
        "kind": "witness",
        "witness": rootplan.witness_reference(candidates[singleton["index"]], row_hashes),
    }
    fixed_mixture = rootplan.mixture_selector(mixture, candidates, row_hashes)
    if mixture["score"] <= singleton["score"]:
        selected = {
            "selector": fixed_mixture,
            "score": float(mixture["score"]),
            "vertex_values": list(map(float, mixture["vertex_values"])),
            "kind": fixed_mixture["kind"],
        }
    else:
        selected = {
            "selector": singleton_selector,
            "score": float(singleton["score"]),
            "vertex_values": list(map(float, singleton["vertex_values"])),
            "kind": "witness",
        }
    return selected, values, singleton, mixture


class TreeBuilder:
    def __init__(
        self,
        support: tuple[int, ...],
        candidates: list[dict[str, Any]],
        row_hashes: dict[int, str],
        max_depth: int,
        denominator: int,
        reduced_cost_tolerance: float,
        maximum_systems: int,
    ) -> None:
        self.support = support
        self.candidates = candidates
        self.row_hashes = row_hashes
        self.max_depth = max_depth
        self.denominator = denominator
        self.reduced_cost_tolerance = reduced_cost_tolerance
        self.maximum_systems = maximum_systems
        self.nodes: list[dict[str, Any]] = []
        self.terminal_scores: list[float] = []
        self.cut_counts = {"dominance": 0, "coordinate": 0}

    def reserve(self) -> tuple[str, int]:
        index = len(self.nodes)
        identifier = f"n{index:06d}"
        self.nodes.append({"node_id": identifier})
        return identifier, index

    def terminal(
        self,
        identifier: str,
        index: int,
        depth: int,
        vertices: tuple[Profile, ...],
        selected: dict[str, Any],
        reason: str,
        singleton: dict[str, Any],
        mixture: dict[str, Any],
    ) -> str:
        self.terminal_scores.append(selected["score"])
        self.nodes[index] = {
            "node_id": identifier,
            "state": "UNRESOLVED",
            "diagnostic": {
                "reason": reason,
                "depth": depth,
                "candidate_selector": selected["selector"],
                "candidate_kind": selected["kind"],
                "candidate_root_upper_log2": selected["score"],
                "candidate_uniform_target_pass": (
                    selected["score"] <= rootplan.UNIFORM_TARGET_LOG2
                ),
                "vertices": serialize_vertices(vertices),
                "vertex_values_log2": selected["vertex_values"],
                "best_singleton_score": float(singleton["score"]),
                "minimax_mixture_score": float(mixture["score"]),
                "mixture_components": len(mixture["indices"]),
                "column_generation_rounds": int(mixture["column_generation_rounds"]),
                "full_column_scans": int(mixture["full_column_scans"]),
            },
        }
        return identifier

    def build(self, constraints: tuple[Constraint, ...], depth: int = 0) -> str:
        identifier, node_index = self.reserve()
        vertices = verifier.enumerate_vertices(
            self.support, constraints, self.maximum_systems
        )
        if not vertices:
            self.nodes[node_index] = {
                "node_id": identifier,
                "state": "EMPTY",
                "emptiness_method": "exact-rational-polytope-v1",
            }
            return identifier
        selected, values, singleton, mixture = select_cell_candidate(
            vertices,
            self.candidates,
            self.row_hashes,
            self.denominator,
            self.reduced_cost_tolerance,
        )
        if selected["score"] <= rootplan.UNIFORM_TARGET_LOG2:
            return self.terminal(
                identifier, node_index, depth, vertices, selected,
                "binary64_candidate_pass_requires_outward_replay", singleton, mixture,
            )
        if depth >= self.max_depth or len(self.support) == 1:
            return self.terminal(
                identifier, node_index, depth, vertices, selected,
                "depth_or_dimension_limit", singleton, mixture,
            )

        chosen = None
        cut_kind = None
        if len(self.support) >= 3:
            winners = np.argmin(values, axis=1)
            conflict_pairs = []
            for left_vertex, right_vertex in itertools.combinations(range(len(vertices)), 2):
                left_winner = int(winners[left_vertex])
                right_winner = int(winners[right_vertex])
                if left_winner == right_winner:
                    continue
                conflict_pairs.append(
                    (
                        -max(
                            selected["vertex_values"][left_vertex],
                            selected["vertex_values"][right_vertex],
                        ),
                        left_winner,
                        right_winner,
                    )
                )
            for _priority, left_winner, right_winner in sorted(conflict_pairs):
                split = dominance_split(
                    self.candidates[left_winner], self.candidates[right_winner], self.support
                )
                if split is None:
                    continue
                replay = split_vertices(
                    self.support, constraints, split, self.maximum_systems
                )
                if replay is not None:
                    chosen = (split, replay)
                    cut_kind = "dominance"
                    break
        if chosen is None:
            split = coordinate_split(self.support, vertices, values)
            if split is not None:
                replay = split_vertices(
                    self.support, constraints, split, self.maximum_systems
                )
                if replay is not None:
                    chosen = (split, replay)
                    cut_kind = "coordinate"
        if chosen is None:
            return self.terminal(
                identifier, node_index, depth, vertices, selected,
                "no_nontrivial_exact_split", singleton, mixture,
            )

        split, replay = chosen
        left_constraints, right_constraints, _left_vertices, _right_vertices = replay
        left = self.build(left_constraints, depth + 1)
        right = self.build(right_constraints, depth + 1)
        coefficients, threshold, child_swap = split
        self.nodes[node_index] = {
            "node_id": identifier,
            "state": "SPLIT",
            "left": left,
            "right": right,
            "split": {
                "coefficients": list(coefficients),
                "threshold": str(threshold),
            },
            "discovery": {
                "kind": cut_kind,
                "canonical_child_swap": child_swap,
                "parent_candidate_selector": selected["selector"],
                "parent_candidate_root_upper_log2": selected["score"],
            },
        }
        self.cut_counts[str(cut_kind)] += 1
        return identifier


def shard_for_root(
    source_row: dict[str, Any],
    atlas_rows: list[dict[str, Any]],
    row_hashes: dict[int, str],
    max_depth: int,
    denominator: int,
    reduced_cost_tolerance: float,
    maximum_systems: int,
    source_plan_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    mask = int(source_row["support_mask"], 16)
    support = verifier.support_of(mask)
    candidates = rootplan.candidate_rows(atlas_rows, support, "all")
    builder = TreeBuilder(
        support,
        candidates,
        row_hashes,
        max_depth,
        denominator,
        reduced_cost_tolerance,
        maximum_systems,
    )
    constraints = verifier.root_constraints(support, MASS, CUTOFF)
    root = builder.build(constraints)
    shard = {
        "schema": SHARD_SCHEMA,
        "manifest_sha256": rootplan.MANIFEST_SHA256,
        "support_mask": verifier.mask_text(mask),
        "active_classes": list(support),
        "root_count": str(source_row["root_count"]),
        "root_node": root,
        "nodes": builder.nodes,
        "aggregation_requested": "collapsed",
        "state": "INCOMPLETE",
        "discovery": {
            "source_plan_sha256": source_plan_sha256,
            "source_root_upper_log2": source_row["diagnostic_root_upper_log2"],
            "max_depth": max_depth,
            "mixture_denominator": denominator,
            "reduced_cost_tolerance": reduced_cost_tolerance,
            "warning": "terminal selectors are binary64 proposals pending outward replay",
        },
    }
    diagnostic = {
        "support_mask": verifier.mask_text(mask),
        "dimension": len(support) - 1,
        "nodes": len(builder.nodes),
        "terminals": sum(node.get("state") == "UNRESOLVED" for node in builder.nodes),
        "splits": sum(node.get("state") == "SPLIT" for node in builder.nodes),
        "cut_counts": builder.cut_counts,
        "worst_terminal_candidate_log2": max(builder.terminal_scores),
        "terminal_uniform_candidate_passes": sum(
            node.get("state") == "UNRESOLVED"
            and node.get("diagnostic", {}).get("candidate_uniform_target_pass", False)
            for node in builder.nodes
        ),
    }
    return shard, diagnostic


def run_self_test() -> None:
    clipped_support = (0, 1)
    constraints = verifier.root_constraints(clipped_support, MASS, CUTOFF)
    vertices = verifier.enumerate_vertices(clipped_support, constraints, 100)
    expected = (
        (Fraction(1), Fraction(MASS - 1), *([Fraction(0)] * 7)),
        (Fraction(MASS - CUTOFF), Fraction(CUTOFF), *([Fraction(0)] * 7)),
    )
    if vertices != expected:
        raise SystemExit("adaptive split self-test: clipped root vertices changed")
    values = np.asarray([[1.0, -3.0], [-3.0, 1.0]])
    split = coordinate_split(clipped_support, vertices, values)
    if split is None or verifier.validate_split(
        {"coefficients": list(split[0]), "threshold": str(split[1])},
        clipped_support,
    )[:2] != split[:2]:
        raise SystemExit("adaptive split self-test: line split is not canonical")
    replay = split_vertices(clipped_support, constraints, split, 100)
    if replay is None or not replay[2] or not replay[3]:
        raise SystemExit("adaptive split self-test: line split failed exact replay")

    left = {
        "affine": {"constant_log2": 0.0, "charge_log2": [0.0, 1.0, 0.0] + [0.0] * 6}
    }
    right = {
        "affine": {"constant_log2": 5.0, "charge_log2": [0.0, 0.0, 1.0] + [0.0] * 6}
    }
    dominance = dominance_split(left, right, (0, 1, 2))
    if dominance is None or math.gcd(*map(abs, dominance[0])) != 1:
        raise SystemExit("adaptive split self-test: dominance cut is not primitive")
    if next(value for value in dominance[0] if value) < 0:
        raise SystemExit("adaptive split self-test: dominance cut sign is not canonical")
    if dominance[0][2] != 0 or any(dominance[0][index] for index in range(3, 9)):
        raise SystemExit("adaptive split self-test: dominance cut gauge is not canonical")
    if any(abs(value).bit_length() > 256 for value in dominance[0]) or abs(dominance[1]).bit_length() > 256:
        raise SystemExit("adaptive split self-test: dominance cut exceeds bit cap")
    replayed_dominance = verifier.validate_split(
        {"coefficients": list(dominance[0]), "threshold": str(dominance[1])},
        (0, 1, 2),
    )
    if replayed_dominance != (dominance[0], dominance[1]):
        raise SystemExit("adaptive split self-test: verifier changed dominance cut")
    nonfinite = {
        "affine": {"constant_log2": math.inf, "charge_log2": [0.0] * 9}
    }
    if dominance_split(nonfinite, right, (0, 1, 2)) is not None:
        raise SystemExit("adaptive split self-test: nonfinite dominance row accepted")
    synthetic_shard = {
        "schema": SHARD_SCHEMA,
        "manifest_sha256": rootplan.MANIFEST_SHA256,
        "support_mask": "0x003",
        "active_classes": [0, 1],
        "root_count": "262123",
        "root_node": "n000000",
        "nodes": [
            {
                "node_id": "n000000",
                "state": "SPLIT",
                "left": "n000001",
                "right": "n000002",
                "split": {
                    "coefficients": list(split[0]),
                    "threshold": str(split[1]),
                },
            },
            {
                "node_id": "n000001",
                "state": "UNRESOLVED",
                "diagnostic": {"candidate_selector": {"kind": "witness"}},
            },
            {
                "node_id": "n000002",
                "state": "UNRESOLVED",
                "diagnostic": {"candidate_selector": {"kind": "witness"}},
            },
        ],
        "aggregation_requested": "collapsed",
        "state": "INCOMPLETE",
    }
    _leaves, _leaf_vertices, _node_vertices, counts = verifier.replay_tree(
        synthetic_shard, clipped_support, {}, 100
    )
    if (
        synthetic_shard["state"] != "INCOMPLETE"
        or counts["split"] != 1
        or counts["unresolved"] != 2
        or counts["certified_leaf"] != 0
    ):
        raise SystemExit("adaptive split self-test: false completion")
    print("physical_weight_clipped_root=PASS")
    print("dimension1_coordinate_interval_split=PASS")
    print("primitive_dominance_cut=PASS")
    print("incomplete_shard_gate=PASS")
    print("status=PASS_LOWDIM_ADAPTIVE_SPLIT_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-plan", type=Path)
    parser.add_argument("--manifest", type=Path, default=rootplan.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=rootplan.DEFAULT_ATLAS)
    parser.add_argument("--top-failures", type=int, default=4)
    parser.add_argument("--support-mask", type=lambda value: int(value, 0), action="append")
    parser.add_argument("--max-depth", type=int, default=2)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--reduced-cost-tolerance", type=float, default=1e-9)
    parser.add_argument("--maximum-vertex-systems", type=int, default=10000)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.source_plan is None or args.output_dir is None or args.report is None:
        parser.error("--source-plan, --output-dir, and --report are required")
    if (
        args.top_failures < 1
        or args.max_depth < 0
        or args.mixture_denominator < 1
        or args.reduced_cost_tolerance < 0
        or args.maximum_vertex_systems < 1
    ):
        parser.error("invalid positive adaptive-planner limit")

    _manifest, row_hashes, manifest_census = rootplan.load_manifest(args.manifest)
    _atlas, atlas_rows = rootplan.load_atlas(args.atlas)
    source, source_rows = load_source_plan(args.source_plan)
    if args.support_mask:
        missing = sorted(set(args.support_mask) - set(source_rows))
        if missing:
            raise ValueError(f"requested masks are absent from source plan: {missing}")
        selected = [source_rows[mask] for mask in args.support_mask]
    else:
        failed = [
            row for row in source_rows.values()
            if not bool(row.get("diagnostic_uniform_target_pass", False))
        ]
        selected = sorted(
            failed,
            key=lambda row: (
                -float(row["diagnostic_contribution_log2"]),
                int(row["support_mask"], 16),
            ),
        )[: args.top_failures]
    if not selected:
        raise ValueError("adaptive split planner selected no failed roots")
    for row in selected:
        mask = int(row["support_mask"], 16)
        if mask not in manifest_census or int(manifest_census[mask]["count"]) != int(row["root_count"]):
            raise ValueError(f"source-plan root disagrees with manifest at {row['support_mask']}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_digest = file_sha256(args.source_plan)
    diagnostics = []
    shard_files = []
    for row in selected:
        shard, diagnostic = shard_for_root(
            row,
            atlas_rows,
            row_hashes,
            args.max_depth,
            args.mixture_denominator,
            args.reduced_cost_tolerance,
            args.maximum_vertex_systems,
            source_digest,
        )
        path = args.output_dir / f"support_{shard['support_mask'][2:]}.json"
        path.write_text(json.dumps(shard, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        shard_files.append({"path": str(path), "sha256": file_sha256(path)})
        diagnostics.append(diagnostic)
        print(
            f"mask={diagnostic['support_mask']} nodes={diagnostic['nodes']} "
            f"splits={diagnostic['splits']} worst={diagnostic['worst_terminal_candidate_log2']:.9f}",
            flush=True,
        )
    report = {
        "schema": REPORT_SCHEMA,
        "status": "DIAGNOSTIC_INCOMPLETE_SHARDS_REQUIRING_OUTWARD_REPLAY",
        "manifest_sha256": rootplan.MANIFEST_SHA256,
        "source_plan": {
            "path": str(args.source_plan),
            "sha256": source_digest,
            "declared_manifest_sha256": source.get("manifest", {}).get("sha256"),
        },
        "configuration": {
            "top_failures": args.top_failures,
            "support_masks": args.support_mask,
            "max_depth": args.max_depth,
            "mixture_denominator": args.mixture_denominator,
            "reduced_cost_tolerance": args.reduced_cost_tolerance,
        },
        "shards": shard_files,
        "diagnostics": diagnostics,
        "warning": "all shards are INCOMPLETE and all candidate terminals require outward replay",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"report={args.report}")
    print("status=DIAGNOSTIC_INCOMPLETE_SHARDS_REQUIRING_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
