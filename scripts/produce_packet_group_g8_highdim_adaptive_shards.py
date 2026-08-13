#!/usr/bin/env python3
"""Produce incomplete adaptive coordinate-slab shards for high-dimensional g=8 supports.

The producer binds itself to the frozen support manifest and affine witness
atlas.  It builds exact, disjoint integer BSP ownership trees.  Binary64
witness scores guide refinement, but they never create CERTIFIED_LEAF nodes.
Every terminal remains UNRESOLVED until an independent outward hardener
replays the proposed witness at the exact vertices.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import gammaln


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
DEFAULT_ATLAS = ROOT / "out" / "g8_support_seed_atlas.json"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_ATLAS_SHA256 = (
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
)
MANIFEST_SCHEMA = "packet-group-g8-support-manifest-v1"
ATLAS_SCHEMA = "permute-conv.packet-group-g8-support-seed-atlas.v1"
SHARD_SCHEMA = "packet-group-g8-support-shard-v1"
BATCH_SCHEMA = "permute-conv.packet-group-g8-highdim-adaptive-batch.v1"
GROUP_BITS = 8
CLASSES = 9
CLASS_MULTIPLICITIES = np.asarray(
    [math.comb(GROUP_BITS, index) for index in range(CLASSES)], dtype=np.float64
)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def write_canonical(path: Path, value: Any) -> str:
    payload = canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return sha256_bytes(payload)


def parse_int_list(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item) for item in text.split(",") if item)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError("integer list must not be empty")
    return values


def parse_mask_list(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item, 0) for item in text.split(",") if item)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "expected comma-separated hexadecimal or decimal masks"
        ) from error
    if not values or any(not 0 < value < (1 << CLASSES) for value in values):
        raise argparse.ArgumentTypeError("support mask lies outside nine classes")
    return values


def load_frozen_inputs(manifest_path: Path, atlas_path: Path):
    manifest_digest = sha256_path(manifest_path)
    atlas_digest = sha256_path(atlas_path)
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError(
            f"manifest digest changed: {manifest_digest}; expected {EXPECTED_MANIFEST_SHA256}"
        )
    if atlas_digest != EXPECTED_ATLAS_SHA256:
        raise ValueError(
            f"atlas digest changed: {atlas_digest}; expected {EXPECTED_ATLAS_SHA256}"
        )
    manifest = json.loads(manifest_path.read_bytes())
    atlas = json.loads(atlas_path.read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unexpected support manifest schema")
    if atlas.get("schema") != ATLAS_SCHEMA:
        raise ValueError("unexpected support atlas schema")
    if manifest["arithmetic"].get("canonical_json") != "sorted-compact-json-v1":
        raise ValueError("unsupported manifest canonical JSON profile")
    source = manifest["witness_sources"][0]
    if source["sha256"] != atlas_digest or source["schema"] != ATLAS_SCHEMA:
        raise ValueError("manifest does not bind the supplied atlas")
    if len(source["rows"]) != len(atlas["rows"]):
        raise ValueError("manifest/atlas row count mismatch")
    row_bindings = []
    for index, (binding, row) in enumerate(zip(source["rows"], atlas["rows"])):
        digest = sha256_bytes(canonical_bytes(row))
        if int(binding["row"]) != index or binding["row_sha256"] != digest:
            raise ValueError(f"atlas row binding failed at row {index}")
        if int(row["ordinal"]) != index:
            raise ValueError(f"atlas ordinal changed at row {index}")
        row_bindings.append(binding)
    return manifest, atlas, source, row_bindings, manifest_digest, atlas_digest


def exact_box_count(active: tuple[int, ...], lower: tuple[int, ...], upper: tuple[int, ...], total: int) -> int:
    """Count positive-support profiles in a coordinate box with fixed sum."""

    size = len(active)
    if size == 0:
        return 0
    local_lower = [lower[index] for index in active]
    local_upper = [upper[index] for index in active]
    if any(low > high for low, high in zip(local_lower, local_upper)):
        return 0
    residual = total - sum(local_lower)
    widths = [high - low for low, high in zip(local_lower, local_upper)]
    result = 0
    for mask in range(1 << size):
        shifted = residual
        parity = 0
        for index, width in enumerate(widths):
            if mask & (1 << index):
                shifted -= width + 1
                parity ^= 1
        term = 0 if shifted < 0 else math.comb(shifted + size - 1, size - 1)
        result += -term if parity else term
    return result


def exact_vertices(
    active: tuple[int, ...], lower: tuple[int, ...], upper: tuple[int, ...], total: int
) -> tuple[tuple[int, ...], ...]:
    """Enumerate exact vertices of a box intersected with one sum hyperplane."""

    vertices: set[tuple[int, ...]] = set()
    for free_class in active:
        fixed = tuple(index for index in active if index != free_class)
        for choices in itertools.product((0, 1), repeat=len(fixed)):
            point = [0] * CLASSES
            for index, take_upper in zip(fixed, choices):
                point[index] = upper[index] if take_upper else lower[index]
            point[free_class] = total - sum(point)
            if lower[free_class] <= point[free_class] <= upper[free_class]:
                vertices.add(tuple(point))
    return tuple(sorted(vertices))


def count_record(
    active: tuple[int, ...], lower: tuple[int, ...], upper: tuple[int, ...], total: int
) -> dict[str, Any]:
    value = exact_box_count(active, lower, upper, total)
    return {
        # The coordinate box is exact for this coordinate-only producer.  The
        # independent verifier deliberately treats it as an upper record so
        # the same leaf remains valid after later affine path refinements.
        "kind": "upper",
        "method": "coordinate-box-ie-v1",
        "value": str(value),
        "lower": [str(value) for value in lower],
        "upper": [str(value) for value in upper],
        "diagnostic_exact_for_coordinate_path": True,
        "diagnostic_total": str(total),
        "diagnostic_active_classes": list(active),
        "diagnostic_minimum_physical_weight": "21",
        "diagnostic_weight_cut_exclusions": "0",
    }


def diagnostic_normalization(vertices: np.ndarray, total: int) -> np.ndarray:
    return (
        gammaln(total + 1.0)
        - np.sum(gammaln(vertices + 1.0), axis=1)
        + vertices @ np.log(CLASS_MULTIPLICITIES)
    ) / math.log(2.0)


@dataclass(frozen=True)
class AtlasBank:
    constants: np.ndarray
    charges: np.ndarray
    references: tuple[dict[str, Any], ...]


def build_atlas_bank(
    atlas: dict[str, Any], source: dict[str, Any], bindings: list[dict[str, Any]]
) -> AtlasBank:
    constants = []
    charges = []
    references = []
    for index, (row, binding) in enumerate(zip(atlas["rows"], bindings)):
        affine = row["affine"]
        constant = float(affine["constant_log2"])
        charge = np.asarray(affine["charge_log2"], dtype=np.float64)
        fugacities = np.asarray(row["inner"]["fugacities"], dtype=np.float64)
        if charge.shape != (CLASSES,) or fugacities.shape != (CLASSES,):
            raise ValueError(f"malformed affine atlas row {index}")
        if not math.isfinite(constant) or not np.all(np.isfinite(charge)):
            raise ValueError(f"non-finite affine atlas row {index}")
        # The producer may use a row on another support only because every
        # frozen atlas fugacity is strictly positive.  A later hardener must
        # independently replay this support-eligibility check.
        if np.any(fugacities <= 0.0):
            raise ValueError(f"atlas row {index} contains a zero fugacity")
        constants.append(constant)
        charges.append(charge)
        references.append(
            {
                "source_id": source["source_id"],
                "source_sha256": source["sha256"],
                "row": index,
                "row_sha256": binding["row_sha256"],
            }
        )
    return AtlasBank(
        constants=np.asarray(constants, dtype=np.float64),
        charges=np.asarray(charges, dtype=np.float64),
        references=tuple(references),
    )


def score_box(
    bank: AtlasBank,
    active: tuple[int, ...],
    lower: tuple[int, ...],
    upper: tuple[int, ...],
    total: int,
    uniform_target: float,
    hardening_reserve: float,
) -> dict[str, Any]:
    vertices = exact_vertices(active, lower, upper, total)
    if not vertices:
        return {
            "candidate_status": "EMPTY_EXACT_BOX",
            "exact_vertex_count": 0,
            "candidate_upper_log2": None,
            "uniform_target_gap_bits": None,
        }
    matrix = np.asarray(vertices, dtype=np.float64)
    normalizations = diagnostic_normalization(matrix, total)
    scores = bank.constants[:, None] - bank.charges @ matrix.T - normalizations[None, :]
    maxima = np.max(scores, axis=1)
    witness_index = int(np.argmin(maxima))
    candidate_scores = scores[witness_index]
    worst_index = int(np.argmax(candidate_scores))
    upper_value = float(candidate_scores[worst_index])
    target_with_reserve = uniform_target - hardening_reserve
    return {
        "candidate_status": (
            "PASS_BINARY64_NEEDS_OUTWARD_REPLAY"
            if upper_value <= target_with_reserve
            else "FAIL_BINARY64_REQUIRES_SPLIT"
        ),
        "candidate_selector": {
            "kind": "witness",
            "witness": bank.references[witness_index],
        },
        "candidate_witness_row": witness_index,
        "candidate_upper_log2": upper_value,
        "uniform_profile_target_log2": uniform_target,
        "hardening_reserve_bits": hardening_reserve,
        "target_with_reserve_log2": target_with_reserve,
        "uniform_target_gap_bits": upper_value - uniform_target,
        "reserved_target_gap_bits": upper_value - target_with_reserve,
        "positive_residual_bits": max(0.0, upper_value - target_with_reserve),
        "worst_vertex": list(vertices[worst_index]),
        "exact_vertex_count": len(vertices),
        "vertex_sha256": sha256_bytes(canonical_bytes([list(row) for row in vertices])),
        "arithmetic": "binary64-discovery-only",
    }


def balanced_threshold(
    active: tuple[int, ...],
    lower: tuple[int, ...],
    upper: tuple[int, ...],
    total: int,
    coordinate: int,
    parent_count: int,
) -> int:
    other_lower = sum(lower[index] for index in active if index != coordinate)
    other_upper = sum(upper[index] for index in active if index != coordinate)
    effective_low = max(lower[coordinate], total - other_upper)
    effective_high = min(upper[coordinate], total - other_lower)
    if effective_low >= effective_high:
        raise ValueError("selected coordinate is fixed on the effective node domain")
    low = effective_low
    high = effective_high - 1
    target = parent_count // 2
    while low < high:
        middle = (low + high) // 2
        left_upper = list(upper)
        left_upper[coordinate] = middle
        count = exact_box_count(active, lower, tuple(left_upper), total)
        if count < target:
            low = middle + 1
        else:
            high = middle
    candidates = {low}
    if low > effective_low:
        candidates.add(low - 1)
    return min(
        candidates,
        key=lambda threshold: (
            abs(
                2
                * exact_box_count(
                    active,
                    lower,
                    tuple(
                        threshold if index == coordinate else value
                        for index, value in enumerate(upper)
                    ),
                    total,
                )
                - parent_count
            ),
            threshold,
        ),
    )


def choose_split(
    active: tuple[int, ...],
    lower: tuple[int, ...],
    upper: tuple[int, ...],
    total: int,
    parent_count: int,
    diagnostic: dict[str, Any],
) -> tuple[int, int] | None:
    worst = diagnostic.get("worst_vertex")
    effective_bounds = {
        index: (
            max(
                lower[index],
                total - sum(upper[other] for other in active if other != index),
            ),
            min(
                upper[index],
                total - sum(lower[other] for other in active if other != index),
            ),
        )
        for index in active
    }
    candidates = [
        index
        for index in active
        if effective_bounds[index][0] < effective_bounds[index][1]
    ]
    if not candidates:
        return None
    if worst is not None:
        candidates.sort(
            key=lambda index: (
                -(
                    (worst[index] - effective_bounds[index][0])
                    / max(1, effective_bounds[index][1] - effective_bounds[index][0])
                ),
                -(effective_bounds[index][1] - effective_bounds[index][0]),
                index,
            )
        )
    else:
        candidates.sort(
            key=lambda index: (
                -(effective_bounds[index][1] - effective_bounds[index][0]),
                index,
            )
        )
    coordinate = candidates[0]
    threshold = balanced_threshold(
        active, lower, upper, total, coordinate, parent_count
    )
    if not effective_bounds[coordinate][0] <= threshold < effective_bounds[coordinate][1]:
        return None
    return coordinate, threshold


@dataclass
class WorkNode:
    node_id: str
    depth: int
    lower: tuple[int, ...]
    upper: tuple[int, ...]
    count: int
    diagnostic: dict[str, Any]


def unresolved_node(work: WorkNode, reason: str, total: int, active: tuple[int, ...]) -> dict[str, Any]:
    diagnostic = dict(work.diagnostic)
    diagnostic["unresolved_reason"] = reason
    if work.count:
        diagnostic["collapsed_candidate_contribution_log2"] = (
            math.log2(work.count) + float(diagnostic["candidate_upper_log2"])
            if diagnostic.get("candidate_upper_log2") is not None
            else None
        )
    return {
        "node_id": work.node_id,
        "state": "UNRESOLVED",
        "count": count_record(active, work.lower, work.upper, total),
        "diagnostic": diagnostic,
    }


def produce_shard(
    manifest_digest: str,
    census_row: dict[str, Any],
    bank: AtlasBank,
    total: int,
    global_count: int,
    probability_target: float,
    max_nodes: int,
    max_depth: int,
    hardening_reserve: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    active = tuple(int(value) for value in census_row["active_classes"])
    dimension = int(census_row["dimension"])
    if dimension not in (6, 7, 8):
        raise ValueError("high-dimensional producer accepts dimensions 6, 7, and 8")
    if int(census_row["weight_cut_exclusions"]) != 0:
        raise ValueError("high-dimensional root unexpectedly needs weight exclusions")
    lower = tuple(1 if index in active else 0 for index in range(CLASSES))
    maximum = total - len(active) + 1
    upper = tuple(maximum if index in active else 0 for index in range(CLASSES))
    root_count = exact_box_count(active, lower, upper, total)
    if root_count != int(census_row["count"]):
        raise ValueError(f"root count mismatch for support {census_row['mask']}")
    uniform_target = probability_target - math.log2(global_count)
    root_diagnostic = score_box(
        bank, active, lower, upper, total, uniform_target, hardening_reserve
    )
    root = WorkNode("r", 0, lower, upper, root_count, root_diagnostic)
    nodes: dict[str, dict[str, Any]] = {}
    pending: list[tuple[float, int, WorkNode]] = []
    serial = 0

    def priority(work: WorkNode) -> float:
        value = work.diagnostic.get("candidate_upper_log2")
        return -(math.log2(work.count) + float(value)) if work.count and value is not None else math.inf

    heapq.heappush(pending, (priority(root), serial, root))
    serial += 1
    allocated_nodes = 1
    split_count = 0
    passing_count = 0
    residual_count = 0
    max_observed_depth = 0
    while pending:
        _priority, _serial, work = heapq.heappop(pending)
        max_observed_depth = max(max_observed_depth, work.depth)
        status = work.diagnostic["candidate_status"]
        if status == "PASS_BINARY64_NEEDS_OUTWARD_REPLAY":
            nodes[work.node_id] = unresolved_node(
                work, "AWAITING_INDEPENDENT_OUTWARD_REPLAY", total, active
            )
            passing_count += 1
            continue
        if work.depth >= max_depth or allocated_nodes + 2 > max_nodes:
            reason = "MAX_DEPTH_REACHED" if work.depth >= max_depth else "NODE_BUDGET_EXHAUSTED"
            nodes[work.node_id] = unresolved_node(work, reason, total, active)
            residual_count += 1
            continue
        split = choose_split(active, work.lower, work.upper, total, work.count, work.diagnostic)
        if split is None:
            nodes[work.node_id] = unresolved_node(
                work, "NO_NONTRIVIAL_COORDINATE_SPLIT", total, active
            )
            residual_count += 1
            continue
        coordinate, threshold = split
        left_upper = list(work.upper)
        left_upper[coordinate] = threshold
        right_lower = list(work.lower)
        right_lower[coordinate] = threshold + 1
        child_specs = (
            (f"{work.node_id}0", work.lower, tuple(left_upper)),
            (f"{work.node_id}1", tuple(right_lower), work.upper),
        )
        children = []
        for child_id, child_lower, child_upper in child_specs:
            child_count = exact_box_count(active, child_lower, child_upper, total)
            if child_count <= 0:
                raise AssertionError("count-balanced split created an empty child")
            child_diagnostic = score_box(
                bank,
                active,
                child_lower,
                child_upper,
                total,
                uniform_target,
                hardening_reserve,
            )
            child = WorkNode(
                child_id,
                work.depth + 1,
                child_lower,
                child_upper,
                child_count,
                child_diagnostic,
            )
            heapq.heappush(pending, (priority(child), serial, child))
            serial += 1
            children.append(child)
        if sum(child.count for child in children) != work.count:
            raise AssertionError("exact split counts do not conserve parent count")
        coefficients = [0] * CLASSES
        coefficients[coordinate] = 1
        nodes[work.node_id] = {
            "node_id": work.node_id,
            "state": "SPLIT",
            "left": children[0].node_id,
            "right": children[1].node_id,
            "split": {
                "coefficients": coefficients,
                "threshold": str(threshold),
            },
            "count": count_record(active, work.lower, work.upper, total),
            "diagnostic": {
                **work.diagnostic,
                "split_policy": "worst-corner-count-balanced-coordinate-v1",
                "split_coordinate": coordinate,
                "left_count": str(children[0].count),
                "right_count": str(children[1].count),
                "count_conservation": str(work.count),
            },
        }
        split_count += 1
        allocated_nodes += 2

    ordered_nodes = [nodes[key] for key in sorted(nodes, key=lambda key: (len(key), key))]
    unresolved = sum(node["state"] == "UNRESOLVED" for node in ordered_nodes)
    shard = {
        "schema": SHARD_SCHEMA,
        "manifest_sha256": manifest_digest,
        "support_mask": census_row["mask"],
        "active_classes": list(active),
        "root_count": str(root_count),
        "root_node": "r",
        "nodes": ordered_nodes,
        "aggregation_requested": "collapsed",
        "state": "INCOMPLETE",
        "discovery": {
            "producer": Path(__file__).name,
            "producer_mode": "binary64-adaptive-coordinate-slab-v1",
            "uniform_profile_target_log2": uniform_target,
            "hardening_reserve_bits": hardening_reserve,
            "max_nodes": max_nodes,
            "max_depth": max_depth,
            "node_counts": {
                "split": split_count,
                "unresolved": unresolved,
                "diagnostic_pass_awaiting_outward": passing_count,
                "diagnostic_residual": residual_count,
            },
            "max_observed_depth": max_observed_depth,
            "completion_blocker": (
                "every terminal is UNRESOLVED; independent outward replay must replace "
                "passing candidates with CERTIFIED_LEAF nodes"
            ),
        },
    }
    summary = {
        "support_mask": census_row["mask"],
        "dimension": dimension,
        "root_count": str(root_count),
        "nodes": len(ordered_nodes),
        "splits": split_count,
        "diagnostic_pass_awaiting_outward": passing_count,
        "diagnostic_residual": residual_count,
        "max_observed_depth": max_observed_depth,
        "state": "INCOMPLETE",
    }
    return shard, summary


def synthetic_self_test() -> None:
    active = (0, 1, 2, 3)
    total = 17
    lower = (1, 1, 1, 1, 0, 0, 0, 0, 0)
    upper = (14, 14, 14, 14, 0, 0, 0, 0, 0)
    root_count = exact_box_count(active, lower, upper, total)
    if root_count != math.comb(total - 1, len(active) - 1):
        raise AssertionError("synthetic root count failed")
    diagnostic = {"worst_vertex": [14, 1, 1, 1, 0, 0, 0, 0, 0]}
    coordinate, threshold = choose_split(
        active, lower, upper, total, root_count, diagnostic
    ) or (None, None)
    if coordinate != 0 or threshold is None:
        raise AssertionError("synthetic split selection failed")
    left_upper = list(upper)
    left_upper[coordinate] = threshold
    right_lower = list(lower)
    right_lower[coordinate] = threshold + 1
    left_count = exact_box_count(active, lower, tuple(left_upper), total)
    right_count = exact_box_count(active, tuple(right_lower), upper, total)
    if left_count <= 0 or right_count <= 0 or left_count + right_count != root_count:
        raise AssertionError("synthetic count conservation failed")
    vertices = exact_vertices(active, lower, upper, total)
    if len(vertices) != 4 or any(sum(vertex) != total for vertex in vertices):
        raise AssertionError("synthetic vertex enumeration failed")
    highdim_active = tuple(range(7))
    highdim_total = 17
    highdim_count = math.comb(highdim_total - 1, len(highdim_active) - 1)
    census_row = {
        "mask": "0x07f",
        "active_classes": list(highdim_active),
        "dimension": 6,
        "count": str(highdim_count),
        "weight_cut_exclusions": "0",
    }
    reference = {
        "source_id": "synthetic-atlas",
        "source_sha256": "0" * 64,
        "row": 0,
        "row_sha256": "1" * 64,
    }
    passing_bank = AtlasBank(
        constants=np.asarray([-1.0e6]),
        charges=np.zeros((1, CLASSES)),
        references=(reference,),
    )
    passing_shard, passing_summary = produce_shard(
        "2" * 64,
        census_row,
        passing_bank,
        highdim_total,
        highdim_count,
        -40.0,
        3,
        2,
        0.0,
    )
    if (
        passing_shard["state"] != "INCOMPLETE"
        or passing_summary["diagnostic_pass_awaiting_outward"] != 1
        or passing_shard["nodes"][0]["state"] != "UNRESOLVED"
    ):
        raise AssertionError("synthetic diagnostic-pass shard failed")
    failing_bank = AtlasBank(
        constants=np.asarray([1.0e6]),
        charges=np.zeros((1, CLASSES)),
        references=(reference,),
    )
    failing_shard, failing_summary = produce_shard(
        "2" * 64,
        census_row,
        failing_bank,
        highdim_total,
        highdim_count,
        -40.0,
        3,
        2,
        0.0,
    )
    if failing_summary["splits"] != 1 or len(failing_shard["nodes"]) != 3:
        raise AssertionError("synthetic node-budget shard failed")
    root = next(node for node in failing_shard["nodes"] if node["node_id"] == "r")
    if root["state"] != "SPLIT" or root["split"]["coefficients"] != [1] + [0] * 8:
        raise AssertionError("synthetic canonical split failed")
    print(f"synthetic_root_count={root_count}")
    print(f"synthetic_split=a{coordinate}<={threshold}")
    print(f"synthetic_child_counts={left_count}+{right_count}")
    print("synthetic_shards=diagnostic-pass,node-budget-residual")
    print("status=SYNTHETIC_STATIC_SMOKE_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dimensions", type=parse_int_list, default=(8, 7, 6))
    parser.add_argument("--masks", type=parse_mask_list)
    parser.add_argument("--shard-start", type=int, default=0)
    parser.add_argument("--max-shards", type=int, default=0)
    parser.add_argument("--max-nodes-per-shard", type=int, default=63)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--hardening-reserve-bits", type=float, default=16.0)
    parser.add_argument("--synthetic-self-test", action="store_true")
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if args.output_dir is None:
        parser.error("--output-dir is required unless --synthetic-self-test is used")
    if (
        args.shard_start < 0
        or args.max_shards < 0
        or args.max_nodes_per_shard < 1
        or args.max_nodes_per_shard % 2 == 0
        or args.max_depth < 0
        or not math.isfinite(args.hardening_reserve_bits)
        or args.hardening_reserve_bits < 0.0
    ):
        parser.error("invalid nonnegative batch or odd node budget")
    dimensions = set(args.dimensions)
    if not dimensions or not dimensions <= {6, 7, 8}:
        parser.error("this lane accepts only dimensions 6, 7, and 8")

    try:
        manifest, atlas, source, bindings, manifest_digest, atlas_digest = (
            load_frozen_inputs(args.manifest, args.atlas)
        )
        bank = build_atlas_bank(atlas, source, bindings)
    except ValueError as error:
        raise SystemExit(f"g=8 high-dimensional producer: {error}") from error
    parameters = manifest["parameters"]
    total = int(parameters["M"])
    global_count = int(manifest["census"]["feasible_profile_count"])
    target_numerator, target_denominator = map(
        int, parameters["probability_target_log2"].split("/")
    )
    probability_target = target_numerator / target_denominator
    selected_masks = set(args.masks or ())
    census_rows = [
        row
        for row in manifest["census"]["feasible_masks"]
        if int(row["dimension"]) in dimensions
        and (not selected_masks or int(row["mask"], 0) in selected_masks)
    ]
    if selected_masks - {int(row["mask"], 0) for row in census_rows}:
        parser.error("a requested mask is absent or outside the selected dimensions")
    census_rows.sort(key=lambda row: (-int(row["dimension"]), int(row["mask"], 0)))
    census_rows = census_rows[args.shard_start :]
    if args.max_shards:
        census_rows = census_rows[: args.max_shards]
    if not census_rows:
        parser.error("batch selection contains no support shards")

    shard_records = []
    for ordinal, census_row in enumerate(census_rows):
        shard, summary = produce_shard(
            manifest_digest,
            census_row,
            bank,
            total,
            global_count,
            probability_target,
            args.max_nodes_per_shard,
            args.max_depth,
            args.hardening_reserve_bits,
        )
        filename = f"support_{int(census_row['mask'], 0):03x}.json"
        path = args.output_dir / filename
        digest = write_canonical(path, shard)
        shard_records.append({**summary, "path": filename, "sha256": digest})
        print(
            f"shard={ordinal + 1}/{len(census_rows)} mask={summary['support_mask']} "
            f"dim={summary['dimension']} nodes={summary['nodes']} "
            f"pass={summary['diagnostic_pass_awaiting_outward']} "
            f"residual={summary['diagnostic_residual']}",
            flush=True,
        )
    batch = {
        "schema": BATCH_SCHEMA,
        "status": "INCOMPLETE_BINARY64_DISCOVERY_BATCH",
        "manifest_sha256": manifest_digest,
        "atlas_sha256": atlas_digest,
        "configuration": {
            "argv": sys.argv,
            "dimensions": sorted(dimensions, reverse=True),
            "requested_masks": (
                [] if args.masks is None else [f"0x{value:03x}" for value in args.masks]
            ),
            "shard_start": args.shard_start,
            "max_shards": args.max_shards,
            "max_nodes_per_shard": args.max_nodes_per_shard,
            "max_depth": args.max_depth,
            "hardening_reserve_bits": args.hardening_reserve_bits,
            "split_policy": "worst-corner-count-balanced-coordinate-v1",
        },
        "shards": shard_records,
        "all_shards_incomplete": True,
        "completion_blocker": "independent outward leaf hardening has not run",
    }
    index_path = args.output_dir / "batch_index.json"
    batch_digest = write_canonical(index_path, batch)
    print(f"batch_index={index_path}")
    print(f"batch_sha256={batch_digest}")
    print("status=INCOMPLETE_BINARY64_DISCOVERY_BATCH")


if __name__ == "__main__":
    main()
