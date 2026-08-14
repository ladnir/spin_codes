#!/usr/bin/env python3
"""Produce an unresolved exact-BSP diagnostic for g=8 reflection blocks.

The five reflection blocks are {0,8}, {1,6}, {2,7}, {3,5}, and {4}.  The
tree peels the exact majority-mass cap of each block.  Every two-class cap then
receives one coordinate split at the exact midpoint.  All split ownership and
vertices are exact.  Frozen affine rows are evaluated only with binary64, so
every terminal remains UNRESOLVED and the output makes no proof claim.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance
from certify_packet_group_g8_support_shard import (
    enumerate_vertices,
    root_constraints,
    split_constraint,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUPPLEMENTARY_ATLAS = ROOT / "out" / "g8_highdim_supplementary_final32.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_reflection_block_caps.json"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_SUPPLEMENTARY_SHA256 = (
    "d50af6e55d5bc1d0aa0a9de65025684fdfe9989cce0014c808a059c2ca356d6c"
)
SCHEMA = "permute-conv.packet-group-g8-reflection-block-cap-diagnostic.v1"
BLOCKS = ((0, 8), (1, 6), (2, 7), (3, 5), (4,))
ACTIVE = tuple(range(9))
CLASSES = 9


@dataclass(frozen=True)
class Cut:
    coefficients: tuple[int, ...]
    threshold: int
    swapped: bool
    raw_coefficients: tuple[int, ...]
    raw_threshold: int
    label: str


@dataclass
class GeometryNode:
    state: str
    region: str
    cut: Cut | None = None
    left: "GeometryNode | None" = None
    right: "GeometryNode | None" = None


def canonical_cut(
    coefficients: Iterable[int], threshold: int, total: int
) -> tuple[tuple[int, ...], int, bool]:
    """Apply the verifier's support-relative primitive integer convention."""

    raw = tuple(int(value) for value in coefficients)
    if len(raw) != CLASSES:
        raise ValueError("a full-support cut needs nine coefficients")
    pivot_value = raw[-1]
    reduced = tuple(value - pivot_value for value in raw)
    canonical_threshold = int(threshold) - pivot_value * total
    if not any(reduced):
        raise ValueError("cut is constant on full support")
    divisor = math.gcd(*(abs(value) for value in reduced))
    reduced = tuple(value // divisor for value in reduced)
    canonical_threshold //= divisor
    swapped = False
    if next(value for value in reduced if value) < 0:
        reduced = tuple(-value for value in reduced)
        canonical_threshold = -canonical_threshold - 1
        swapped = True
    if max(abs(value).bit_length() for value in reduced) > 256:
        raise ValueError("canonical coefficient exceeds 256 bits")
    if abs(canonical_threshold).bit_length() > 256:
        raise ValueError("canonical threshold exceeds 256 bits")
    return reduced, canonical_threshold, swapped


def make_cut(
    classes: tuple[int, ...], threshold: int, total: int, label: str
) -> Cut:
    raw = tuple(1 if index in classes else 0 for index in range(CLASSES))
    coefficients, canonical_threshold, swapped = canonical_cut(raw, threshold, total)
    return Cut(
        coefficients=coefficients,
        threshold=canonical_threshold,
        swapped=swapped,
        raw_coefficients=raw,
        raw_threshold=threshold,
        label=label,
    )


def unresolved(region: str) -> GeometryNode:
    return GeometryNode(state="UNRESOLVED", region=region)


def split_node(
    cut: Cut,
    raw_low: GeometryNode,
    raw_high: GeometryNode,
    region: str,
) -> GeometryNode:
    left, right = (raw_high, raw_low) if cut.swapped else (raw_low, raw_high)
    return GeometryNode(
        state="SPLIT",
        region=region,
        cut=cut,
        left=left,
        right=right,
    )


def build_geometry(total: int) -> GeometryNode:
    """Build five disjoint majority caps and one midpoint split per pair cap."""

    majority_floor = total // 2

    def cap_leaf(block_index: int) -> GeometryNode:
        block = BLOCKS[block_index]
        cap_name = f"cap_g{block_index}"
        if len(block) == 1:
            return unresolved(cap_name)
        coordinate_cut = make_cut(
            (block[0],),
            majority_floor,
            total,
            f"{cap_name}:a{block[0]}<={majority_floor}",
        )
        return split_node(
            coordinate_cut,
            unresolved(f"{cap_name}_a{block[0]}_at_most_half"),
            unresolved(f"{cap_name}_a{block[0]}_above_half"),
            cap_name,
        )

    continuation = unresolved("central_all_block_masses_at_most_half")
    for block_index in reversed(range(len(BLOCKS))):
        block_cut = make_cut(
            BLOCKS[block_index],
            majority_floor,
            total,
            f"g{block_index}_mass<={majority_floor}",
        )
        continuation = split_node(
            block_cut,
            continuation,
            cap_leaf(block_index),
            f"peel_g{block_index}",
        )
    return continuation


def raw_value(cut: Cut, profile: tuple[int, ...]) -> int:
    return sum(left * right for left, right in zip(cut.raw_coefficients, profile))


def canonical_value(cut: Cut, profile: tuple[int, ...]) -> int:
    return sum(left * right for left, right in zip(cut.coefficients, profile))


def route(root: GeometryNode, profile: tuple[int, ...]) -> str:
    node = root
    while node.state == "SPLIT":
        assert node.cut is not None and node.left is not None and node.right is not None
        node = node.left if canonical_value(node.cut, profile) <= node.cut.threshold else node.right
    return node.region


def semantic_region(profile: tuple[int, ...], total: int) -> str:
    majority_floor = total // 2
    for block_index, block in enumerate(BLOCKS):
        if sum(profile[index] for index in block) > majority_floor:
            if len(block) == 1:
                return f"cap_g{block_index}"
            side = "at_most_half" if profile[block[0]] <= majority_floor else "above_half"
            return f"cap_g{block_index}_a{block[0]}_{side}"
    return "central_all_block_masses_at_most_half"


def positive_compositions(total: int, slots: int) -> Iterable[tuple[int, ...]]:
    def visit(remaining: int, positions: int, prefix: tuple[int, ...]):
        if positions == 1:
            yield prefix + (remaining,)
            return
        for value in range(1, remaining - positions + 2):
            yield from visit(remaining - value, positions - 1, prefix + (value,))

    if total >= slots:
        yield from visit(total, slots, ())


def run_self_test() -> None:
    # M=18 is the smallest convenient mass at which all five strict-majority
    # caps and both sides of every pair-coordinate split contain profiles.
    total = 18
    geometry = build_geometry(total)
    counts: dict[str, int] = {}
    seen = 0
    for profile in positive_compositions(total, CLASSES):
        actual = route(geometry, profile)
        expected = semantic_region(profile, total)
        if actual != expected:
            raise AssertionError(f"canonical ownership mismatch: {profile}: {actual} != {expected}")
        counts[actual] = counts.get(actual, 0) + 1
        seen += 1
    if seen != math.comb(total - 1, CLASSES - 1):
        raise AssertionError("synthetic full-support census mismatch")
    expected_regions = {
        "central_all_block_masses_at_most_half",
        "cap_g4",
        *(
            f"cap_g{i}_a{block[0]}_{side}"
            for i, block in enumerate(BLOCKS[:4])
            for side in ("at_most_half", "above_half")
        ),
    }
    if set(counts) != expected_regions or any(value <= 0 for value in counts.values()):
        raise AssertionError(f"synthetic geometry has an empty terminal: {counts}")

    # For every cut, the canonical side must match the raw side after applying
    # the declared exchange.  This also exercises the pivot-containing G0 cut.
    stack = [geometry]
    while stack:
        node = stack.pop()
        if node.state != "SPLIT":
            continue
        assert node.cut is not None and node.left is not None and node.right is not None
        for profile in positive_compositions(total, CLASSES):
            raw_low = raw_value(node.cut, profile) <= node.cut.raw_threshold
            canonical_low = canonical_value(node.cut, profile) <= node.cut.threshold
            if canonical_low != (not raw_low if node.cut.swapped else raw_low):
                raise AssertionError("canonical child exchange is inconsistent")
        stack.extend((node.left, node.right))


def geometry_statistics(root: GeometryNode) -> dict[str, int]:
    counts = {"nodes": 0, "splits": 0, "leaves": 0, "maximum_depth": 0}

    def visit(node: GeometryNode, depth: int) -> None:
        counts["nodes"] += 1
        counts["maximum_depth"] = max(counts["maximum_depth"], depth)
        if node.state == "SPLIT":
            counts["splits"] += 1
            assert node.left is not None and node.right is not None
            visit(node.left, depth + 1)
            visit(node.right, depth + 1)
        else:
            counts["leaves"] += 1

    visit(root, 0)
    return counts


def evaluate_geometry(
    root: GeometryNode,
    bank: coordinate.AtlasBank,
    total: int,
    uniform_target: float,
    hardening_reserve: float,
    root_count: int,
    maximum_systems: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    nodes: list[dict[str, Any]] = []
    leaf_diagnostics: list[dict[str, Any]] = []
    root_constraints_value = root_constraints(ACTIVE, total, 21)

    def visit(
        node: GeometryNode,
        identifier: str,
        constraints: tuple,
    ) -> None:
        systems = math.comb(len(constraints), len(ACTIVE) - 1)
        if systems > maximum_systems:
            raise ValueError(
                f"node {identifier} needs {systems} vertex systems; cap is {maximum_systems}"
            )
        vertices = enumerate_vertices(ACTIVE, constraints, maximum_systems)
        if not vertices:
            raise ValueError(f"structured geometry produced empty node {identifier}")
        common = {
            "node_id": identifier,
            "state": node.state,
            "region": node.region,
            "count": {
                "kind": "upper",
                "method": "root-census-upper-v1",
                "value": str(root_count),
            },
            "exact_vertex_count": len(vertices),
            "vertex_reconstruction_systems": systems,
        }
        if node.state == "UNRESOLVED":
            diagnostic, _scores, _leaders = dominance.score_vertices(
                bank, vertices, total, uniform_target, hardening_reserve
            )
            common.update(
                {
                    "reason": "AWAITING_INDEPENDENT_OUTWARD_REPLAY",
                    "diagnostic": diagnostic,
                }
            )
            nodes.append(common)
            leaf_diagnostics.append({"node_id": identifier, "region": node.region, **diagnostic})
            return

        assert node.cut is not None and node.left is not None and node.right is not None
        common.update(
            {
                "split": {
                    "coefficients": list(node.cut.coefficients),
                    "threshold": str(node.cut.threshold),
                },
                "left": identifier + "0",
                "right": identifier + "1",
                "selection": {
                    "kind": "fixed-reflection-block-majority-v1",
                    "label": node.cut.label,
                    "raw_coefficients": list(node.cut.raw_coefficients),
                    "raw_threshold": str(node.cut.raw_threshold),
                    "canonical_child_swap": node.cut.swapped,
                },
            }
        )
        nodes.append(common)
        visit(
            node.left,
            identifier + "0",
            constraints
            + (split_constraint(node.cut.coefficients, node.cut.threshold, ACTIVE, True),),
        )
        visit(
            node.right,
            identifier + "1",
            constraints
            + (split_constraint(node.cut.coefficients, node.cut.threshold, ACTIVE, False),),
        )

    visit(root, "r", root_constraints_value)
    maximum = max(float(row["candidate_upper_log2"]) for row in leaf_diagnostics)
    collapsed = maximum + math.log2(root_count)
    return nodes, {
        "leaf_diagnostics": leaf_diagnostics,
        "maximum_leaf_upper_log2": maximum,
        "collapsed_support_contribution_log2": collapsed,
        "collapsed_probability_target_gap_bits": collapsed - (-40.0),
        "aggregation": "collapsed-root-count-v1",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=coordinate.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=coordinate.DEFAULT_ATLAS)
    parser.add_argument(
        "--supplementary-atlas", type=Path, default=DEFAULT_SUPPLEMENTARY_ATLAS
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--hardening-reserve-bits", type=float, default=8.0)
    parser.add_argument("--maximum-vertex-systems", type=int, default=20_000)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        print("reflection-block canonical ownership tests passed")
        return
    if args.maximum_vertex_systems < 1:
        raise ValueError("maximum-vertex-systems must be positive")
    if not math.isfinite(args.hardening_reserve_bits) or args.hardening_reserve_bits < 0:
        raise ValueError("hardening reserve must be finite and nonnegative")

    manifest, atlas, source, bindings, manifest_digest, atlas_digest = coordinate.load_frozen_inputs(
        args.manifest, args.atlas
    )
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError("reflection-block producer manifest pin changed")
    supplementary_digest = coordinate.sha256_path(args.supplementary_atlas)
    if supplementary_digest != EXPECTED_SUPPLEMENTARY_SHA256:
        raise ValueError(
            f"supplementary atlas changed: {supplementary_digest}; "
            f"expected {EXPECTED_SUPPLEMENTARY_SHA256}"
        )
    supplementary = dominance.load_supplementary_atlas(
        args.supplementary_atlas, manifest_digest, atlas_digest
    )
    bank = dominance.combined_bank(
        coordinate.build_atlas_bank(atlas, source, bindings), [supplementary]
    )

    total = int(manifest["parameters"]["M"])
    if total != 1 << 18:
        raise ValueError("the exact verifier geometry is frozen to M=2^18")
    root_census = next(
        row for row in manifest["census"]["feasible_masks"] if row["mask"] == "0x1ff"
    )
    root_count = int(root_census["count"])
    global_count = int(manifest["census"]["feasible_profile_count"])
    probability_target = float(
        manifest["parameters"]["probability_target_log2"].split("/")[0]
    )
    uniform_target = probability_target - math.log2(global_count)

    geometry = build_geometry(total)
    statistics = geometry_statistics(geometry)
    nodes, evaluation = evaluate_geometry(
        geometry,
        bank,
        total,
        uniform_target,
        args.hardening_reserve_bits,
        root_count,
        args.maximum_vertex_systems,
    )
    derived_worst_case_systems = math.comb(
        len(root_constraints(ACTIVE, total, 21)) + statistics["maximum_depth"],
        len(ACTIVE) - 1,
    )
    artifact = {
        "schema": SCHEMA,
        "status": "INCOMPLETE_DIAGNOSTIC_AWAITING_OUTWARD_REPLAY",
        "scope": "full support 0x1ff; frozen-witness binary64 discovery, not a proof",
        "manifest_sha256": manifest_digest,
        "base_atlas_binding": {
            "path": str(args.atlas),
            "source_id": source["source_id"],
            "sha256": atlas_digest,
        },
        "supplementary_atlas_binding": supplementary[3],
        "support_mask": "0x1ff",
        "active_classes": list(ACTIVE),
        "root_count": str(root_count),
        "root_node": "r",
        "geometry": {
            "kind": "reflection-block-majority-caps-v1",
            "blocks": [list(block) for block in BLOCKS],
            "cap_rule": f"block_mass>={total // 2 + 1}",
            "within_pair_rule": f"first_coordinate<={total // 2}",
            "caps_are_pairwise_disjoint": True,
            "maximum_vertex_systems_per_node": args.maximum_vertex_systems,
            "derived_worst_case_systems": derived_worst_case_systems,
            **statistics,
        },
        "nodes": nodes,
        "evaluation": evaluation,
    }
    digest = coordinate.write_canonical(args.output, artifact)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": digest,
                "nodes": statistics["nodes"],
                "leaves": statistics["leaves"],
                "maximum_leaf_upper_log2": evaluation["maximum_leaf_upper_log2"],
                "collapsed_probability_target_gap_bits": evaluation[
                    "collapsed_probability_target_gap_bits"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
