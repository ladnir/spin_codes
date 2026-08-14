#!/usr/bin/env python3
"""Produce an incomplete ordered-band diagnostic for full-support g=8.

Dominance scores choose one class order and exact target-profile masses choose
the split thresholds.  The resulting proof geometry is independent of those
binary64 choices: every cell is an integer slab

    lower <= sum(profile[j] for j in band) <= upper.

The script counts each slab exactly and enumerates its exact integral vertices.
It never emits a certified leaf.  An independent outward replay would still be
required before any witness selector could become a proof leaf.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUPPLEMENTARY_ATLAS = ROOT / "out" / "g8_highdim_supplementary_final32.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_ordered_band_diagnostic.json"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_SUPPLEMENTARY_SHA256 = (
    "d50af6e55d5bc1d0aa0a9de65025684fdfe9989cce0014c808a059c2ca356d6c"
)
SCHEMA = "permute-conv.packet-group-g8-full-support-ordered-band-diagnostic.v1"
SPLIT_POLICY = "dominance-ordered-single-prefix-target-mass-v1"
COUNT_METHOD = "two-band-positive-composition-scan-v1"
CLASSES = 9


@dataclass(frozen=True)
class BandChoice:
    order: tuple[int, ...]
    width: int
    witness_pair: tuple[int, int]
    alignment: int
    root_leaders: tuple[int, ...]

    @property
    def band(self) -> tuple[int, ...]:
        return self.order[: self.width]


def full_support_vertices(total: int, classes: int = CLASSES) -> tuple[tuple[int, ...], ...]:
    if total < classes:
        return ()
    result = []
    for dominant in range(classes):
        point = [1] * classes
        point[dominant] = total - classes + 1
        result.append(tuple(point))
    return tuple(result)


def slab_vertices(
    total: int,
    classes: int,
    band: tuple[int, ...],
    lower: int,
    upper: int,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate vertices of the positive simplex intersected with one slab."""

    band_set = frozenset(band)
    if not band or len(band_set) != len(band) or any(not 0 <= j < classes for j in band):
        raise ValueError("band must contain distinct in-range classes")
    complement = tuple(j for j in range(classes) if j not in band_set)
    if not complement:
        raise ValueError("band must be proper")
    q = len(band)
    r = len(complement)
    minimum = q
    maximum = total - r
    if not minimum <= lower <= upper <= maximum:
        raise ValueError("slab bounds lie outside the full-support simplex")

    vertices: set[tuple[int, ...]] = set()
    if lower == minimum:
        for dominant in complement:
            point = [1] * classes
            point[dominant] = total - classes + 1
            vertices.add(tuple(point))
    if upper == maximum:
        for dominant in band:
            point = [1] * classes
            point[dominant] = total - classes + 1
            vertices.add(tuple(point))

    boundary_values = set()
    if lower > minimum:
        boundary_values.add(lower)
    if upper < maximum:
        boundary_values.add(upper)
    for mass in sorted(boundary_values):
        left_value = mass - q + 1
        right_value = total - mass - r + 1
        for left in band:
            for right in complement:
                point = [1] * classes
                point[left] = left_value
                point[right] = right_value
                vertices.add(tuple(point))
    return tuple(sorted(vertices))


def slab_summand(total: int, q: int, mass: int, classes: int = CLASSES) -> int:
    r = classes - q
    if mass < q or total - mass < r:
        return 0
    return math.comb(mass - 1, q - 1) * math.comb(total - mass - 1, r - 1)


def exact_interval_counts(
    total: int,
    q: int,
    intervals: Iterable[tuple[int, int]],
    classes: int = CLASSES,
) -> dict[tuple[int, int], int]:
    """Count all requested slabs in one linear scan and O(intervals) memory."""

    requested = tuple(sorted(set(intervals)))
    r = classes - q
    minimum = q
    maximum = total - r
    for lower, upper in requested:
        if not minimum <= lower <= upper <= maximum:
            raise ValueError("requested interval lies outside feasible band masses")
    endpoints = {minimum - 1}
    for lower, upper in requested:
        endpoints.add(lower - 1)
        endpoints.add(upper)
    prefix: dict[int, int] = {minimum - 1: 0}
    running = 0
    for mass in range(minimum, maximum + 1):
        running += slab_summand(total, q, mass, classes)
        if mass in endpoints:
            prefix[mass] = running
    return {(lower, upper): prefix[upper] - prefix[lower - 1] for lower, upper in requested}


def score_matrix(bank: coordinate.AtlasBank, profiles: np.ndarray, total: int) -> np.ndarray:
    normalizations = coordinate.diagnostic_normalization(profiles, total)
    return bank.constants[:, None] - bank.charges @ profiles.T - normalizations[None, :]


def choose_dominance_order(bank: coordinate.AtlasBank, total: int) -> BandChoice:
    """Choose a prefix whose root corners separate two leading witnesses."""

    roots = np.asarray(full_support_vertices(total), dtype=np.float64)
    scores = score_matrix(bank, roots, total)
    leaders = np.argmin(scores, axis=0)
    distinct = tuple(sorted(int(value) for value in set(leaders.tolist())))
    if len(distinct) < 2:
        # A single root leader supplies no pairwise plane.  Charge order from
        # that witness is still deterministic, but the zero alignment records
        # that dominance did not expose an incompatibility.
        witness = distinct[0]
        order = tuple(sorted(range(CLASSES), key=lambda j: (bank.charges[witness, j], j)))
        return BandChoice(order, CLASSES // 2, (witness, witness), 0, tuple(int(x) for x in leaders))

    best_key: tuple[Any, ...] | None = None
    best: BandChoice | None = None
    for first_position, first in enumerate(distinct):
        for second in distinct[first_position + 1 :]:
            difference = scores[first] - scores[second]
            order = tuple(sorted(range(CLASSES), key=lambda j: (float(difference[j]), j)))
            for width in range(1, CLASSES):
                left = frozenset(order[:width])
                forward = sum(
                    1
                    for j in range(CLASSES)
                    if (j in left and leaders[j] == first) or (j not in left and leaders[j] == second)
                )
                reverse = sum(
                    1
                    for j in range(CLASSES)
                    if (j in left and leaders[j] == second) or (j not in left and leaders[j] == first)
                )
                if reverse > forward:
                    oriented_order = tuple(reversed(order))
                    oriented_pair = (second, first)
                    alignment = reverse
                    oriented_width = CLASSES - width
                else:
                    oriented_order = order
                    oriented_pair = (first, second)
                    alignment = forward
                    oriented_width = width
                balance = min(oriented_width, CLASSES - oriented_width)
                key = (alignment, balance, -oriented_pair[0], -oriented_pair[1], tuple(-j for j in oriented_order))
                if best_key is None or key > best_key:
                    best_key = key
                    best = BandChoice(
                        oriented_order,
                        oriented_width,
                        oriented_pair,
                        alignment,
                        tuple(int(x) for x in leaders),
                    )
    assert best is not None
    return best


def choose_target_boundaries(
    supplementary_artifact: dict[str, Any],
    bank: coordinate.AtlasBank,
    band: tuple[int, ...],
    total: int,
    maximum_boundaries: int,
) -> tuple[tuple[int, ...], dict[str, Any]]:
    """Choose exact target masses, prioritizing changes in the leading witness."""

    if maximum_boundaries < 0:
        raise ValueError("maximum_boundaries must be nonnegative")
    q = len(band)
    minimum = q
    maximum = total - (CLASSES - q)
    rows = supplementary_artifact.get("rows", [])
    anchors = []
    for row_index, row in enumerate(rows):
        profile = tuple(int(value) for value in row.get("profile", []))
        if len(profile) != CLASSES or sum(profile) != total or any(value < 1 for value in profile):
            raise ValueError(f"supplementary target profile is malformed at row {row_index}")
        anchors.append((sum(profile[j] for j in band), row_index, profile))
    if not anchors or maximum_boundaries == 0:
        return (), {"anchor_count": len(anchors), "leader_transition_masses": []}

    matrix = np.asarray([profile for _, _, profile in anchors], dtype=np.float64)
    leaders = np.argmin(score_matrix(bank, matrix, total), axis=0)
    annotated = sorted(
        (mass, row_index, int(leader))
        for (mass, row_index, _profile), leader in zip(anchors, leaders)
    )
    transitions = []
    for left, right in zip(annotated, annotated[1:]):
        if left[0] < right[0] and left[2] != right[2] and minimum <= left[0] < maximum:
            transitions.append(left[0])

    selected: list[int] = []
    for mass in sorted(set(transitions)):
        if len(selected) == maximum_boundaries:
            break
        selected.append(mass)
    candidate_masses = sorted(
        set(mass for mass, _row, _leader in annotated if minimum <= mass < maximum)
    )
    remaining = maximum_boundaries - len(selected)
    if remaining > 0 and candidate_masses:
        # Exact rank quantiles of target masses fill gaps left by leader
        # transitions.  No floating value is converted into a proof threshold.
        for rank in range(1, remaining + 1):
            index = (rank * len(candidate_masses)) // (remaining + 1)
            selected.append(candidate_masses[min(index, len(candidate_masses) - 1)])
    selected = sorted(set(selected))[:maximum_boundaries]
    return tuple(selected), {
        "anchor_count": len(annotated),
        "anchor_band_masses": [mass for mass, _row, _leader in annotated],
        "anchor_leading_witnesses": [leader for _mass, _row, leader in annotated],
        "leader_transition_masses": sorted(set(transitions)),
    }


def canonical_band_split(
    band: tuple[int, ...], threshold: int, total: int
) -> tuple[list[int], int, bool]:
    """Return the verifier's support-relative canonical cut and swap flag."""

    coefficients = [1 if j in band else 0 for j in range(CLASSES)]
    pivot_value = coefficients[-1]
    coefficients = [value - pivot_value for value in coefficients]
    coefficients[-1] = 0
    canonical_threshold = threshold - pivot_value * total
    first = next((value for value in coefficients if value), 0)
    swapped = False
    if first < 0:
        coefficients = [-value for value in coefficients]
        canonical_threshold = -canonical_threshold - 1
        swapped = True
    if not first:
        raise ValueError("band split is constant on full support")
    return coefficients, canonical_threshold, swapped


def make_intervals(minimum: int, maximum: int, boundaries: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    valid = tuple(sorted(set(boundaries)))
    if any(not minimum <= value < maximum for value in valid):
        raise ValueError("boundary lies outside feasible mass range")
    lower = minimum
    result = []
    for threshold in valid:
        result.append((lower, threshold))
        lower = threshold + 1
    result.append((lower, maximum))
    return tuple(result)


def score_slab(
    bank: coordinate.AtlasBank,
    vertices: tuple[tuple[int, ...], ...],
    total: int,
    uniform_target: float,
    hardening_reserve: float,
) -> dict[str, Any]:
    matrix = np.asarray(vertices, dtype=np.float64)
    scores = score_matrix(bank, matrix, total)
    maxima = np.max(scores, axis=1)
    witness = int(np.argmin(maxima))
    worst = int(np.argmax(scores[witness]))
    upper = float(scores[witness, worst])
    target_with_reserve = uniform_target - hardening_reserve
    return {
        "candidate_status": (
            "PASS_BINARY64_NEEDS_OUTWARD_REPLAY"
            if upper <= target_with_reserve
            else "FAIL_BINARY64_DISCOVERY_ONLY"
        ),
        "candidate_selector": {"kind": "witness", "witness": bank.references[witness]},
        "candidate_upper_log2": upper,
        "uniform_profile_target_log2": uniform_target,
        "hardening_reserve_bits": hardening_reserve,
        "positive_residual_bits": max(0.0, upper - target_with_reserve),
        "worst_vertex": list(vertices[worst]),
        "diagnostic_arithmetic": "binary64-discovery-only",
    }


def build_tree(
    band: tuple[int, ...],
    intervals: tuple[tuple[int, int], ...],
    counts: dict[tuple[int, int], int],
    total: int,
    bank: coordinate.AtlasBank,
    uniform_target: float,
    hardening_reserve: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    leaves: list[dict[str, Any]] = []

    def count_record(lower: int, upper: int) -> dict[str, Any]:
        return {
            "kind": "exact",
            "method": COUNT_METHOD,
            "value": str(counts[(lower, upper)]),
            "parameters": {
                "total": str(total),
                "class_count": CLASSES,
                "band_classes": list(band),
                "band_size": len(band),
                "lower_mass": str(lower),
                "upper_mass": str(upper),
                "summand": "C(m-1,q-1)*C(M-m-1,8-q)",
            },
        }

    def recurse(first: int, last: int, node_id: str) -> dict[str, Any]:
        lower = intervals[first][0]
        upper = intervals[last][1]
        if first == last:
            vertices = slab_vertices(total, CLASSES, band, lower, upper)
            diagnostic = score_slab(
                bank, vertices, total, uniform_target, hardening_reserve
            )
            node = {
                "state": "UNRESOLVED",
                "reason": "AWAITING_INDEPENDENT_OUTWARD_REPLAY",
                "cell": {
                    "kind": "ordered-band-mass-slab-v1",
                    "band_classes": list(band),
                    "lower_mass": str(lower),
                    "upper_mass": str(upper),
                },
                "count": count_record(lower, upper),
                "exact_vertices": [list(vertex) for vertex in vertices],
                "exact_vertex_count": len(vertices),
                "diagnostic": diagnostic,
            }
            leaves.append({"node_id": node_id, **node})
            return node

        middle = (first + last) // 2
        threshold = intervals[middle][1]
        coefficients, canonical_threshold, swapped = canonical_band_split(
            band, threshold, total
        )
        low_child = recurse(first, middle, node_id + ("1" if swapped else "0"))
        high_child = recurse(middle + 1, last, node_id + ("0" if swapped else "1"))
        left, right = (high_child, low_child) if swapped else (low_child, high_child)
        return {
            "state": "SPLIT",
            "cell": {
                "kind": "ordered-band-mass-slab-v1",
                "band_classes": list(band),
                "lower_mass": str(lower),
                "upper_mass": str(upper),
            },
            "count": count_record(lower, upper),
            "split": {
                "coefficients": coefficients,
                "threshold": str(canonical_threshold),
                "selection_source": "binary64 dominance ordering plus exact supplementary target masses",
                "raw_band_classes": list(band),
                "raw_band_threshold": str(threshold),
                "canonical_child_swap": swapped,
            },
            "left": left,
            "right": right,
        }

    return recurse(0, len(intervals) - 1, "r"), leaves


def all_positive_compositions(total: int, classes: int) -> Iterable[tuple[int, ...]]:
    def recurse(remaining: int, slots: int, prefix: tuple[int, ...]):
        if slots == 1:
            yield prefix + (remaining,)
            return
        for value in range(1, remaining - slots + 2):
            yield from recurse(remaining - value, slots - 1, prefix + (value,))

    if total >= classes:
        yield from recurse(total, classes, ())


def run_self_test() -> None:
    total = 17
    classes = 5
    band = (0, 3)
    q = len(band)
    minimum = q
    maximum = total - (classes - q)
    intervals = make_intervals(minimum, maximum, (5, 10))
    counts = exact_interval_counts(total, q, intervals, classes)
    brute = {interval: 0 for interval in intervals}
    for profile in all_positive_compositions(total, classes):
        mass = sum(profile[j] for j in band)
        owners = [interval for interval in intervals if interval[0] <= mass <= interval[1]]
        if len(owners) != 1:
            raise AssertionError("ordered bands do not own the profile exactly once")
        brute[owners[0]] += 1
    if counts != brute:
        raise AssertionError(f"exact count mismatch: {counts} != {brute}")
    if sum(counts.values()) != math.comb(total - 1, classes - 1):
        raise AssertionError("partition count does not equal the full-support root")
    for lower, upper in intervals:
        vertices = slab_vertices(total, classes, band, lower, upper)
        if len(vertices) > classes + 2 * q * (classes - q):
            raise AssertionError("slab vertex bound failed")
        for vertex in vertices:
            if sum(vertex) != total or min(vertex) < 1:
                raise AssertionError("invalid exact slab vertex")
            mass = sum(vertex[j] for j in band)
            if not lower <= mass <= upper:
                raise AssertionError("slab vertex violates ownership interval")

    coefficients, threshold, swapped = canonical_band_split((0, 4), 7, total)
    if coefficients != [1, 0, 0, 0, 1, 0, 0, 0, 0] or threshold != 7 or swapped:
        raise AssertionError("non-pivot band cut canonicalization failed")
    coefficients, threshold, swapped = canonical_band_split((0, 8), 7, total)
    if coefficients != [0, 1, 1, 1, 1, 1, 1, 1, 0] or threshold != 9 or not swapped:
        raise AssertionError("pivot band cut canonicalization failed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=coordinate.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=coordinate.DEFAULT_ATLAS)
    parser.add_argument(
        "--supplementary-atlas", type=Path, default=DEFAULT_SUPPLEMENTARY_ATLAS
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-boundaries", type=int, default=15)
    parser.add_argument("--hardening-reserve-bits", type=float, default=8.0)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        print("ordered-band synthetic partition tests passed")
        return
    if args.max_boundaries < 0 or args.max_boundaries > 63:
        raise ValueError("max-boundaries must lie in [0,63]")
    if not math.isfinite(args.hardening_reserve_bits) or args.hardening_reserve_bits < 0:
        raise ValueError("hardening reserve must be finite and nonnegative")

    manifest, atlas, source, bindings, manifest_digest, atlas_digest = coordinate.load_frozen_inputs(
        args.manifest, args.atlas
    )
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError("ordered-band producer manifest pin changed")
    supplementary_digest = coordinate.sha256_path(args.supplementary_atlas)
    if supplementary_digest != EXPECTED_SUPPLEMENTARY_SHA256:
        raise ValueError(
            f"supplementary atlas changed: {supplementary_digest}; "
            f"expected {EXPECTED_SUPPLEMENTARY_SHA256}"
        )
    supplementary_artifact = json.loads(args.supplementary_atlas.read_bytes())
    supplementary = dominance.load_supplementary_atlas(
        args.supplementary_atlas, manifest_digest, atlas_digest
    )
    base_bank = coordinate.build_atlas_bank(atlas, source, bindings)
    bank = dominance.combined_bank(base_bank, [supplementary])

    total = int(manifest["parameters"]["M"])
    root_census = next(
        row for row in manifest["census"]["feasible_masks"] if row["mask"] == "0x1ff"
    )
    root_count = int(root_census["count"])
    global_count = int(manifest["census"]["feasible_profile_count"])
    probability_target = float(manifest["parameters"]["probability_target_log2"].split("/")[0])
    uniform_target = probability_target - math.log2(global_count)

    choice = choose_dominance_order(bank, total)
    boundaries, boundary_diagnostic = choose_target_boundaries(
        supplementary_artifact,
        bank,
        choice.band,
        total,
        args.max_boundaries,
    )
    minimum = len(choice.band)
    maximum = total - (CLASSES - len(choice.band))
    intervals = make_intervals(minimum, maximum, boundaries)

    # Internal tree ranges are unions of adjacent leaf intervals.  Counting
    # all of them in one scan keeps runtime linear in M, independent of depth.
    all_intervals = set(intervals)
    def add_tree_ranges(first: int, last: int) -> None:
        all_intervals.add((intervals[first][0], intervals[last][1]))
        if first == last:
            return
        middle = (first + last) // 2
        add_tree_ranges(first, middle)
        add_tree_ranges(middle + 1, last)
    add_tree_ranges(0, len(intervals) - 1)
    counts = exact_interval_counts(total, len(choice.band), all_intervals)
    if counts[(minimum, maximum)] != root_count:
        raise ValueError("ordered-band root count disagrees with the frozen census")

    tree, leaves = build_tree(
        choice.band,
        intervals,
        counts,
        total,
        bank,
        uniform_target,
        args.hardening_reserve_bits,
    )
    leaf_count_sum = sum(int(leaf["count"]["value"]) for leaf in leaves)
    if leaf_count_sum != root_count:
        raise AssertionError("ordered-band leaves do not partition the full-support root")

    supplementary_source = supplementary[3]
    artifact = {
        "schema": SCHEMA,
        "status": "INCOMPLETE_DIAGNOSTIC_AWAITING_OUTWARD_REPLAY",
        "scope": "full support 0x1ff only; discovery geometry, not a coverage proof",
        "manifest_binding": {
            "path": str(args.manifest),
            "sha256": manifest_digest,
        },
        "base_atlas_binding": {
            "path": str(args.atlas),
            "sha256": atlas_digest,
            "source_id": source["source_id"],
        },
        "supplementary_atlas_binding": supplementary_source,
        "split_policy": SPLIT_POLICY,
        "geometry": {
            "kind": "single-ordered-prefix-two-band-slabs-v1",
            "class_order": list(choice.order),
            "band_classes": list(choice.band),
            "band_size": len(choice.band),
            "complement_size": CLASSES - len(choice.band),
            "feasible_band_mass": [str(minimum), str(maximum)],
            "exact_boundaries": [str(value) for value in boundaries],
            "vertex_bound": {
                "formula": "s+2*q*(s-q)",
                "value": CLASSES + 2 * len(choice.band) * (CLASSES - len(choice.band)),
            },
            "count_method": COUNT_METHOD,
            "count_time": "O(M) binomial-integer operations for all nodes",
            "count_extra_memory": "O(node_count) big integers",
        },
        "dominance_selection": {
            "binary64_discovery_only": True,
            "witness_pair_indices": list(choice.witness_pair),
            "witness_pair_bindings": [
                bank.references[index] for index in choice.witness_pair
            ],
            "root_leading_witness_indices": list(choice.root_leaders),
            "pair_alignment_count": choice.alignment,
            "threshold_integrality": (
                "every threshold is copied from an exact supplementary target profile mass; "
                "no float-to-integer quantization occurs"
            ),
            **boundary_diagnostic,
        },
        "root": {
            "support_mask": "0x1ff",
            "active_classes": list(range(CLASSES)),
            "count": str(root_count),
            "tree": tree,
        },
        "summary": {
            "boundary_count": len(boundaries),
            "leaf_count": len(leaves),
            "node_count": 2 * len(leaves) - 1,
            "leaf_count_sum": str(leaf_count_sum),
            "maximum_exact_vertex_count": max(leaf["exact_vertex_count"] for leaf in leaves),
            "binary64_candidate_passes": sum(
                leaf["diagnostic"]["candidate_status"] == "PASS_BINARY64_NEEDS_OUTWARD_REPLAY"
                for leaf in leaves
            ),
        },
    }
    digest = coordinate.write_canonical(args.output, artifact)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": digest,
                "band": list(choice.band),
                "boundaries": len(boundaries),
                "leaves": len(leaves),
                "maximum_exact_vertices": artifact["summary"]["maximum_exact_vertex_count"],
                "candidate_passes": artifact["summary"]["binary64_candidate_passes"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
