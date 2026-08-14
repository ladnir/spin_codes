#!/usr/bin/env python3
"""Run one bounded adaptive cumulative-prefix diagnostic wave for g=8."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np

import evaluate_packet_group_g8_cumulative_boxes as h2
import plan_packet_group_g8_lowdim_root_cover as minimax
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import prototype_packet_group_g8_cumulative_boxes as geometry


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARENT = ROOT / "out" / "g8_full_support_cumulative_h2_diagnostic.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_adaptive_cumulative_wave1.json"
PARENT_SHA256 = "47308693e73546e249f21c1a2edef4a820c4017820b6d9f8f670e95d720151b3"
SCHEMA = "permute-conv.packet-group-g8-adaptive-cumulative-wave.v1"
PARENT_CELLS = (4, 162, 10, 161, 20, 35, 160, 163)
DIMENSION = 8
VERTEX_CAP_PER_CHILD = 320
VERTEX_CAP_WAVE = 5120


def compositions(total: int):
    if total == 0:
        yield ()
        return
    for first in range(1, total + 1):
        for rest in compositions(total - first):
            yield (first, *rest)


def close_bounds(
    lower: tuple[int, ...], upper: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    local_lower = list(lower)
    local_upper = list(upper)
    for index in range(1, len(local_lower)):
        local_lower[index] = max(local_lower[index], local_lower[index - 1])
    for index in range(len(local_upper) - 2, -1, -1):
        local_upper[index] = min(local_upper[index], local_upper[index + 1])
    if any(low > high for low, high in zip(local_lower, local_upper)):
        return None
    return tuple(local_lower), tuple(local_upper)


def child_bounds(
    lower: tuple[int, ...], upper: tuple[int, ...], k: int, t: int, left: bool
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    local_lower = list(lower)
    local_upper = list(upper)
    if left:
        local_upper[k] = min(local_upper[k], t)
    else:
        local_lower[k] = max(local_lower[k], t + 1)
    return close_bounds(tuple(local_lower), tuple(local_upper))


def exact_chain_count(lower: tuple[int, ...], upper: tuple[int, ...]) -> int:
    """Count bounded nondecreasing integer chains by exact prefix sums."""

    closed = close_bounds(lower, upper)
    if closed is None:
        return 0
    lower, upper = closed
    previous_low = lower[0]
    previous_high = upper[0]
    previous = [1] * (previous_high - previous_low + 1)
    for low, high in zip(lower[1:], upper[1:]):
        prefixes = list(itertools.accumulate(previous))
        total = prefixes[-1]
        if low > previous_high:
            current = [total] * (high - low + 1)
        else:
            start = low - previous_low
            current = prefixes[start:]
            current.extend([total] * (high - previous_high))
        previous = current
        previous_low = low
        previous_high = high
    return sum(previous)


def chain_vertices(
    lower: tuple[int, ...], upper: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Enumerate exact vertices from anchored constant blocks."""

    closed = close_bounds(lower, upper)
    if closed is None:
        return ()
    lower, upper = closed
    result: set[tuple[int, ...]] = set()
    for lengths in compositions(len(lower)):
        blocks = []
        start = 0
        for length in lengths:
            end = start + length
            anchors = sorted(
                {value for index in range(start, end) for value in (lower[index], upper[index])}
            )
            blocks.append((start, end, anchors))
            start = end
        for chosen in itertools.product(*(block[2] for block in blocks)):
            row = []
            for (_start, _end, _anchors), value in zip(blocks, chosen):
                row.extend([value] * (_end - _start))
            vertex = tuple(row)
            if (
                all(low <= value <= high for value, low, high in zip(vertex, lower, upper))
                and all(left <= right for left, right in zip(vertex, vertex[1:]))
            ):
                result.add(vertex)
    return tuple(sorted(result))


def profile_from_chain(chain: tuple[int, ...]) -> tuple[int, ...]:
    shifted = []
    previous = 0
    for value in chain:
        shifted.append(value - previous)
        previous = value
    shifted.append(geometry.MASS - 9 - previous)
    return tuple(value + 1 for value in shifted)


def chain_from_profile(profile: tuple[int, ...]) -> tuple[int, ...]:
    running = 0
    result = []
    for value in profile[:-1]:
        running += value - 1
        result.append(running)
    return tuple(result)


def parse_profile(row: list[Any]) -> tuple[int, ...]:
    values = tuple(Fraction(value) for value in row)
    if any(value.denominator != 1 for value in values):
        raise ValueError("full-support cumulative parent has a nonintegral vertex")
    return tuple(value.numerator for value in values)


def serialize_vertices(vertices: tuple[tuple[int, ...], ...]) -> list[list[int]]:
    return [list(profile_from_chain(vertex)) for vertex in vertices]


def vertex_digest(vertices: tuple[tuple[int, ...], ...]) -> str:
    return coordinate.sha256_bytes(
        coordinate.canonical_bytes(serialize_vertices(vertices))
    )


def score_matrix(
    vertices: tuple[tuple[int, ...], ...], bank: coordinate.AtlasBank
) -> np.ndarray:
    profiles = np.asarray([profile_from_chain(vertex) for vertex in vertices], dtype=np.float64)
    normalizations = coordinate.diagnostic_normalization(profiles, geometry.MASS)
    return (
        bank.constants[:, None] - bank.charges @ profiles.T - normalizations[None, :]
    ).T


def reference_key(reference: dict[str, Any]) -> tuple[Any, ...]:
    return (
        reference["source_sha256"],
        int(reference["row"]),
        reference["row_sha256"],
    )


def leading_witnesses(values: np.ndarray, bank: coordinate.AtlasBank) -> tuple[int, ...]:
    result = []
    for row in values:
        minimum = float(np.min(row))
        tied = np.flatnonzero(row == minimum)
        result.append(min(map(int, tied), key=lambda index: reference_key(bank.references[index])))
    return tuple(result)


def transition_candidates(
    vertices: tuple[tuple[int, ...], ...],
    values: np.ndarray,
    leaders: tuple[int, ...],
    bank: coordinate.AtlasBank,
) -> dict[tuple[int, int], dict[str, Any]]:
    proposals: dict[tuple[int, int], dict[str, Any]] = {}
    profiles = tuple(profile_from_chain(vertex) for vertex in vertices)
    for k in range(DIMENSION):
        groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for index, vertex in enumerate(vertices):
            groups[vertex[:k] + vertex[k + 1 :]].append(index)
        for indices in groups.values():
            indices.sort(key=lambda index: (vertices[index][k], index))
            for left_index, right_index in zip(indices, indices[1:]):
                left_value = vertices[left_index][k]
                right_value = vertices[right_index][k]
                p = leaders[left_index]
                q = leaders[right_index]
                if left_value >= right_value or p == q:
                    continue
                delta_constant = bank.constants[p] - bank.constants[q]
                delta_charge = bank.charges[p] - bank.charges[q]
                d_left = float(delta_constant - delta_charge @ np.asarray(profiles[left_index]))
                d_right = float(delta_constant - delta_charge @ np.asarray(profiles[right_index]))
                thresholds: set[int] = set()
                crossing = None
                if (
                    math.isfinite(d_left)
                    and math.isfinite(d_right)
                    and d_left != d_right
                    and ((d_left <= 0.0 <= d_right) or (d_right <= 0.0 <= d_left))
                ):
                    ratio = -d_left / (d_right - d_left)
                    crossing = left_value + ratio * (right_value - left_value)
                    if math.isfinite(crossing):
                        thresholds.update((math.floor(crossing), math.ceil(crossing) - 1))
                if not thresholds:
                    thresholds.add((left_value + right_value) // 2)
                for threshold in thresholds:
                    if not left_value <= threshold < right_value:
                        continue
                    diagnostic = {
                        "coordinate": k,
                        "threshold": str(threshold),
                        "left_chain_vertex": list(vertices[left_index]),
                        "right_chain_vertex": list(vertices[right_index]),
                        "left_leading_witness": bank.references[p],
                        "right_leading_witness": bank.references[q],
                        "affine_difference_left": d_left,
                        "affine_difference_right": d_right,
                        "crossing_coordinate_binary64": crossing,
                    }
                    key = (k, threshold)
                    encoded = coordinate.canonical_bytes(diagnostic)
                    old = proposals.get(key)
                    if old is None or encoded < coordinate.canonical_bytes(old):
                        proposals[key] = diagnostic
    return proposals


def child_diagnostic(
    bounds: tuple[tuple[int, ...], tuple[int, ...]],
    bank: coordinate.AtlasBank,
    denominator: int,
) -> dict[str, Any] | None:
    lower, upper = bounds
    vertices = chain_vertices(lower, upper)
    if not vertices or len(vertices) > VERTEX_CAP_PER_CHILD:
        return None
    count = exact_chain_count(lower, upper)
    if count <= 0:
        return None
    values = score_matrix(vertices, bank)
    mixed = minimax.optimize_minimax_mixture(
        values,
        minimax.UNIFORM_TARGET_LOG2,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-9,
        weight_tolerance=1e-12,
        denominator=denominator,
    )
    upper_score = float(mixed["score"])
    worst_index = int(np.argmax(np.asarray(mixed["vertex_values"])))
    return {
        "lower": list(lower),
        "upper": list(upper),
        "exact_count": str(count),
        "count": {
            "kind": "exact",
            "method": "bounded-monotone-chain-prefix-dp-v1",
            "parameters": {
                "lower": list(lower),
                "upper": list(upper),
                "free_mass": str(geometry.MASS - 9),
            },
            "value": str(count),
        },
        "exact_vertex_count": len(vertices),
        "exact_vertices": serialize_vertices(vertices),
        "exact_vertex_sha256": vertex_digest(vertices),
        "candidate_upper_log2": upper_score,
        "contribution_log2": math.log2(count) + upper_score,
        "worst_vertex": list(profile_from_chain(vertices[worst_index])),
        "selector": h2.selector(mixed, bank),
        "component_count": len(mixed["indices"]),
        "candidate_status": (
            "PASS_BINARY64_REQUIRES_OUTWARD_REPLAY"
            if upper_score <= minimax.UNIFORM_TARGET_LOG2
            else "FAIL_BINARY64"
        ),
        "diagnostic": {
            "binary64_lp_upper_log2": float(mixed["binary64_lp_score"]),
            "rationalization_loss_bits": float(mixed["rationalization_loss_bits"]),
            "column_generation_rounds": int(mixed["column_generation_rounds"]),
            "full_column_scans": int(mixed["full_column_scans"]),
        },
    }


def parent_bounds(cell: dict[str, Any]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    intervals = tuple(tuple(map(int, row)) for row in cell["bin_intervals"])
    return tuple(row[0] for row in intervals), tuple(row[1] for row in intervals)


def evaluate_parent(
    cell: dict[str, Any], bank: coordinate.AtlasBank, denominator: int
) -> dict[str, Any]:
    lower, upper = parent_bounds(cell)
    vertices = tuple(chain_from_profile(parse_profile(row)) for row in cell["vertices"])
    values = score_matrix(vertices, bank)
    leaders = leading_witnesses(values, bank)
    proposals = transition_candidates(vertices, values, leaders, bank)
    parent_count = int(cell["exact_count"])
    parent_contribution = float(cell["mixture"]["contribution_log2"])
    candidates = []
    rejected_vertex_cap = 0
    for (k, threshold), transition in sorted(proposals.items()):
        left_bounds = child_bounds(lower, upper, k, threshold, True)
        right_bounds = child_bounds(lower, upper, k, threshold, False)
        if left_bounds is None or right_bounds is None:
            continue
        left = child_diagnostic(left_bounds, bank, denominator)
        right = child_diagnostic(right_bounds, bank, denominator)
        if left is None or right is None:
            rejected_vertex_cap += 1
            continue
        if int(left["exact_count"]) + int(right["exact_count"]) != parent_count:
            raise RuntimeError("adaptive child counts do not sum to the parent")
        maximum_contribution = max(left["contribution_log2"], right["contribution_log2"])
        maximum_bound = max(left["candidate_upper_log2"], right["candidate_upper_log2"])
        minimum_count = min(int(left["exact_count"]), int(right["exact_count"]))
        total_vertices = left["exact_vertex_count"] + right["exact_vertex_count"]
        improvement = parent_contribution - maximum_contribution
        candidates.append(
            {
                "key": (
                    maximum_contribution,
                    maximum_bound,
                    -minimum_count,
                    total_vertices,
                    k,
                    threshold,
                ),
                "improvement": improvement,
                "k": k,
                "threshold": threshold,
                "transition": transition,
                "left": left,
                "right": right,
            }
        )
    positive = [candidate for candidate in candidates if candidate["improvement"] > 0.0]
    result = {
        "parent_cell": int(cell["cell"]),
        "parent_bin_indices": cell["bin_indices"],
        "parent_exact_count": cell["exact_count"],
        "parent_candidate_upper_log2": float(cell["mixture"]["candidate_upper_log2"]),
        "parent_contribution_log2": parent_contribution,
        "transition_candidate_count": len(proposals),
        "evaluated_candidate_count": len(candidates),
        "rejected_candidate_count": len(proposals) - len(candidates),
        "rejected_by_empty_or_vertex_cap": rejected_vertex_cap,
    }
    if not positive:
        result.update(
            {
                "state": "UNRESOLVED",
                "reason": "NO_POSITIVE_DIAGNOSTIC_IMPROVEMENT",
                "improvement_bits": 0.0,
            }
        )
        return result
    chosen = min(positive, key=lambda candidate: candidate["key"])
    coefficients = [1 if index <= chosen["k"] else 0 for index in range(9)]
    profile_threshold = chosen["threshold"] + chosen["k"] + 1
    result.update(
        {
            "state": "SPLIT_DIAGNOSTIC_ONLY",
            "improvement_bits": chosen["improvement"],
            "maximum_child_contribution_log2": chosen["key"][0],
            "maximum_child_bound_log2": chosen["key"][1],
            "split": {
                "k": chosen["k"],
                "t": str(chosen["threshold"]),
                "coefficients": coefficients,
                "threshold": str(profile_threshold),
            },
            "diagnostic": {"selected_transition": chosen["transition"]},
            "left": chosen["left"],
            "right": chosen["right"],
        }
    )
    return result


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    parent_digest = coordinate.sha256_path(args.parent)
    if parent_digest != PARENT_SHA256:
        raise ValueError("h2 parent diagnostic digest changed")
    parent = json.loads(args.parent.read_bytes())
    if parent.get("schema") != h2.SCHEMA:
        raise ValueError("unexpected h2 parent schema")
    bank, sources, bindings = h2.build_bank(
        args.manifest, args.atlas, args.supplementary
    )
    selected = []
    for index in PARENT_CELLS:
        cell = parent["cells"][index]
        if int(cell["cell"]) != index:
            raise RuntimeError("h2 cell index is not canonical")
        selected.append(evaluate_parent(cell, bank, args.mixture_denominator))
    improvements = [float(row["improvement_bits"]) for row in selected]
    gate_count = sum(value >= 16.0 for value in improvements)
    child_vertices = sum(
        int(row[side]["exact_vertex_count"])
        for row in selected
        if row["state"] == "SPLIT_DIAGNOSTIC_ONLY"
        for side in ("left", "right")
    )
    if child_vertices > VERTEX_CAP_WAVE:
        raise RuntimeError("adaptive wave exceeded its conservative vertex cap")
    return {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_BINARY64_ADAPTIVE_WAVE_REQUIRES_OUTWARD_REPLAY",
        **bindings,
        "parent_artifact": {
            "path": str(args.parent),
            "sha256": parent_digest,
            "schema": parent["schema"],
        },
        "witness_sources": sources,
        "parent_cells": list(PARENT_CELLS),
        "unchanged_parent_cells": [
            index for index in range(len(parent["cells"])) if index not in PARENT_CELLS
        ],
        "selection_policy": "witness-transition-prefix-split-v1",
        "geometry": {
            "method": "bounded-monotone-chain-intervals-v1",
            "count_method": "bounded-monotone-chain-prefix-dp-v1",
            "vertex_method": "anchored-constant-block-enumeration-v1",
            "vertex_cap_per_child": VERTEX_CAP_PER_CHILD,
            "vertex_cap_wave": VERTEX_CAP_WAVE,
        },
        "refinements": selected,
        "stop_gate": {
            "required_improvement_bits": 16.0,
            "required_parent_count": 4,
            "observed_parent_count": gate_count,
            "passes": gate_count >= 4,
        },
        "summary": {
            "selected_parent_count": len(selected),
            "accepted_split_count": sum(
                row["state"] == "SPLIT_DIAGNOSTIC_ONLY" for row in selected
            ),
            "actual_child_vertex_incidences": child_vertices,
            "minimum_improvement_bits": min(improvements),
            "maximum_improvement_bits": max(improvements),
            "runtime_seconds": time.perf_counter() - started,
        },
        "scope_limit": (
            "Exact integer ownership, counts, and vertices only. Binary64 transition "
            "selection and rational-mixture scoring require independent outward replay. "
            "This artifact is not a certificate, completion claim, or probability proof."
        ),
    }


def run_self_test() -> None:
    parent = geometry.build(tuple(range(9)), 2, 10000)["cells"][4]
    lower, upper = parent_bounds(parent)
    expected_vertices = {
        chain_from_profile(parse_profile(row)) for row in parent["vertices"]
    }
    observed_vertices = set(chain_vertices(lower, upper))
    if observed_vertices != expected_vertices:
        raise SystemExit("adaptive anchored-block vertex self-test failed")
    count = exact_chain_count(lower, upper)
    if count != int(parent["exact_count"]):
        raise SystemExit("adaptive monotone-chain count self-test failed")
    k = 6
    threshold = (lower[k] + upper[k]) // 2
    left = child_bounds(lower, upper, k, threshold, True)
    right = child_bounds(lower, upper, k, threshold, False)
    if (
        left is None
        or right is None
        or exact_chain_count(*left) + exact_chain_count(*right) != count
        or max(len(chain_vertices(*left)), len(chain_vertices(*right))) > VERTEX_CAP_PER_CHILD
    ):
        raise SystemExit("adaptive split partition self-test failed")
    print("parent_cell=4")
    print(f"parent_vertices={len(observed_vertices)}")
    print(f"left_vertices={len(chain_vertices(*left))}")
    print(f"right_vertices={len(chain_vertices(*right))}")
    print("exact_child_count_sum=PASS")
    print("status=PASS_G8_ADAPTIVE_CUMULATIVE_WAVE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=h2.DEFAULT_ATLAS)
    parser.add_argument("--supplementary", type=Path, default=h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--parent", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.mixture_denominator < 1:
        parser.error("mixture denominator must be positive")
    report = build(args)
    digest = coordinate.write_canonical(args.output, report)
    print(f"parents={report['summary']['selected_parent_count']}")
    print(f"accepted_splits={report['summary']['accepted_split_count']}")
    for row in report["refinements"]:
        print(
            f"cell={row['parent_cell']} improvement_bits={row['improvement_bits']:.12f} "
            f"state={row['state']}"
        )
    gate = report["stop_gate"]
    print(
        f"stop_gate_count={gate['observed_parent_count']} "
        f"stop_gate_pass={str(gate['passes']).upper()}"
    )
    print(f"child_vertex_incidences={report['summary']['actual_child_vertex_incidences']}")
    print(f"runtime_seconds={report['summary']['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=DIAGNOSTIC_BINARY64_ADAPTIVE_WAVE_REQUIRES_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
