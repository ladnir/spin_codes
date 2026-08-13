#!/usr/bin/env python3
"""Build a diagnostic dominance BSP inside one frozen exact g=4 cell.

The global support mesh is not changed.  Every BSP split is an exact rational
half-space cut of the selected root simplex.  Discovery uses binary64 witness
values, but the emitted tree contains only exact hyperplanes, leaf owners, and
digests; ``certify_packet_group_g4_cell_bsp.py`` independently replays it with
the outward hardener.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import numpy as np

import probe_packet_group_g4_stellar_cover as cover


DIMENSION = 4
CLASSES = 5
SUPPORT_MASK = 31
SUPPORT = (0, 1, 2, 3, 4)
CHART_INDICES = (0, 1, 2, 3)
DEPENDENT_INDEX = 4
BSP_SCHEMA = "packet-group-g4-cell-dominance-bsp-v1"
Profile = tuple[Fraction, ...]
Constraint = tuple[tuple[Fraction, ...], Fraction]
Affine = tuple[Fraction, tuple[Fraction, ...]]


def configure_support(support_mask: int) -> None:
    """Select the affine chart for one exact positive-support stratum."""

    global DIMENSION, SUPPORT_MASK, SUPPORT, CHART_INDICES, DEPENDENT_INDEX
    support = tuple(index for index in range(CLASSES) if support_mask & (1 << index))
    if len(support) < 2:
        raise ValueError("g4 cell BSP requires a positive-dimensional support")
    SUPPORT_MASK = support_mask
    SUPPORT = support
    DIMENSION = len(support) - 1
    CHART_INDICES = support[:-1]
    DEPENDENT_INDEX = support[-1]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_profile(raw: Iterable[Any]) -> Profile:
    profile = tuple(Fraction(str(value)) for value in raw)
    if len(profile) != CLASSES or sum(profile) != cover.M:
        raise ValueError("g4 cell BSP: malformed profile")
    return profile


def solve_linear(
    matrix: Iterable[Iterable[Fraction]], right: Iterable[Fraction]
) -> tuple[Fraction, ...] | None:
    rows = [list(map(Fraction, row)) + [Fraction(value)] for row, value in zip(matrix, right)]
    n = len(rows)
    if n == 0 or any(len(row) != n + 1 for row in rows):
        raise ValueError("g4 cell BSP: malformed linear system")
    for column in range(n):
        pivot = next((row for row in range(column, n) if rows[row][column]), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(n):
            if row == column or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                left - factor * right_value
                for left, right_value in zip(rows[row], rows[column])
            ]
    return tuple(row[-1] for row in rows)


def simplex_constraints(vertices: tuple[Profile, ...]) -> tuple[Constraint, ...]:
    if len(vertices) != DIMENSION + 1:
        raise ValueError("g4 cell BSP: root has the wrong simplex dimension")
    augmented = [
        [vertex[index] for index in CHART_INDICES] + [Fraction(1)]
        for vertex in vertices
    ]
    constraints = []
    for owner in range(DIMENSION + 1):
        target = [Fraction(int(index == owner)) for index in range(DIMENSION + 1)]
        barycentric = solve_linear(augmented, target)
        if barycentric is None:
            raise ValueError("g4 cell BSP: degenerate root simplex")
        # lambda_owner(x) >= 0.
        constraints.append(
            (tuple(-value for value in barycentric[:DIMENSION]), barycentric[-1])
        )
    return tuple(constraints)


def profile_from_chart(point: tuple[Fraction, ...]) -> Profile:
    profile = [Fraction(0)] * CLASSES
    for index, value in zip(CHART_INDICES, point):
        profile[index] = value
    profile[DEPENDENT_INDEX] = Fraction(cover.M) - sum(point, Fraction(0))
    return tuple(profile)


def enumerate_vertices(constraints: tuple[Constraint, ...]) -> tuple[Profile, ...]:
    vertices: set[Profile] = set()
    for active in itertools.combinations(range(len(constraints)), DIMENSION):
        solution = solve_linear(
            [constraints[index][0] for index in active],
            [constraints[index][1] for index in active],
        )
        if solution is None:
            continue
        if all(
            sum(coefficient * value for coefficient, value in zip(left, solution)) <= right
            for left, right in constraints
        ):
            vertices.add(profile_from_chart(solution))
    return tuple(sorted(vertices))


def affine_value(affine: Affine, profile: Profile) -> Fraction:
    constant, coefficients = affine
    return constant + sum(
        coefficient * value for coefficient, value in zip(coefficients, profile)
    )


def affine_constraint(affine: Affine, nonpositive: bool) -> Constraint:
    constant, coefficients = affine
    chart = tuple(
        coefficients[index] - coefficients[DEPENDENT_INDEX]
        for index in CHART_INDICES
    )
    chart_constant = constant + coefficients[DEPENDENT_INDEX] * cover.M
    if nonpositive:
        return chart, -chart_constant
    return tuple(-value for value in chart), chart_constant


def primitive_affine(affine: Affine) -> Affine:
    values = [affine[0], *affine[1]]
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in values]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("g4 cell BSP: zero affine cut")
    integers = [value // divisor for value in integers]
    return Fraction(integers[0]), tuple(Fraction(value) for value in integers[1:])


def witness_affine(left: cover.Witness, right: cover.Witness) -> Affine | None:
    if left.subtract_normalization != right.subtract_normalization:
        return None
    # h = A_left - A_right.  The nonpositive child is owned locally by left.
    constant = Fraction.from_float(left.constant_log2) - Fraction.from_float(
        right.constant_log2
    )
    coefficients = tuple(
        Fraction.from_float(right_charge) - Fraction.from_float(left_charge)
        for left_charge, right_charge in zip(left.charge, right.charge)
    )
    if constant == 0 and not any(coefficients):
        return None
    return primitive_affine((constant, coefficients))


def serialized_vertices(vertices: tuple[Profile, ...]) -> list[list[str]]:
    return [[fraction_text(value) for value in profile] for profile in vertices]


@dataclass(frozen=True)
class CutChoice:
    affine: Affine
    negative_constraints: tuple[Constraint, ...]
    positive_constraints: tuple[Constraint, ...]
    negative_vertices: tuple[Profile, ...]
    positive_vertices: tuple[Profile, ...]
    witness_pair: tuple[str, str] | None
    predicted_worst_owner_score: float


class Builder:
    def __init__(
        self,
        witnesses: list[cover.Witness],
        target: float,
        safety_bits: float,
        max_depth: int,
        top_per_vertex: int,
    ) -> None:
        self.witnesses = witnesses
        self.evaluator = cover.Evaluator(witnesses)
        self.target = target
        self.safety_bits = safety_bits
        self.max_depth = max_depth
        self.top_per_vertex = top_per_vertex
        self.nodes: list[dict[str, Any]] = []
        self.leaf_scores: list[float] = []
        self.leaf_statuses: list[str] = []

    def values(self, vertices: tuple[Profile, ...]) -> np.ndarray:
        return np.asarray(
            [self.evaluator.evaluate(profile)[1] for profile in vertices],
            dtype=np.float64,
        )

    @staticmethod
    def best_owner(values: np.ndarray) -> tuple[int, float]:
        scores = np.max(values, axis=0)
        owner = int(np.argmin(scores))
        return owner, float(scores[owner])

    def choose_cut(
        self,
        constraints: tuple[Constraint, ...],
        vertices: tuple[Profile, ...],
        values: np.ndarray,
    ) -> CutChoice | None:
        finite = np.flatnonzero(np.all(np.isfinite(values), axis=0))
        candidates: set[int] = set()
        for vertex in range(len(vertices)):
            order = finite[np.lexsort((finite, values[vertex, finite]))]
            candidates.update(map(int, order[: self.top_per_vertex]))

        choices: list[CutChoice] = []
        for left_index, right_index in itertools.combinations(sorted(candidates), 2):
            affine = witness_affine(
                self.witnesses[left_index], self.witnesses[right_index]
            )
            if affine is None:
                continue
            signs = [affine_value(affine, profile) for profile in vertices]
            if not any(value < 0 for value in signs) or not any(value > 0 for value in signs):
                continue
            negative_constraints = constraints + (affine_constraint(affine, True),)
            positive_constraints = constraints + (affine_constraint(affine, False),)
            negative_vertices = enumerate_vertices(negative_constraints)
            positive_vertices = enumerate_vertices(positive_constraints)
            if len(negative_vertices) < DIMENSION + 1 or len(positive_vertices) < DIMENSION + 1:
                continue
            negative_score = self.best_owner(self.values(negative_vertices))[1]
            positive_score = self.best_owner(self.values(positive_vertices))[1]
            choices.append(
                CutChoice(
                    affine=affine,
                    negative_constraints=negative_constraints,
                    positive_constraints=positive_constraints,
                    negative_vertices=negative_vertices,
                    positive_vertices=positive_vertices,
                    witness_pair=(
                        self.witnesses[left_index].reference,
                        self.witnesses[right_index].reference,
                    ),
                    predicted_worst_owner_score=max(negative_score, positive_score),
                )
            )
        if choices:
            return min(
                choices,
                key=lambda choice: (
                    choice.predicted_worst_owner_score,
                    abs(len(choice.negative_vertices) - len(choice.positive_vertices)),
                    choice.witness_pair,
                ),
            )

        # Geometry-only fallback.  It remains sound because the verifier does
        # not rely on a cut being a dominance plane.
        coordinate = max(
            SUPPORT,
            key=lambda index: max(p[index] for p in vertices) - min(p[index] for p in vertices),
        )
        low = min(profile[coordinate] for profile in vertices)
        high = max(profile[coordinate] for profile in vertices)
        if low == high:
            return None
        coefficients = [Fraction(0)] * CLASSES
        coefficients[coordinate] = 1
        affine = primitive_affine((-(low + high) / 2, tuple(coefficients)))
        negative_constraints = constraints + (affine_constraint(affine, True),)
        positive_constraints = constraints + (affine_constraint(affine, False),)
        negative_vertices = enumerate_vertices(negative_constraints)
        positive_vertices = enumerate_vertices(positive_constraints)
        if len(negative_vertices) < DIMENSION + 1 or len(positive_vertices) < DIMENSION + 1:
            return None
        return CutChoice(
            affine=affine,
            negative_constraints=negative_constraints,
            positive_constraints=positive_constraints,
            negative_vertices=negative_vertices,
            positive_vertices=positive_vertices,
            witness_pair=None,
            predicted_worst_owner_score=max(
                self.best_owner(self.values(negative_vertices))[1],
                self.best_owner(self.values(positive_vertices))[1],
            ),
        )

    def add_leaf(
        self,
        depth: int,
        vertices: tuple[Profile, ...],
        owner: int,
        owner_score: float,
        envelope_maximum: float,
        status: str,
    ) -> str:
        identifier = f"n{len(self.nodes):06d}"
        rendered = serialized_vertices(vertices)
        self.nodes.append(
            {
                "node_id": identifier,
                "kind": "leaf",
                "depth": depth,
                "status": status,
                "owner_witness": self.witnesses[owner].reference,
                "vertex_count": len(vertices),
                "vertices_sha256": compact_hash(rendered),
                "diagnostic_owner_maximum_log2": owner_score,
                "diagnostic_envelope_maximum_at_vertices_log2": envelope_maximum,
            }
        )
        self.leaf_scores.append(owner_score)
        self.leaf_statuses.append(status)
        return identifier

    def build_node(
        self,
        constraints: tuple[Constraint, ...],
        depth: int,
    ) -> str:
        vertices = enumerate_vertices(constraints)
        if len(vertices) < DIMENSION + 1:
            raise RuntimeError("g4 cell BSP: a recursive node lost full dimension")
        values = self.values(vertices)
        owner, owner_score = self.best_owner(values)
        envelope_by_vertex = np.min(values, axis=1)
        envelope_maximum = float(np.max(envelope_by_vertex))
        threshold = self.target - self.safety_bits
        if owner_score <= threshold:
            return self.add_leaf(
                depth, vertices, owner, owner_score, envelope_maximum, "covered"
            )
        if envelope_maximum > threshold:
            return self.add_leaf(
                depth,
                vertices,
                owner,
                owner_score,
                envelope_maximum,
                "pointwise_atlas_failure",
            )
        if depth >= self.max_depth:
            return self.add_leaf(
                depth, vertices, owner, owner_score, envelope_maximum, "depth_cap"
            )
        choice = self.choose_cut(constraints, vertices, values)
        if choice is None:
            return self.add_leaf(
                depth, vertices, owner, owner_score, envelope_maximum, "no_cut"
            )

        identifier = f"n{len(self.nodes):06d}"
        row: dict[str, Any] = {
            "node_id": identifier,
            "kind": "split",
            "depth": depth,
            "hyperplane": {
                "constant": fraction_text(choice.affine[0]),
                "coefficients": [fraction_text(value) for value in choice.affine[1]],
                "nonpositive_child_semantics": "constant+dot(coefficients,profile)<=0",
            },
            "discovery_witness_pair": (
                None if choice.witness_pair is None else list(choice.witness_pair)
            ),
            "diagnostic_predicted_worst_owner_score": choice.predicted_worst_owner_score,
        }
        self.nodes.append(row)
        negative = self.build_node(choice.negative_constraints, depth + 1)
        positive = self.build_node(choice.positive_constraints, depth + 1)
        row["nonpositive_child"] = negative
        row["nonnegative_child"] = positive
        return identifier


def select_cell(stratum: dict[str, Any], requested: str | None) -> dict[str, Any]:
    rows = stratum.get("covered_cell_ledger", [])
    if requested is not None:
        return next(row for row in rows if str(row.get("cell_id")) == requested)
    return max(rows, key=lambda row: float(row["cell_local_contribution_log2_upperish"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--support-mask", type=lambda value: int(value, 0), default=31)
    parser.add_argument("--cell-id")
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--top-witnesses-per-vertex", type=int, default=2)
    parser.add_argument("--safety-bits", type=float, default=1e-5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_depth < 0 or args.top_witnesses_per_vertex < 1 or args.safety_bits < 0:
        parser.error("invalid BSP depth, candidate count, or safety margin")
    configure_support(args.support_mask)

    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    if int(ledger.get("group_bits", -1)) != cover.GROUP_BITS:
        raise ValueError("g4 cell BSP: source ledger is not g=4")
    stratum = next(
        row for row in ledger["support_strata"] if int(row["support_mask"]) == SUPPORT_MASK
    )
    cell = select_cell(stratum, args.cell_id)
    anchors = tuple(parse_profile(row) for row in stratum["anchor_profiles"])
    indices = tuple(map(int, cell["anchor_indices"]))
    root_vertices = tuple(anchors[index] for index in indices)
    constraints = simplex_constraints(root_vertices)
    if enumerate_vertices(constraints) != tuple(sorted(root_vertices)):
        raise RuntimeError("g4 cell BSP: root half-spaces do not reconstruct the simplex")

    witnesses, sources, _atlas_anchors, duplicates = cover.load_witnesses(args.atlas)
    target = float(cell["target_log2"])
    builder = Builder(
        witnesses,
        target,
        args.safety_bits,
        args.max_depth,
        args.top_witnesses_per_vertex,
    )
    root_node = builder.build_node(constraints, 0)
    cell_maximum = max(builder.leaf_scores)
    count = int(cell["integer_profile_count_upper"])
    contribution = math.nextafter(cell_maximum + math.log2(count), math.inf)
    leaf_counts = {
        status: builder.leaf_statuses.count(status) for status in sorted(set(builder.leaf_statuses))
    }
    report = {
        "schema": BSP_SCHEMA,
        "status": "DIAGNOSTIC_G4_SINGLE_CELL_DOMINANCE_BSP",
        "group_bits": cover.GROUP_BITS,
        "source_ledger": {
            "path": str(args.ledger.resolve()),
            "sha256": file_sha256(args.ledger),
        },
        "sources": sources,
        "duplicate_atlas_anchor_rows": duplicates,
        "support_mask": SUPPORT_MASK,
        "cell_id": str(cell["cell_id"]),
        "root_anchor_indices": list(indices),
        "root_anchor_profiles": serialized_vertices(root_vertices),
        "root_integer_profile_count_upper": count,
        "uniform_target_log2": target,
        "safety_bits": args.safety_bits,
        "root_node": root_node,
        "nodes": builder.nodes,
        "complete_exact_partition_claimed": True,
        "discovery": {
            "method": "greedy exact-dyadic dominance BSP with coordinate fallback",
            "max_depth": args.max_depth,
            "top_witnesses_per_vertex": args.top_witnesses_per_vertex,
            "witnesses": len(witnesses),
        },
        "diagnostic": {
            "root_maximum_vertex_log2": float(cell["maximum_vertex_log2"]),
            "root_union_contribution_log2": float(
                cell["cell_local_contribution_log2_upperish"]
            ),
            "bsp_cell_maximum_log2": cell_maximum,
            "bsp_cell_union_contribution_log2": contribution,
            "improvement_bits": float(cell["maximum_vertex_log2"]) - cell_maximum,
            "leaf_count": len(builder.leaf_scores),
            "leaf_status_counts": leaf_counts,
            "all_leaves_pass_uniform_target": all(
                status == "covered" for status in builder.leaf_statuses
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["diagnostic"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
