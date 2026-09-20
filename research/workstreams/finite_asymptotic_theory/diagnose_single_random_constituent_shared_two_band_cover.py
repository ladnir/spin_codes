#!/usr/bin/env python3
"""Adaptive convex cover of every dense two-band composition."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from evaluate_single_random_constituent_highprob_bands import band_majorant
from evaluate_single_random_constituent_highprob_renyi import LOG2, spectrum_caps
from evaluate_single_random_constituent_shared_two_band_dense import (
    evaluate_with_witness,
    optimize_counts,
)


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_s22.json"


def counts_from_point(point: tuple[float, float], outer_rows: int) -> tuple[float, float, float]:
    inactive, defects = point
    return inactive, defects, outer_rows - inactive - defects


def midpoint(
    left: tuple[float, float], right: tuple[float, float]
) -> tuple[float, float]:
    return ((left[0] + right[0]) / 2.0, (left[1] + right[1]) / 2.0)


def split_longest(
    triangle: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
) -> tuple[
    tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
]:
    vertices = list(triangle)
    edges = [(0, 1), (1, 2), (2, 0)]
    first, second = max(
        edges,
        key=lambda pair: (
            (vertices[pair[0]][0] - vertices[pair[1]][0]) ** 2
            + (vertices[pair[0]][1] - vertices[pair[1]][1]) ** 2
        ),
    )
    third = 3 - first - second
    middle = midpoint(vertices[first], vertices[second])
    return (
        (vertices[first], middle, vertices[third]),
        (middle, vertices[second], vertices[third]),
    )


def point_in_triangle(
    point: tuple[int, int],
    triangle: list[list[float]] | tuple[tuple[float, float], ...],
) -> bool:
    x, y = point
    (x1, y1), (x2, y2), (x3, y3) = triangle
    denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    if denominator == 0.0:
        return False
    first = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator
    second = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator
    third = 1.0 - first - second
    return min(first, second, third) >= -1e-12


def integer_points(
    triangle: list[list[float]] | tuple[tuple[float, float], ...],
) -> set[tuple[int, int]]:
    minimum_inactive = math.ceil(min(vertex[0] for vertex in triangle))
    maximum_inactive = math.floor(max(vertex[0] for vertex in triangle))
    minimum_defects = math.ceil(min(vertex[1] for vertex in triangle))
    maximum_defects = math.floor(max(vertex[1] for vertex in triangle))
    return {
        (inactive, defects)
        for inactive in range(minimum_inactive, maximum_inactive + 1)
        for defects in range(minimum_defects, maximum_defects + 1)
        if point_in_triangle((inactive, defects), triangle)
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-occupation", type=int, default=160)
    parser.add_argument("--target-pointwise-margin", type=float, default=80.0)
    parser.add_argument("--maximum-depth", type=int, default=18)
    parser.add_argument("--maximum-cells", type=int, default=20000)
    parser.add_argument(
        "--resolve-unresolved-lattice",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--low-upper", type=int, default=79)
    parser.add_argument("--central-upper", type=int, default=432)
    parser.add_argument("--prior", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")
    outer_rows = args.output_bits // args.block_bits
    if not 1 <= args.minimum_occupation <= outer_rows:
        parser.error("invalid minimum occupation")
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator
    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=args.block_bits,
        dimension=args.dimension,
    )
    support = [weight for weight in range(1, args.block_bits + 1) if int(caps[weight])]
    low = band_majorant(caps, min(support), args.low_upper, block_bits=args.block_bits)
    central = band_majorant(
        caps, args.low_upper + 1, args.central_upper, block_bits=args.block_bits
    )
    high = band_majorant(
        caps, args.central_upper + 1, max(support), block_bits=args.block_bits
    )
    defect_probability = min(
        float(low["value_probability"]),
        1.0 - float(high["value_probability"]),
    )
    defect_log_majorant = math.log(2.0) + max(
        float(low["log_majorant"]), float(high["log_majorant"])
    )
    values = np.asarray([0.0, defect_probability, float(central["value_probability"])])
    log_majorants = np.asarray(
        [0.0, defect_log_majorant, float(central["log_majorant"])]
    )

    inactive_maximum = outer_rows - args.minimum_occupation
    all_central = (0.0, 0.0)
    all_defect = (0.0, float(outer_rows))
    sparse_defect = (float(inactive_maximum), float(args.minimum_occupation))
    sparse_central = (float(inactive_maximum), 0.0)
    accepted: list[dict[str, object]] = []
    if args.prior is None:
        pending = [
            ((all_central, all_defect, sparse_defect), 0),
            ((all_central, sparse_defect, sparse_central), 0),
        ]
    else:
        prior = json.loads(args.prior.read_text(encoding="utf-8"))
        accepted = list(prior["accepted_cells"])
        pending = [
            (
                tuple(tuple(float(value) for value in vertex) for vertex in cell["vertices"]),
                int(cell["depth"]),
            )
            for cell in prior["unresolved_cells"]
        ]
    unresolved: list[dict[str, object]] = []
    processed = 0

    while pending:
        triangle, depth = pending.pop()
        processed += 1
        if processed > args.maximum_cells:
            unresolved.extend(
                {"vertices": cell[0], "depth": cell[1], "reason": "cell budget"}
                for cell in pending
            )
            break
        centroid = (
            sum(vertex[0] for vertex in triangle) / 3.0,
            sum(vertex[1] for vertex in triangle) / 3.0,
        )
        center_counts = counts_from_point(centroid, outer_rows)
        witness = optimize_counts(
            center_counts,
            values=values,
            log_majorants=log_majorants,
            block_bits=args.block_bits,
            distance_cutoff=distance_cutoff,
            memory_bits=args.memory_bits,
            output_bits=args.output_bits,
            temperatures=(0.6, 0.9),
            scales=(1.5, 4.0),
            maximum_iterations=500,
        )
        probabilities = np.asarray(witness["reference_type_probabilities"])
        surprisal = float(witness["surprisal"])
        vertex_rows = []
        for vertex in triangle:
            counts = counts_from_point(vertex, outer_rows)
            value, raw, _ = evaluate_with_witness(
                counts,
                probabilities,
                surprisal,
                values=values,
                log_majorants=log_majorants,
                block_bits=args.block_bits,
                distance_cutoff=distance_cutoff,
                memory_bits=args.memory_bits,
                output_bits=args.output_bits,
            )
            vertex_rows.append(
                {
                    "counts": list(counts),
                    "log2_upper": value / LOG2,
                    "margin_bits": -value / LOG2,
                    "raw_reference_log2_upper": raw / LOG2,
                }
            )
        passes = all(
            row["raw_reference_log2_upper"] < 0.0
            and row["margin_bits"] >= args.target_pointwise_margin
            for row in vertex_rows
        )
        if passes:
            accepted.append(
                {
                    "vertices": [list(vertex) for vertex in triangle],
                    "depth": depth,
                    "witness_probabilities": probabilities.tolist(),
                    "witness_surprisal": surprisal,
                    "vertices_evaluated": vertex_rows,
                }
            )
        elif depth >= args.maximum_depth:
            unresolved.append(
                {
                    "vertices": [list(vertex) for vertex in triangle],
                    "depth": depth,
                    "reason": "depth limit",
                    "vertices_evaluated": vertex_rows,
                }
            )
        else:
            left, right = split_longest(triangle)
            pending.append((right, depth + 1))
            pending.append((left, depth + 1))
        if processed % 10 == 0:
            print(
                f"processed,{processed},accepted,{len(accepted)},pending,{len(pending)},unresolved,{len(unresolved)}",
                flush=True,
            )

    exception_rows = []
    unresolved_points = sorted(
        {
            point
            for cell in unresolved
            for point in integer_points(cell["vertices"])
        }
    )
    if args.resolve_unresolved_lattice:
        for index, (inactive, defects) in enumerate(unresolved_points, start=1):
            counts = (inactive, defects, outer_rows - inactive - defects)
            row = optimize_counts(
                counts,
                values=values,
                log_majorants=log_majorants,
                block_bits=args.block_bits,
                distance_cutoff=distance_cutoff,
                memory_bits=args.memory_bits,
                output_bits=args.output_bits,
                temperatures=(0.5, 0.75, 1.0),
                scales=(1.4, 3.0, 6.0),
                maximum_iterations=1000,
            )
            exception_rows.append(row)
            print(
                f"exception,{index},{len(unresolved_points)},counts,{counts},margin,{row['margin_bits']:.9f}",
                flush=True,
            )

    lattice_exceptions_pass = (
        args.resolve_unresolved_lattice
        and len(exception_rows) == len(unresolved_points)
        and all(
            bool(row["raw_negative"])
            and float(row["margin_bits"]) >= args.target_pointwise_margin
            for row in exception_rows
        )
    )
    lattice_count_upper = (outer_rows + 1) * (outer_rows + 2) // 2
    aggregate_margin_lower = args.target_pointwise_margin - math.log2(lattice_count_upper)
    payload = {
        "schema": "single-random-constituent-shared-two-band-cover-v1",
        "status": "BINARY64_CONVEX_CELL_DIAGNOSTIC",
        "parameters": {
            "outer_rows": outer_rows,
            "minimum_occupation": args.minimum_occupation,
            "target_pointwise_margin_bits": args.target_pointwise_margin,
            "lattice_count_upper": lattice_count_upper,
            "aggregate_margin_lower_bits_if_complete": aggregate_margin_lower,
            "memory_bits": args.memory_bits,
            "distance_cutoff": distance_cutoff,
        },
        "spectrum_event": event,
        "bands": {"low": low, "central": central, "high": high},
        "merged_defect": {
            "value_probability": defect_probability,
            "log2_majorant": defect_log_majorant / LOG2,
        },
        "claim": {
            "complete": not unresolved and not pending,
            "complete_integer_lattice": (
                not pending and (not unresolved or lattice_exceptions_pass)
            ),
            "processed_cells": processed,
            "accepted_cells": len(accepted),
            "unresolved_cells": len(unresolved),
            "unresolved_integer_points": len(unresolved_points),
            "lattice_exceptions_pass": lattice_exceptions_pass,
            "worst_accepted_vertex_margin_bits": min(
                row["margin_bits"]
                for cell in accepted
                for row in cell["vertices_evaluated"]
            ) if accepted else None,
        },
        "accepted_cells": accepted,
        "unresolved_cells": unresolved,
        "lattice_exceptions": exception_rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Completeness relies on the stated fixed-witness convexity lemma.",
        ],
        "prior_receipt": None if args.prior is None else str(args.prior),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
