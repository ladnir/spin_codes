#!/usr/bin/env python3
"""Audit elementary region-packing moves over a parameter grid.

This is finite floating-point evidence, not an outward-rounded proof.  It
checks every group profile at one small occupation, so every elementary edge
whose endpoints have that occupation is present.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import epoch_transfers
from analyze_riffle_ldpcsplitstate_shared_groups import (
    averaged_epoch_transfers,
    profile_region_transfer,
    profiles_of_occupation,
)
from audit_riffle_ldpcsplitstate_region_packing import packed_neighbors


ROOT = Path("constructions/riffle_ldpcsplitstate_g4_t256_s64")
DEFAULT_ACTIVATION = ROOT / "receipts/zero_state_activation_table.json"
DEFAULT_OUTPUT = ROOT / "receipts/region_packing_parameter_grid.json"

MOVE_NAMES = {
    (-2, 1, 0, 0): "1+1->2",
    (-1, -1, 1, 0): "1+2->3",
    (-1, 0, -1, 1): "1+3->4",
    (0, -2, 0, 1): "2+2->4",
    (1, -1, -1, 1): "2+3->1+4",
    (0, 1, -2, 1): "3+3->2+4",
}


def parse_floats(value: str) -> list[float]:
    return [float(item) for item in value.split(",") if item]


def activation_bounds(path: Path, maximum: int) -> list[float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = [0.0] * (maximum + 1)
    result[0] = 1.0
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > maximum:
            break
        result[weight] = float(row["maximum_distinct_conditioned_upper_bound"])
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    profiles = list(profiles_of_occupation(args.occupation))
    bounds = activation_bounds(args.activation, args.occupation)
    rows: list[dict[str, object]] = []
    total_edges = 0
    total_violations = 0
    largest_excess = -float("inf")
    largest_ratio = 0.0
    largest_witness = None
    by_move = {
        name: {"edges": 0, "entrywise_violations": 0, "maximum_entrywise_ratio": 0.0}
        for name in MOVE_NAMES.values()
    }

    for probability in args.probabilities:
        for surprisal in args.surprisals:
            z = math.exp(-surprisal)
            epoch = epoch_transfers(
                z=z,
                distance=args.constituent_distance,
                moment_order=args.live_moment_order,
                activation_upper=bounds,
                maximum=args.occupation,
            )
            averaged = averaged_epoch_transfers(
                epoch, probability, args.occupation
            )
            matrices = {
                profile: profile_region_transfer(averaged, profile)
                for profile in profiles
            }
            point_edges = 0
            point_violations = 0
            point_largest_ratio = 0.0
            point_largest_witness = None
            for source, source_matrix in matrices.items():
                for target in packed_neighbors(source):
                    if target not in matrices:
                        continue
                    point_edges += 1
                    delta = tuple(target[i] - source[i] for i in range(4))
                    move_row = by_move[MOVE_NAMES[delta]]
                    move_row["edges"] += 1
                    target_matrix = matrices[target]
                    excess = float(np.max(source_matrix - target_matrix))
                    ratio = float(np.max(source_matrix / target_matrix))
                    violation_mask = source_matrix > (
                        args.absolute_tolerance
                        + (1.0 + args.relative_tolerance) * target_matrix
                    )
                    if bool(np.any(violation_mask)):
                        point_violations += 1
                        move_row["entrywise_violations"] += 1
                    if ratio > move_row["maximum_entrywise_ratio"]:
                        move_row["maximum_entrywise_ratio"] = ratio
                    if ratio > point_largest_ratio:
                        point_largest_ratio = ratio
                        point_largest_witness = {
                            "source_profile_n1_n2_n3_n4": list(source),
                            "target_profile_n1_n2_n3_n4": list(target),
                            "source_region_matrix": source_matrix.tolist(),
                            "target_region_matrix": target_matrix.tolist(),
                        }
                    if excess > largest_excess:
                        largest_excess = excess
                    if ratio > largest_ratio:
                        largest_ratio = ratio
                        largest_witness = {
                            "probability": probability,
                            "surprisal": surprisal,
                            "source_profile_n1_n2_n3_n4": list(source),
                            "target_profile_n1_n2_n3_n4": list(target),
                            "maximum_entrywise_ratio": ratio,
                            "maximum_absolute_excess": excess,
                        }
            total_edges += point_edges
            total_violations += point_violations
            rows.append(
                {
                    "probability": probability,
                    "surprisal": surprisal,
                    "elementary_edges": point_edges,
                    "entrywise_violations": point_violations,
                    "maximum_entrywise_ratio": point_largest_ratio,
                    "maximum_ratio_witness": point_largest_witness,
                }
            )
            print(
                f"p,{probability:.6g},surprisal,{surprisal:.6g},"
                f"edges,{point_edges},violations,{point_violations},"
                f"max_ratio,{point_largest_ratio:.12g}",
                flush=True,
            )

    return {
        "schema": "riffle-ldpcsplitstate-region-packing-grid-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupation": args.occupation,
            "probabilities": args.probabilities,
            "surprisals": args.surprisals,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
        },
        "summary": {
            "parameter_points": len(rows),
            "total_elementary_edges": total_edges,
            "entrywise_violations": total_violations,
            "largest_absolute_excess": largest_excess,
            "largest_entrywise_ratio": largest_ratio,
            "largest_ratio_witness": largest_witness,
            "by_elementary_move": by_move,
        },
        "rows": rows,
        "scope": (
            "Exact combinatorial region averaging evaluated in float64 on a "
            "finite parameter grid. This is diagnostic evidence, not an "
            "outward-rounded proof or a statement outside the recorded grid."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occupation", type=int, default=8)
    parser.add_argument(
        "--probabilities", type=parse_floats, default=parse_floats("0.1,0.3,0.5,0.7,0.9")
    )
    parser.add_argument(
        "--surprisals", type=parse_floats, default=parse_floats("0.01,0.03,0.1,0.3,1")
    )
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-300)
    parser.add_argument("--relative-tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
