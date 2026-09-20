#!/usr/bin/env python3
"""Audit when elementary packing order appears across epoch depth.

The calculation keeps the 64 packet slots per epoch and varies the number of
uniformly averaged epochs.  It enumerates every profile and elementary edge
at a small fixed occupation.  Results are float64 diagnostics, not proofs.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import epoch_transfers
from analyze_riffle_ldpcsplitstate_shared_groups import (
    averaged_epoch_transfers,
    profiles_of_occupation,
)
from audit_riffle_ldpcsplitstate_region_packing import packed_neighbors


ROOT = Path("constructions/riffle_ldpcsplitstate_g4_t256_s64")
DEFAULT_ACTIVATION = ROOT / "receipts/zero_state_activation_table.json"
DEFAULT_OUTPUT = ROOT / "receipts/region_packing_epoch_depth.json"


def parse_ints(value: str) -> list[int]:
    return [int(item) for item in value.split(",") if item]


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


def placement_count(slot_count: int, counts: tuple[int, int, int, int]) -> int:
    groups = sum(counts)
    result = math.factorial(slot_count) // math.factorial(slot_count - groups)
    for count in counts:
        result //= math.factorial(count)
    return result


def region_transfer(
    averaged: list[np.ndarray],
    profile: tuple[int, int, int, int],
    epochs: int,
    slots_per_epoch: int,
) -> np.ndarray:
    zero = (0, 0, 0, 0)
    additions = list(itertools.product(*(range(count + 1) for count in profile)))
    additions = [counts for counts in additions if sum(counts) <= slots_per_epoch]
    factors = {
        counts: placement_count(slots_per_epoch, counts) for counts in additions
    }
    masses = {
        counts: sum((width + 1) * count for width, count in enumerate(counts))
        for counts in additions
    }
    coefficients: dict[tuple[int, int, int, int], np.ndarray] = {
        zero: np.eye(2)
    }
    for _ in range(epochs):
        following: dict[tuple[int, int, int, int], np.ndarray] = {}
        for used, prefix in coefficients.items():
            for added in additions:
                target = tuple(used[index] + added[index] for index in range(4))
                if any(target[index] > profile[index] for index in range(4)):
                    continue
                contribution = prefix @ averaged[masses[added]] * factors[added]
                following[target] = following.get(target, np.zeros((2, 2))) + contribution
        coefficients = following

    total_slots = epochs * slots_per_epoch
    normalization = placement_count(total_slots, profile)
    return coefficients[profile] / normalization


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    profiles = list(profiles_of_occupation(args.occupation))
    bounds = activation_bounds(args.activation, args.occupation)
    epoch = epoch_transfers(
        z=math.exp(-args.surprisal),
        distance=args.constituent_distance,
        moment_order=args.live_moment_order,
        activation_upper=bounds,
        maximum=args.occupation,
    )
    averaged = averaged_epoch_transfers(epoch, args.probability, args.occupation)
    rows = []
    for depth in args.epoch_depths:
        matrices = {
            profile: region_transfer(averaged, profile, depth, args.slots_per_epoch)
            for profile in profiles
        }
        edges = 0
        violations = 0
        maximum_ratio = 0.0
        witness = None
        for source, source_matrix in matrices.items():
            for target in packed_neighbors(source):
                if target not in matrices:
                    continue
                edges += 1
                target_matrix = matrices[target]
                ratio = float(np.max(source_matrix / target_matrix))
                if bool(
                    np.any(
                        source_matrix
                        > args.absolute_tolerance
                        + (1.0 + args.relative_tolerance) * target_matrix
                    )
                ):
                    violations += 1
                if ratio > maximum_ratio:
                    maximum_ratio = ratio
                    witness = {
                        "source_profile_n1_n2_n3_n4": list(source),
                        "target_profile_n1_n2_n3_n4": list(target),
                        "source_region_matrix": source_matrix.tolist(),
                        "target_region_matrix": target_matrix.tolist(),
                    }
        rows.append(
            {
                "epochs": depth,
                "total_packet_slots": depth * args.slots_per_epoch,
                "elementary_edges": edges,
                "entrywise_violations": violations,
                "maximum_entrywise_ratio": maximum_ratio,
                "maximum_ratio_witness": witness,
            }
        )
        print(
            f"epochs,{depth},edges,{edges},violations,{violations},"
            f"max_ratio,{maximum_ratio:.12g}",
            flush=True,
        )
    return {
        "schema": "riffle-ldpcsplitstate-region-packing-epoch-depth-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupation": args.occupation,
            "slots_per_epoch": args.slots_per_epoch,
            "probability": args.probability,
            "surprisal": args.surprisal,
            "epoch_depths": args.epoch_depths,
        },
        "rows": rows,
        "scope": (
            "Exact combinatorial coefficient averaging in float64 for all "
            "profiles at the recorded occupation. This is diagnostic evidence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occupation", type=int, default=8)
    parser.add_argument("--slots-per-epoch", type=int, default=64)
    parser.add_argument("--epoch-depths", type=parse_ints, default=parse_ints("1,2,4,8,16,32"))
    parser.add_argument("--probability", type=float, default=0.5)
    parser.add_argument("--surprisal", type=float, default=math.exp(-3.546242271905631))
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
