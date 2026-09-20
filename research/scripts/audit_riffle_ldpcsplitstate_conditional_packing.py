#!/usr/bin/env python3
"""Test packing after conditioning on an arbitrary region background.

For a fixed placement of every background packet, this script averages only
the one or two packet groups changed by an elementary packing move over the
remaining packet slots.  Entrywise domination here would be a sufficient
local lemma for the full region-packing induction.  The randomized audit is
diagnostic and is not a proof.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import epoch_transfers
from analyze_riffle_ldpcsplitstate_shared_groups import averaged_epoch_transfers


ROOT = Path("constructions/riffle_ldpcsplitstate_g4_t256_s64")
DEFAULT_ACTIVATION = ROOT / "receipts/zero_state_activation_table.json"
DEFAULT_OUTPUT = ROOT / "receipts/conditional_region_packing_audit.json"
MOVES = [
    ("1+1->2", (1, 1), (2,)),
    ("1+2->3", (1, 2), (3,)),
    ("1+3->4", (1, 3), (4,)),
    ("2+2->4", (2, 2), (4,)),
    ("2+3->1+4", (2, 3), (1, 4)),
    ("3+3->2+4", (3, 3), (2, 4)),
]


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


def product_with_additions(
    averaged: list[np.ndarray], background_mass: list[int], additions: dict[int, int]
) -> np.ndarray:
    result = np.eye(2)
    for epoch, mass in enumerate(background_mass):
        result = result @ averaged[mass + additions.get(epoch, 0)]
    return result


def insertion_kernel(
    averaged: list[np.ndarray],
    background_mass: list[int],
    free_slots: list[int],
    widths: tuple[int, ...],
) -> np.ndarray:
    total_free = sum(free_slots)
    result = np.zeros((2, 2))
    if len(widths) == 1:
        for epoch, capacity in enumerate(free_slots):
            if capacity:
                result += capacity * product_with_additions(
                    averaged, background_mass, {epoch: widths[0]}
                )
        return result / total_free

    first, second = widths
    for left, left_capacity in enumerate(free_slots):
        if not left_capacity:
            continue
        for right, right_capacity in enumerate(free_slots):
            placements = left_capacity * (
                right_capacity - (1 if left == right else 0)
            )
            if placements <= 0:
                continue
            additions = {left: first}
            additions[right] = additions.get(right, 0) + second
            result += placements * product_with_additions(
                averaged, background_mass, additions
            )
    return result / (total_free * (total_free - 1))


def sampled_background(rng: random.Random, total_mass: int) -> tuple[list[int], list[int]]:
    masses = [0] * 32
    used = [0] * 32
    remaining = total_mass
    while remaining:
        width = rng.randint(1, min(4, remaining))
        candidates = [epoch for epoch in range(32) if used[epoch] < 64]
        epoch = rng.choice(candidates)
        masses[epoch] += width
        used[epoch] += 1
        remaining -= width
    return masses, [64 - count for count in used]


def canonical_backgrounds(total_mass: int) -> list[tuple[list[int], list[int]]]:
    if total_mass == 0:
        return [([0] * 32, [64] * 32)]
    result = []
    for occupied_epochs in (1, 2, 4, 8, 16, 32):
        masses = [0] * 32
        used = [0] * 32
        remaining = total_mass
        index = 0
        while remaining:
            width = min(4, remaining)
            epoch = index % occupied_epochs
            masses[epoch] += width
            used[epoch] += 1
            remaining -= width
            index += 1
        result.append((masses, [64 - count for count in used]))
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    maximum_epoch_mass = args.occupation
    bounds = activation_bounds(args.activation, maximum_epoch_mass)
    z = math.exp(-args.surprisal)
    epoch = epoch_transfers(
        z=z,
        distance=args.constituent_distance,
        moment_order=args.live_moment_order,
        activation_upper=bounds,
        maximum=maximum_epoch_mass,
    )
    averaged = averaged_epoch_transfers(
        epoch, args.probability, maximum_epoch_mass
    )
    rng = random.Random(args.seed)
    by_move = {
        name: {"backgrounds": 0, "entrywise_violations": 0, "maximum_entrywise_ratio": 0.0}
        for name, _, _ in MOVES
    }
    maximum_ratio = 0.0
    maximum_witness = None
    total_backgrounds = 0
    for name, source_widths, target_widths in MOVES:
        background_mass = args.occupation - sum(source_widths)
        backgrounds = canonical_backgrounds(background_mass)
        backgrounds.extend(
            sampled_background(rng, background_mass)
            for _ in range(args.random_trials)
        )
        total_backgrounds += len(backgrounds)
        for background_index, (masses, free) in enumerate(backgrounds):
            source = insertion_kernel(averaged, masses, free, source_widths)
            target = insertion_kernel(averaged, masses, free, target_widths)
            ratio = float(np.max(source / target))
            violation = bool(
                np.any(
                    source
                    > args.absolute_tolerance
                    + (1.0 + args.relative_tolerance) * target
                )
            )
            row = by_move[name]
            row["backgrounds"] += 1
            row["entrywise_violations"] += int(violation)
            row["maximum_entrywise_ratio"] = max(
                row["maximum_entrywise_ratio"], ratio
            )
            if ratio > maximum_ratio:
                maximum_ratio = ratio
                maximum_witness = {
                    "background_index": background_index,
                    "background_total_mass": background_mass,
                    "move": name,
                    "background_mass_by_epoch": masses,
                    "free_slots_by_epoch": free,
                    "source_region_matrix": source.tolist(),
                    "target_region_matrix": target.tolist(),
                    "maximum_entrywise_ratio": ratio,
                }

    violations = sum(row["entrywise_violations"] for row in by_move.values())
    return {
        "schema": "riffle-ldpcsplitstate-conditional-region-packing-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "probability": args.probability,
            "surprisal": args.surprisal,
            "occupation": args.occupation,
            "random_trials": args.random_trials,
            "canonical_backgrounds_per_move": len(
                canonical_backgrounds(args.occupation - sum(MOVES[0][1]))
            ),
            "seed": args.seed,
        },
        "summary": {
            "backgrounds_across_moves": total_backgrounds,
            "move_background_pairs": total_backgrounds,
            "entrywise_violations": violations,
            "maximum_entrywise_ratio": maximum_ratio,
            "maximum_ratio_witness": maximum_witness,
            "by_move": by_move,
        },
        "scope": (
            "Floating-point audit of a sufficient conditional insertion "
            "inequality on canonical and pseudorandom backgrounds. It is not "
            "a proof and does not enumerate every feasible background."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probability", type=float, default=0.5)
    parser.add_argument("--surprisal", type=float, default=math.exp(-3.546242271905631))
    parser.add_argument("--occupation", type=int, default=128)
    parser.add_argument("--random-trials", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260829)
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
    summary = payload["summary"]
    print(
        f"backgrounds,{summary['backgrounds_across_moves']},"
        f"move_background_pairs,{summary['move_background_pairs']},"
        f"violations,{summary['entrywise_violations']},"
        f"max_ratio,{summary['maximum_entrywise_ratio']:.12g}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
