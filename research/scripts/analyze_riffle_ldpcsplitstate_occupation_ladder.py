#!/usr/bin/env python3
"""Occupation ladder through a Bernoulli-conditioned two-state transfer."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from analyze_riffle_ldpcsplitstate_sparse_prefix import (
    log_choose,
    modeled_log2_multiplicity,
)


DEFAULT_ACTIVATION = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/occupation_ladder_4_8_16.json"
)


def deterministic_live_moment(
    z: float,
    distance: int,
    order: int,
    step_bits: int = 256,
    state_bits: int = 64,
) -> float:
    denominator = (1 << state_bits) - 1
    if order > 1:
        weights = np.arange(distance, step_bits + 1, dtype=np.float64)
        equality_rows = [np.ones_like(weights)]
        equality_values = [1.0]
        falling = np.ones_like(weights)
        normalization = 1.0
        for degree in range(1, order + 1):
            falling *= weights - (degree - 1)
            normalization *= step_bits - (degree - 1)
            equality_rows.append(falling / normalization)
            equality_values.append((1 << (state_bits - degree)) / denominator)
        # At the useful Chernoff tilts z^distance can be much smaller than
        # HiGHS' absolute feasibility/optimality tolerances.  Normalize the
        # objective so its largest coefficient is one; otherwise the solver
        # can return visibly nonmonotone and spuriously tiny moment bounds.
        objective_scale = z**distance
        normalized_objective = (z**weights) / objective_scale
        result = linprog(
            c=-normalized_objective,
            A_eq=np.asarray(equality_rows),
            b_eq=np.asarray(equality_values),
            bounds=(0.0, None),
            method="highs",
        )
        if not result.success:
            raise RuntimeError(f"moment LP failed: {result.message}")
        return -float(result.fun) * objective_scale
    mean_weight = step_bits * (1 << (state_bits - 1)) / denominator
    low_fraction = (step_bits - mean_weight) / (step_bits - distance)
    return low_fraction * z**distance + (1.0 - low_fraction) * z**step_bits


def epoch_transfers(
    *,
    z: float,
    distance: int,
    moment_order: int,
    activation_upper: list[float],
    maximum: int,
) -> list[np.ndarray]:
    denominator = (1 << 64) - 1
    live_moment = deterministic_live_moment(z, distance, moment_order)
    transfers = []
    for weight in range(maximum + 1):
        matrix = np.zeros((2, 2), dtype=np.float64)
        if weight == 0:
            matrix[0, 0] = 1.0
            matrix[1, 1] = live_moment
            transfers.append(matrix)
            continue
        output_factor = z**weight
        nonactivation = activation_upper[weight]
        if nonactivation == 0.0:
            matrix[0, 1] = output_factor
        else:
            # Entrywise upper bound: activation and nonactivation are each
            # bounded separately.  Their sum may exceed the exact row sum.
            matrix[0, 0] = nonactivation * output_factor
            matrix[0, 1] = output_factor
        sparse_floor = z ** max(0, distance - weight)
        coset_bound = live_moment + (1.0 - output_factor) / denominator
        matrix[1, 0] = sparse_floor / denominator
        matrix[1, 1] = min(sparse_floor, coset_bound)
        transfers.append(matrix)
    return transfers


def singleton_region_transfers(
    epoch: list[np.ndarray], maximum_packets: int
) -> list[np.ndarray]:
    coefficients = [np.zeros((2, 2)) for _ in range(maximum_packets + 1)]
    coefficients[0] = np.eye(2)
    for epoch_index in range(32):
        following = [np.zeros((2, 2)) for _ in range(maximum_packets + 1)]
        used_maximum = min(maximum_packets, 64 * epoch_index)
        for used in range(used_maximum + 1):
            if not np.any(coefficients[used]):
                continue
            for added in range(min(64, maximum_packets - used) + 1):
                following[used + added] += (
                    coefficients[used]
                    @ epoch[added]
                    * math.comb(64, added)
                )
        coefficients = following
    return [
        coefficients[packets] / math.comb(2048, packets)
        for packets in range(maximum_packets + 1)
    ]


def bernoulli_region_transfer(
    region: list[np.ndarray], active_blocks: int, probability: float
) -> np.ndarray:
    result = np.zeros((2, 2))
    for active in range(active_blocks + 1):
        coefficient = (
            math.comb(active_blocks, active)
            * probability**active
            * (1.0 - probability) ** (active_blocks - active)
        )
        result += coefficient * region[active]
    return result


def regular_outer_density_log2(probability: float, minimum_distance: int = 38) -> float:
    weights = list(range(minimum_distance, 256 - minimum_distance + 1, 2))
    normalization = (
        math.log2((1 << 128) - 2)
        - math.log2(sum(math.comb(256, weight) for weight in weights))
    )
    return max(
        normalization
        - weight * math.log2(probability)
        - (256 - weight) * math.log2(1.0 - probability)
        for weight in weights
    )


def matrix_power_moment(matrix: np.ndarray, regions: int) -> float:
    power = np.eye(2)
    scale_log = 0.0
    for _ in range(regions):
        power = power @ matrix
        maximum = float(np.max(power))
        power /= maximum
        scale_log += math.log(maximum)
    terminal = float(np.sum(power[0]))
    return scale_log + math.log(terminal)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    activation = json.loads(args.activation.read_text(encoding="utf-8"))
    activation_upper = [0.0] * (max(args.occupations) + 1)
    activation_upper[0] = 1.0
    for row in activation["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight >= len(activation_upper):
            break
        activation_upper[weight] = float(
            row["maximum_distinct_conditioned_upper_bound"]
        )

    nominal_probability = args.outer_weight / 256
    local_multiplicity_log2 = modeled_log2_multiplicity(args.outer_weight)
    target_distance = math.floor(args.relative_distance * (1 << 21))
    best = {
        occupation: None for occupation in args.occupations
    }
    grid_rows = []
    for log_surprisal in np.linspace(-10.0, 1.0, 133):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        epoch = epoch_transfers(
            z=z,
            distance=args.constituent_distance,
            moment_order=args.live_moment_order,
            activation_upper=activation_upper,
            maximum=max(args.occupations),
        )
        region = singleton_region_transfers(epoch, max(args.occupations))
        for occupation in args.occupations:
            for placement_logit in np.linspace(-8.0, 5.0, 105):
                probability = 1.0 / (1.0 + math.exp(-float(placement_logit)))
                average_region = bernoulli_region_transfer(
                    region, occupation, probability
                )
                log_moment = matrix_power_moment(average_region, 256)
                if args.outer_mode == "fixed-w38":
                    placement_conditioning_bits = -occupation * (
                        log_choose(256, args.outer_weight) / math.log(2.0)
                        + args.outer_weight * math.log2(probability)
                        + (256 - args.outer_weight) * math.log2(1.0 - probability)
                    )
                    local_outer_bits = occupation * local_multiplicity_log2
                else:
                    placement_conditioning_bits = (
                        occupation * regular_outer_density_log2(probability)
                    )
                    local_outer_bits = 0.0
                outer_set_log2 = (
                    log_choose(2048, occupation) / math.log(2.0)
                    + 2 * occupation
                )
                multiplicity_log2 = outer_set_log2 + local_outer_bits
                log2_bound = (
                    log_moment / math.log(2.0)
                    + placement_conditioning_bits
                    + multiplicity_log2
                    + target_distance * surprisal / math.log(2.0)
                )
                row = {
                    "occupation": occupation,
                    "log_surprisal": float(log_surprisal),
                    "z": z,
                    "placement_logit": float(placement_logit),
                    "placement_probability": probability,
                    "log2_inner_moment": log_moment / math.log(2.0),
                    "placement_conditioning_bits": placement_conditioning_bits,
                    "modeled_multiplicity_bits": multiplicity_log2,
                    "aggregate_log2_bound": log2_bound,
                    "aggregate_margin_bits": -log2_bound,
                    "average_region_matrix": average_region.tolist(),
                }
                if (
                    best[occupation] is None
                    or log2_bound < best[occupation]["aggregate_log2_bound"]
                ):
                    best[occupation] = row
            grid_rows.append(best[occupation])

    return {
        "schema": "riffle-ldpcsplitstate-occupation-ladder-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupations": list(args.occupations),
            "outer_weight": args.outer_weight,
            "outer_mode": args.outer_mode,
            "nominal_outer_region_probability": nominal_probability,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
        },
        "best": {str(key): value for key, value in best.items()},
        "grid": grid_rows,
        "scope": (
            "Floating-point diagnostic for active outer blocks in distinct "
            "four-block groups. Fixed-w38 mode uses one coefficient bound per "
            "block. Regular-envelope mode uses a pointwise Bernoulli-density "
            "envelope for every modeled even weight from 38 through 218. "
            "Epoch inputs "
            "are singleton packets before packet shuffling. Activation uses "
            "the worst split-weight upper table, and live moments use the "
            "requested deterministic factorial-moment bound. Configurations "
            "with two active blocks in one packet group are not included."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--occupations", type=int, nargs="+", default=(1, 2, 4, 8, 16))
    parser.add_argument("--outer-weight", type=int, default=38)
    parser.add_argument(
        "--outer-mode",
        choices=("fixed-w38", "regular-envelope"),
        default="fixed-w38",
    )
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument(
        "--live-moment-order", type=int, choices=(1, 2, 3), default=1
    )
    parser.add_argument("--relative-distance", type=float, default=0.09)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for occupation in args.occupations:
        row = payload["best"][str(occupation)]
        print(
            f"occupation,{occupation},margin_bits,{row['aggregate_margin_bits']:.6f},"
            f"log_surprisal,{row['log_surprisal']:.6f},"
            f"placement_probability,{row['placement_probability']:.6f},"
            f"conditioning_bits,{row['placement_conditioning_bits']:.6f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
