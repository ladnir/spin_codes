#!/usr/bin/env python3
"""Exact shared-group transfer diagnostic for MultiBlockStripeFresh-32.

For one q-block group, fix the active block lanes A subseteq Z_q.  A fresh
shift in each outer-coordinate region samples delta uniformly from Z_q and
uses A+delta.  The exact averaged region transfer is

    R_A(z) = q^-1 sum_delta product_r L_{1[r in A+delta]}(z),

where L_k is the exact FieldCheckpoint transfer for one epoch containing k
fair candidates.  Fresh shifts make the 256 region transfers independent, so
the complete inner moment is e_0^T R_A(z)^256 1.

The script compares this exact law with the iid-lane model previously used as
a proposed domination bound.  All reported proof comparisons are numerical;
they are not outward-rounded certificates.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_fieldcheckpoint_regular_envelope import (
    regular_region_matrices,
    spectrum_density_envelope_log,
)
from analyze_riffle_multiblockstripe_fieldcheckpoint import (
    independent_lane_region_matrices,
    matrix_power_log_moment,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import (
    LOG2,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/shared_group_small_profiles_delta09.json"
)


def fixed_group_profile_transfer(
    lane_matrices: list[np.ndarray],
    patterns: tuple[tuple[int, ...], ...],
    q: int,
) -> np.ndarray:
    """Return the exact region transfer for fixed active-lane patterns.

    Each pattern belongs to one q-block group.  The construction samples one
    independent cyclic shift per group and outer-coordinate region.  This
    routine enumerates those q^m shifts exactly.  It is intended for small m;
    the formula, rather than this exponential implementation, defines the
    general fixed-profile law.
    """
    averaged = np.zeros((2, 2), dtype=np.float64)
    shift_count = q ** len(patterns)
    for shifts in itertools.product(range(q), repeat=len(patterns)):
        occupations = [0] * q
        for pattern, shift in zip(patterns, shifts):
            for lane in pattern:
                occupations[(lane + shift) % q] += 1
        matrix = np.eye(2)
        for lane in range(q):
            matrix = matrix @ lane_matrices[occupations[lane]]
        averaged += matrix
    return averaged / shift_count


def rotated_pattern_transfer(
    lane_matrices: list[np.ndarray], pattern: tuple[int, ...], q: int
) -> np.ndarray:
    """Return R_A for one fixed binary cyclic pattern A."""
    return fixed_group_profile_transfer(lane_matrices, (pattern,), q)


def canonical_patterns(q: int, weight: int):
    """Yield a rotation-covering set by requiring lane zero to be active."""
    if weight < 1 or weight > q:
        return
    for tail in itertools.combinations(range(1, q), weight - 1):
        yield (0, *tail)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    q = args.q
    if args.region_bits // args.epoch_bits != q:
        raise ValueError("this diagnostic requires one q-lane slice per epoch")
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // q
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, _ = spectrum_density_envelope_log(args.outer_bits, dimension, spectrum)
    regular_block_log_mass = log_two_power_minus_one(dimension) + log_eta

    pattern_records: dict[tuple[int, ...], dict[str, object]] = {}
    iid_best = {weight: math.inf for weight in range(1, args.maximum_weight + 1)}
    iid_best_tilt = {weight: math.nan for weight in range(1, args.maximum_weight + 1)}
    maximum_entrywise_excess = {
        weight: (-math.inf, None, None, None)
        for weight in range(1, args.maximum_weight + 1)
    }

    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        lanes = regular_region_matrices(
            args.state_bits,
            args.epoch_bits,
            1,
            args.maximum_weight,
            z,
        )
        iid = independent_lane_region_matrices(
            q=q,
            lane_matrices=lanes,
            maximum_occupation=args.maximum_weight,
        )
        for weight in range(1, args.maximum_weight + 1):
            candidate = (
                matrix_power_log_moment(iid[weight], args.outer_bits)
                + distance * surprisal
            )
            if candidate < iid_best[weight]:
                iid_best[weight] = candidate
                iid_best_tilt[weight] = log_surprisal

            for pattern in canonical_patterns(q, weight):
                transfer = rotated_pattern_transfer(lanes, pattern, q)
                difference = transfer - iid[weight]
                flat_index = int(np.argmax(difference))
                entrywise_excess = float(difference.reshape(-1)[flat_index])
                if entrywise_excess > maximum_entrywise_excess[weight][0]:
                    maximum_entrywise_excess[weight] = (
                        entrywise_excess,
                        pattern,
                        log_surprisal,
                        divmod(flat_index, 2),
                    )

                candidate = (
                    matrix_power_log_moment(transfer, args.outer_bits)
                    + distance * surprisal
                )
                record = pattern_records.setdefault(
                    pattern,
                    {
                        "weight": weight,
                        "best_log_upper": math.inf,
                        "best_log_surprisal": math.nan,
                    },
                )
                if candidate < float(record["best_log_upper"]):
                    record["best_log_upper"] = candidate
                    record["best_log_surprisal"] = log_surprisal
        print(f"grid,{grid_index + 1},{grid_count},{log_surprisal:.6f}", flush=True)

    weight_rows = []
    for weight in range(1, args.maximum_weight + 1):
        records = [record | {"pattern": list(pattern)} for pattern, record in pattern_records.items() if int(record["weight"]) == weight]
        worst = max(records, key=lambda record: float(record["best_log_upper"]))
        iid_log = iid_best[weight]
        exact_subset_outer_log = math.log(groups) + weight * regular_block_log_mass
        orbit = len(
            {
                tuple(sorted((lane + shift) % q for lane in worst["pattern"]))
                for shift in range(q)
            }
        )
        orbit_outer_log = exact_subset_outer_log + math.log(orbit)
        excess, excess_pattern, excess_tilt, excess_entry = maximum_entrywise_excess[weight]
        weight_rows.append(
            {
                "weight": weight,
                "canonical_pattern_count": len(records),
                "worst_pattern": worst["pattern"],
                "worst_pattern_best_log_surprisal": worst["best_log_surprisal"],
                "worst_pattern_inner_log2_upper": min(0.0, float(worst["best_log_upper"])) / LOG2,
                "iid_inner_log2_upper": min(0.0, iid_log) / LOG2,
                "worst_minus_iid_inner_log2": (float(worst["best_log_upper"]) - iid_log) / LOG2,
                "worst_pattern_orbit_size": orbit,
                "worst_orbit_pointwise_log2_upper": (
                    orbit_outer_log + min(0.0, float(worst["best_log_upper"]))
                ) / LOG2,
                "worst_orbit_lambda_bits_lower_float": -(
                    orbit_outer_log + min(0.0, float(worst["best_log_upper"]))
                ) / LOG2,
                "maximum_entrywise_excess_over_iid": excess,
                "entrywise_counterexample_pattern": list(excess_pattern) if excess_pattern else None,
                "entrywise_counterexample_log_surprisal": excess_tilt,
                "entrywise_counterexample_entry": list(excess_entry) if excess_entry else None,
            }
        )

    return {
        "schema": "riffle-multiblockstripe-shared-group-small-profiles-v1",
        "candidate": "Riffle MultiBlockStripeFresh-32 FieldCheckpoint",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "q": q,
            "groups": groups,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "maximum_weight": args.maximum_weight,
        },
        "exact_law": (
            "R_A(z)=q^-1 sum_delta product_r L_{1[r in A+delta]}(z); "
            "complete moment=e0^T R_A(z)^256 1"
        ),
        "weight_rows": weight_rows,
        "scope": (
            "Exhaustive rotation-covering single-group patterns through the "
            "requested weight. Floating-point comparisons are diagnostics, "
            "not outward-rounded certificates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--region-bits", type=int, default=8192)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--maximum-weight", type=int, default=3)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-2.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in payload["weight_rows"]:
        print(
            f"weight,{row['weight']},pattern,{row['worst_pattern']},"
            f"worst_minus_iid_inner_log2,{row['worst_minus_iid_inner_log2']:.9f},"
            f"orbit_lambda,{row['worst_orbit_lambda_bits_lower_float']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
