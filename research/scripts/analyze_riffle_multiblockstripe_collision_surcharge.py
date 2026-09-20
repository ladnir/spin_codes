#!/usr/bin/env python3
"""One shared-pair surcharge inside an iid-singleton background.

Fix one q-block group with active lane pattern {0,d}.  Add s active singleton
groups, each of which chooses an independent uniform lane in every region.
The script computes the exact floating-point region transfer by an exponential
generating function, averages the fresh cyclic shift of the pair, and compares
the resulting 256-region moment with a+2 fully iid lane choices.

This is a diagnostic for a possible composable collision surcharge.  It is not
an outward-rounded certificate and does not cover two or more defect groups.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_fieldcheckpoint_regular_envelope import regular_region_matrices
from analyze_riffle_multiblockstripe_fieldcheckpoint import (
    independent_lane_region_matrices,
    matrix_power_log_moment,
)
from analyze_riffle_striped_random_outer import LOG2


DEFAULT_OUTPUT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/one_pair_collision_surcharge_delta09.json"
)
DEFAULT_TILT_RECEIPT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/q32_independent_lane_low64_delta09.json"
)


def polynomial_matrix_product(
    left: np.ndarray, right: np.ndarray, degree: int
) -> np.ndarray:
    """Multiply two 2x2 matrix polynomials, truncated through degree."""
    result = np.zeros((degree + 1, 2, 2), dtype=np.float64)
    for row in range(2):
        for column in range(2):
            for inner in range(2):
                result[:, row, column] += np.convolve(
                    left[:, row, inner], right[:, inner, column]
                )[: degree + 1]
    return result


def defects_with_singletons_region(
    lane_matrices: list[np.ndarray],
    q: int,
    patterns: tuple[tuple[int, ...], ...],
    singletons: int,
) -> np.ndarray:
    """Exact EGF average for fixed defect patterns plus iid singletons."""
    inverse_factorials = np.asarray(
        [1.0 / math.factorial(count) for count in range(singletons + 1)]
    )
    factors = []
    for base in range(len(patterns) + 1):
        factor = np.zeros((singletons + 1, 2, 2), dtype=np.float64)
        for count in range(singletons + 1):
            factor[count] = lane_matrices[base + count] * inverse_factorials[count]
        factors.append(factor)

    averaged = np.zeros((2, 2), dtype=np.float64)
    for shifts in itertools.product(range(q), repeat=len(patterns)):
        occupations = [0] * q
        for pattern, shift in zip(patterns, shifts):
            for lane in pattern:
                occupations[(lane + shift) % q] += 1
        product = np.zeros((singletons + 1, 2, 2), dtype=np.float64)
        product[0] = np.eye(2)
        for lane in range(q):
            product = polynomial_matrix_product(
                product, factors[occupations[lane]], singletons
            )
        averaged += product[singletons]
    averaged /= q ** len(patterns)
    averaged *= math.factorial(singletons) / (q**singletons)
    return averaged


def pair_with_singletons_region(
    lane_matrices: list[np.ndarray], q: int, separation: int, singletons: int
) -> np.ndarray:
    return defects_with_singletons_region(
        lane_matrices, q, ((0, separation),), singletons
    )


def parse_int_list(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    q = args.q
    outer_bits = 256
    tilt_rows = json.loads(args.tilt_receipt.read_text(encoding="utf-8"))[
        "results"
    ][0]["independent_per_coordinate_lane_model"]["occupation_rows"]
    recorded_tilts = {
        int(row["active_regular_outer_blocks"]): float(row["best_log_surprisal"])
        for row in tilt_rows
    }
    rows = []
    for total in parse_int_list(args.total_occupations):
        if total < 2:
            raise ValueError("total occupation must be at least two")
        singletons = total - 2
        # Use the iid model's recorded low-occupation tilt schedule unless the
        # caller asks for a fixed tilt.
        log_surprisal = (
            args.log_surprisal
            if args.log_surprisal is not None
            else recorded_tilts[total]
        )
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        lanes = regular_region_matrices(
            args.state_bits,
            args.epoch_bits,
            1,
            total,
            z,
        )
        iid = independent_lane_region_matrices(
            q=q, lane_matrices=lanes, maximum_occupation=total
        )[total]
        iid_log_moment = matrix_power_log_moment(iid, outer_bits)

        separation_rows = []
        for separation in range(1, q // 2 + 1):
            defect = pair_with_singletons_region(
                lanes, q, separation, singletons
            )
            defect_log_moment = matrix_power_log_moment(defect, outer_bits)
            ratios = [
                float(defect[row, column] / iid[row, column])
                for row in range(2)
                for column in range(2)
                if iid[row, column] > 0.0
            ]
            separation_rows.append(
                {
                    "separation": separation,
                    "scalar_log2_surcharge": (
                        defect_log_moment - iid_log_moment
                    )
                    / LOG2,
                    "maximum_entrywise_log2_surcharge": math.log2(max(ratios)),
                }
            )
        worst_scalar = max(
            separation_rows, key=lambda row: float(row["scalar_log2_surcharge"])
        )
        worst_entrywise = max(
            separation_rows,
            key=lambda row: float(row["maximum_entrywise_log2_surcharge"]),
        )
        two_pair = None
        if total >= 4:
            two_pair_singletons = total - 4
            two_pair_region = defects_with_singletons_region(
                lanes, q, ((0, 1), (0, 1)), two_pair_singletons
            )
            two_pair_log_moment = matrix_power_log_moment(
                two_pair_region, outer_bits
            )
            two_pair_ratios = [
                float(two_pair_region[row, column] / iid[row, column])
                for row in range(2)
                for column in range(2)
                if iid[row, column] > 0.0
            ]
            two_pair = {
                "defect_groups": 2,
                "pattern_per_group": [0, 1],
                "singleton_groups": two_pair_singletons,
                "scalar_log2_surcharge": (
                    two_pair_log_moment - iid_log_moment
                )
                / LOG2,
                "scalar_log2_surcharge_per_defect": (
                    two_pair_log_moment - iid_log_moment
                )
                / (2 * LOG2),
                "maximum_entrywise_log2_surcharge": math.log2(
                    max(two_pair_ratios)
                ),
            }
        selected_patterns: list[tuple[str, tuple[int, ...]]] = []
        if total >= 3:
            selected_patterns.append(("triple_contiguous", (0, 1, 2)))
        if total >= 4:
            selected_patterns.append(("quad_contiguous", (0, 1, 2, 3)))
        if total >= 8:
            selected_patterns.extend(
                (
                    ("eight_contiguous", tuple(range(8))),
                    ("eight_evenly_spaced", tuple(range(0, q, q // 8))),
                )
            )
        if total >= 16:
            selected_patterns.extend(
                (
                    ("sixteen_contiguous", tuple(range(16))),
                    ("sixteen_alternating", tuple(range(0, q, 2))),
                )
            )
        if total >= q:
            selected_patterns.append(("full_group", tuple(range(q))))

        selected_group_rows = []
        for name, pattern in selected_patterns:
            background_singletons = total - len(pattern)
            region = defects_with_singletons_region(
                lanes, q, (pattern,), background_singletons
            )
            log_moment = matrix_power_log_moment(region, outer_bits)
            ratios = [
                float(region[row, column] / iid[row, column])
                for row in range(2)
                for column in range(2)
                if iid[row, column] > 0.0
            ]
            defects = len(pattern) - 1
            selected_group_rows.append(
                {
                    "name": name,
                    "group_weight": len(pattern),
                    "singleton_groups": background_singletons,
                    "defects": defects,
                    "scalar_log2_surcharge": (log_moment - iid_log_moment) / LOG2,
                    "scalar_log2_surcharge_per_defect": (
                        (log_moment - iid_log_moment) / (defects * LOG2)
                        if defects
                        else 0.0
                    ),
                    "maximum_entrywise_log2_surcharge": math.log2(max(ratios)),
                }
            )
        rows.append(
            {
                "total_active_blocks": total,
                "singleton_groups": singletons,
                "log_surprisal": log_surprisal,
                "worst_scalar": worst_scalar,
                "worst_entrywise": worst_entrywise,
                "two_adjacent_pair_defects": two_pair,
                "selected_single_group_patterns": selected_group_rows,
                "separation_rows": separation_rows,
            }
        )
        print(
            f"total,{total},tilt,{log_surprisal:.6f},"
            f"scalar_surcharge,{worst_scalar['scalar_log2_surcharge']:.9f},"
            f"scalar_sep,{worst_scalar['separation']},"
            f"entry_surcharge,{worst_entrywise['maximum_entrywise_log2_surcharge']:.9f},"
            f"two_pair_scalar,{None if two_pair is None else two_pair['scalar_log2_surcharge']}",
            flush=True,
        )
    return {
        "schema": "riffle-multiblockstripe-one-pair-surcharge-v1",
        "candidate": "Riffle MultiBlockStripeFresh-32 FieldCheckpoint",
        "parameters": {
            "q": q,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "outer_bits": outer_bits,
            "total_occupations": parse_int_list(args.total_occupations),
        },
        "rows": rows,
        "scope": (
            "Exact floating-point EGF for one same-group active pair plus iid "
            "singleton groups. Multiple defect groups and outward rounding "
            "are not covered."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--total-occupations", default="2,4,8,16,32,48,64")
    parser.add_argument("--log-surprisal", type=float)
    parser.add_argument("--tilt-receipt", type=Path, default=DEFAULT_TILT_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
