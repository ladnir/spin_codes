#!/usr/bin/env python3
"""Low mixed regular/all-one envelope for SplitState-PreAddMul."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    log_choose,
    log_matrix_power_moments_batch,
    log_two_power_minus_one,
    logsumexp,
    modeled_even_floor_spectrum_logs,
    regular_region_log_matrices,
    spectrum_density_envelope_log,
    splitstate_impulse_matrices,
)


DEFAULT_ACTIVATION = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/preaddmul_mixed_low_total64.json"
)


def mixed_region_logs(
    active_support_logs: np.ndarray, regular: int, all_one: int
) -> np.ndarray:
    terms = []
    for regular_ones in range(regular + 1):
        terms.append(
            active_support_logs[all_one + regular_ones]
            + log_choose(regular, regular_ones)
            - regular * LOG2
        )
    stacked = np.stack(terms)
    result = np.empty((2, 2), dtype=np.float64)
    for row in range(2):
        for column in range(2):
            result[row, column] = logsumexp(stacked[:, row, column])
    return result


def load_nonactivation(path: Path, step_bits: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = np.ones(step_bits + 1, dtype=np.float64)
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > step_bits:
            continue
        key = "uniform_support_average_distinct_upper_bound"
        if key not in row:
            key = "uniform_256_support_average_distinct_upper_bound"
        result[weight] = float(row[key])
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    outer_dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // outer_dimension
    output_bits = 2 * args.message_bits
    target_distance = math.floor(args.relative_distance * output_bits)
    region_bits = output_bits // args.outer_bits
    epochs = region_bits // args.step_bits
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, outer_dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(outer_dimension) + log_eta
    nonactivation = load_nonactivation(args.activation, args.step_bits)
    pairs = [
        (regular, all_one)
        for regular in range(args.maximum_total_blocks + 1)
        for all_one in range(1, args.maximum_total_blocks - regular + 1)
    ]
    best = np.full(len(pairs), math.inf, dtype=np.float64)
    best_tilt = np.full(len(pairs), math.nan, dtype=np.float64)

    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        impulses = splitstate_impulse_matrices(
            z=z,
            step_bits=args.step_bits,
            state_bits=args.state_bits,
            constituent_distance=args.constituent_distance,
            live_moment_order=args.live_moment_order,
            live_model="support-averaged-preaddmul",
            nonactivation=nonactivation,
        )
        with np.errstate(divide="ignore"):
            active_epoch_logs = np.log(impulses)
        active_regions = regular_region_log_matrices(
            active_epoch_logs,
            args.step_bits,
            epochs,
            args.maximum_total_blocks,
        )
        matrices = np.stack(
            [
                mixed_region_logs(active_regions, regular, all_one)
                for regular, all_one in pairs
            ]
        )
        moments = log_matrix_power_moments_batch(matrices, args.outer_bits)
        candidates = moments + target_distance * surprisal
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    logs = []
    for index, (regular, all_one) in enumerate(pairs):
        outer_log = (
            log_choose(outer_blocks, regular)
            + log_choose(outer_blocks - regular, all_one)
            + regular * regular_log_mass
        )
        inner_log = min(0.0, float(best[index]))
        total = outer_log + inner_log
        logs.append(total)
        rows.append(
            {
                "regular_active_blocks": regular,
                "all_one_active_blocks": all_one,
                "best_log_surprisal": float(best_tilt[index]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": total / LOG2,
                "pointwise_margin_bits": -total / LOG2,
            }
        )
    aggregate = float(logsumexp(np.asarray(logs)))
    dominant = sorted(
        rows, key=lambda row: row["pointwise_log2_upper"], reverse=True
    )[:30]
    return {
        "schema": "riffle-splitstate-preaddmul-mixed-low-v1",
        "candidate": "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "target_distance": target_distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "maximum_total_blocks": args.maximum_total_blocks,
            "epochs_per_region": epochs,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Every mixed occupation with at least one all-one block and total "
            "occupation at most the recorded cap. The outer regular spectrum "
            "is modeled, arithmetic is nearest binary64, and no outward "
            "rounding is applied."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=32)
    parser.add_argument("--constituent-distance", type=int, default=20)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument("--maximum-total-blocks", type=int, default=64)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-8.0, -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0),
    )
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"aggregate_margin_bits,{payload['aggregate_margin_bits']:.9f}")
    for row in payload["dominant_rows"][:10]:
        print(
            f"dominant,{row['regular_active_blocks']},"
            f"{row['all_one_active_blocks']},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()

