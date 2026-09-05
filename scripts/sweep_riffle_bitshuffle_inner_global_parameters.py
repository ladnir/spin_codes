#!/usr/bin/env python3
"""Compare FieldCheckpoint and SplitState under one bit-transpose model.

This is a one-active outer-block landscape, not a complete distance proof.
Both inner families use the same modeled outer spectrum, independent outer
coordinate permutations, and independent uniform bit permutation per region.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import region_matrices as field_region_matrices
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
    weight_conditioned_log_moments,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/"
    "receipts/inner_global_parameter_landscape.json"
)
LOG2 = math.log(2.0)
OUTER_POINTS = ((128, 22), (256, 38), (512, 56))
SPLIT_POINTS = ((128, 32, 20), (256, 64, 40), (512, 128, 80))
FIELD_EPOCHS = (256, 512, 1024)


def split_epoch_matrices(
    *, z: float, step_bits: int, state_bits: int, constituent_distance: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return the zero-input and one-impulse SplitState upper matrices."""
    denominator = math.ldexp(1.0, state_bits) - 1.0
    mean_weight = step_bits * math.ldexp(1.0, state_bits - 1) / denominator
    low_fraction = (step_bits - mean_weight) / (
        step_bits - constituent_distance
    )
    live_moment = (
        low_fraction * z**constituent_distance
        + (1.0 - low_fraction) * z**step_bits
    )
    zero = np.asarray(((1.0, 0.0), (0.0, live_moment)), dtype=np.float64)

    sparse_floor = z ** max(0, constituent_distance - 1)
    coset_bound = live_moment + (1.0 - z) / denominator
    impulse = np.asarray(
        (
            (0.0, z),
            (sparse_floor / denominator, min(sparse_floor, coset_bound)),
        ),
        dtype=np.float64,
    )
    return zero, impulse


def marked_region(
    zero_epoch: np.ndarray, impulse_epoch: np.ndarray, epochs: int
) -> tuple[np.ndarray, np.ndarray]:
    """Average one impulse uniformly over the epochs of one region."""
    powers = [np.eye(2)]
    for _ in range(epochs):
        powers.append(powers[-1] @ zero_epoch)
    inactive = powers[epochs]
    active = np.zeros((2, 2), dtype=np.float64)
    for epoch in range(epochs):
        active += powers[epoch] @ impulse_epoch @ powers[epochs - 1 - epoch]
    active /= epochs
    return inactive, active


def evaluate_one_active(
    *,
    outer_bits: int,
    minimum_distance: int,
    distance: int,
    outer_blocks: int,
    matrix_factory,
    grid: np.ndarray,
) -> dict[str, object]:
    spectrum = modeled_even_floor_spectrum_logs(outer_bits, minimum_distance)
    best = np.full(outer_bits + 1, math.inf)
    best_tilt = np.full(outer_bits + 1, math.nan)
    for log_surprisal in grid:
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        inactive, active = matrix_factory(z)
        moments = weight_conditioned_log_moments(
            outer_bits=outer_bits, zero_row=inactive, active_row=active
        )
        values = moments + distance * surprisal
        improved = values < best
        best[improved] = values[improved]
        best_tilt[improved] = log_surprisal

    contributions = []
    rows = []
    for weight in range(1, outer_bits + 1):
        if not math.isfinite(float(spectrum[weight])):
            continue
        inner = min(0.0, float(best[weight]))
        total = math.log(outer_blocks) + float(spectrum[weight]) + inner
        contributions.append(total)
        rows.append((total, weight, inner, float(best_tilt[weight])))
    aggregate = float(logsumexp(np.asarray(contributions)))
    dominant = max(rows)
    return {
        "one_active_margin_bits": -aggregate / LOG2,
        "dominant_outer_weight": dominant[1],
        "dominant_inner_log2": dominant[2] / LOG2,
        "dominant_log_surprisal": dominant[3],
    }


def split_logical_xor_proxy(step_bits: int, state_bits: int) -> float:
    """Scale the recorded t=256,s=64 nested-pair logical accounting."""
    if step_bits != 4 * state_bits:
        raise ValueError("the calibrated SplitState proxy requires t=4s")
    structural = 4.2421875 + 2.25 + 1.0
    # The recorded coefficient-specific 64-bit field circuit uses 799.375
    # XORs. Scale the random-linear circuit proxy as s^2/log2(s).
    field_map = 799.375 * (state_bits / 64.0) ** 2
    field_map *= 6.0 / math.log2(state_bits)
    return structural + field_map / step_bits


def field_logical_xor_proxy(epoch_bits: int) -> float:
    """Use the recorded 64-bit coefficient-specific field circuit count."""
    return 1.0 + 799.375 / epoch_bits


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    grid = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    rows = []
    for outer_bits, outer_distance in OUTER_POINTS:
        outer_dimension = outer_bits // 2
        if args.message_bits % outer_dimension:
            continue
        outer_blocks = args.message_bits // outer_dimension
        region_bits = output_bits // outer_bits

        for step_bits, state_bits, constituent_distance in SPLIT_POINTS:
            if region_bits % step_bits:
                continue
            epochs = region_bits // step_bits

            def split_factory(z: float, *, _epochs=epochs, _t=step_bits, _s=state_bits, _d=constituent_distance):
                zero, impulse = split_epoch_matrices(
                    z=z,
                    step_bits=_t,
                    state_bits=_s,
                    constituent_distance=_d,
                )
                return marked_region(zero, impulse, _epochs)

            result = evaluate_one_active(
                outer_bits=outer_bits,
                minimum_distance=outer_distance,
                distance=distance,
                outer_blocks=outer_blocks,
                matrix_factory=split_factory,
                grid=grid,
            )
            rows.append(
                {
                    "inner_family": "SplitState",
                    "outer_bits": outer_bits,
                    "outer_minimum_distance": outer_distance,
                    "outer_blocks": outer_blocks,
                    "region_bits": region_bits,
                    "step_bits": step_bits,
                    "state_bits": state_bits,
                    "constituent_distance": constituent_distance,
                    "epochs_per_region": epochs,
                    "logical_inner_xor_per_output_proxy": split_logical_xor_proxy(
                        step_bits, state_bits
                    ),
                    **result,
                }
            )

        for epoch_bits in FIELD_EPOCHS:
            state_bits = 64
            if region_bits % epoch_bits or epoch_bits % state_bits:
                continue
            epochs = region_bits // epoch_bits

            def field_factory(z: float, *, _epochs=epochs, _epoch_bits=epoch_bits):
                return field_region_matrices(64, _epoch_bits, _epochs, z)

            result = evaluate_one_active(
                outer_bits=outer_bits,
                minimum_distance=outer_distance,
                distance=distance,
                outer_blocks=outer_blocks,
                matrix_factory=field_factory,
                grid=grid,
            )
            rows.append(
                {
                    "inner_family": "FieldCheckpoint",
                    "outer_bits": outer_bits,
                    "outer_minimum_distance": outer_distance,
                    "outer_blocks": outer_blocks,
                    "region_bits": region_bits,
                    "epoch_bits": epoch_bits,
                    "state_bits": state_bits,
                    "epochs_per_region": epochs,
                    "logical_inner_xor_per_output_proxy": field_logical_xor_proxy(
                        epoch_bits
                    ),
                    **result,
                }
            )
        print(f"outer_complete,{outer_bits}", flush=True)

    feasible = [row for row in rows if row["one_active_margin_bits"] >= args.target_margin]
    proof_best = max(rows, key=lambda row: row["one_active_margin_bits"])
    cost_best = min(
        feasible,
        key=lambda row: row["logical_inner_xor_per_output_proxy"],
    ) if feasible else None
    return {
        "schema": "riffle-bitshuffle-inner-global-parameter-landscape-v1",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "relative_distance": args.relative_distance,
            "distance": distance,
            "target_margin_bits": args.target_margin,
            "outer_points_B_d": [list(point) for point in OUTER_POINTS],
            "split_points_t_s_d": [list(point) for point in SPLIT_POINTS],
            "field_epoch_bits": list(FIELD_EPOCHS),
        },
        "proof_best": proof_best,
        "lowest_inner_xor_proxy_above_target": cost_best,
        "rows": rows,
        "scope": (
            "One-active outer-block floating-point landscape. Outer spectra "
            "are modeled. SplitState distances at t=128 and t=512 extrapolate "
            "the recorded t=256,s=64,d=40 relative distance and are assumptions. "
            "The XOR proxy excludes memory traffic, routing, outer encoding, "
            "circuit sharing, and instruction-level effects."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--target-margin", type=float, default=40.0)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=0.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for row in payload["rows"]:
        interval = row.get("step_bits", row.get("epoch_bits"))
        print(
            f"row,{row['inner_family']},B,{row['outer_bits']},"
            f"interval,{interval},s,{row['state_bits']},"
            f"margin,{row['one_active_margin_bits']:.6f},"
            f"xor_proxy,{row['logical_inner_xor_per_output_proxy']:.6f}",
            flush=True,
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
