#!/usr/bin/env python3
"""Exact modeled one-active outer block for SplitState-PreAddMul."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    load_outer_spectrum_logs,
    random_linear_extension_spectrum_logs,
    logsumexp,
    log_choose,
    modeled_even_floor_spectrum_logs,
    load_nonzero_spectrum,
    splitstate_impulse_matrices,
    support_averaged_exact_spectrum_moments,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    weight_conditioned_log_moments,
)


DEFAULT_ACTIVATION = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/preaddmul_support_averaged_moment3_one_active.json"
)


def marked_region(
    zero_epoch: np.ndarray, impulse_epoch: np.ndarray, epochs: int
) -> tuple[np.ndarray, np.ndarray]:
    powers = [np.eye(zero_epoch.shape[0])]
    for _ in range(epochs):
        powers.append(powers[-1] @ zero_epoch)
    inactive = powers[epochs]
    active = np.zeros_like(zero_epoch, dtype=np.float64)
    for epoch in range(epochs):
        active += powers[epoch] @ impulse_epoch @ powers[epochs - 1 - epoch]
    active /= epochs
    return inactive, active


def refreshed_three_state_impulse_matrices(
    *, z: float, step_bits: int, state_bits: int, live_spectrum: np.ndarray
) -> np.ndarray:
    """Exact one-impulse envelope with zero/live/punctured-live states.

    State 1 is uniform on F*, while state 2 is bounded by a uniform F*
    distribution with one nonzero value omitted.  A zero-input epoch maps
    either live class back to uniform F*, because multiplication by the fresh
    nonzero scalar is transitive on F*.  For a weight-one input, B(X) is
    nonzero, so conditioning the next state on being live omits exactly zero
    and B(X).
    """
    denominator = math.ldexp(1.0, state_bits) - 1.0
    punctured_factor = denominator / (denominator - 1.0)
    moments = support_averaged_exact_spectrum_moments(
        z=z,
        step_bits=step_bits,
        spectrum=live_spectrum,
    )
    matrices = np.zeros((2, 3, 3), dtype=np.float64)

    # Zero input.  The current state determines the emitted A(Q); the fresh
    # scalar then makes every nonzero next state exactly uniform on F*.
    matrices[0, 0, 0] = 1.0
    matrices[0, 1, 1] = float(moments[0])
    matrices[0, 2, 1] = punctured_factor * float(moments[0])

    # Weight-one input.  Distinct nonzero columns of B make B(X) nonzero.
    # Termination has probability 1/|F*| independently of the emitted word.
    output_factor = z
    matrices[1, 0, 2] = output_factor
    uniform_moment = float(moments[1])
    punctured_moment = punctured_factor * uniform_moment
    matrices[1, 1, 0] = uniform_moment / denominator
    matrices[1, 1, 2] = uniform_moment * (1.0 - 1.0 / denominator)
    matrices[1, 2, 0] = punctured_moment / denominator
    matrices[1, 2, 2] = punctured_moment * (1.0 - 1.0 / denominator)
    return matrices


def weight_conditioned_log_moments_general(
    *, outer_bits: int, zero_row: np.ndarray, active_row: np.ndarray
) -> np.ndarray:
    """General-state version of the exact row-weight coefficient DP."""
    states = zero_row.shape[0]
    if zero_row.shape != (states, states) or active_row.shape != zero_row.shape:
        raise ValueError("row transfers must be equal-size square matrices")

    def log_entries(matrix: np.ndarray) -> np.ndarray:
        result = np.full(matrix.shape, -math.inf)
        positive = matrix > 0.0
        result[positive] = np.log(matrix[positive])
        return result

    def log_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        # left[b,i,k] + right[k,j], reduced over k.
        return np.logaddexp.reduce(
            left[:, :, :, np.newaxis] + right[np.newaxis, np.newaxis, :, :],
            axis=2,
        )

    log_zero = log_entries(zero_row)
    log_active = log_entries(active_row)
    coefficients = np.full(
        (outer_bits + 1, states, states), -math.inf, dtype=np.float64
    )
    for state in range(states):
        coefficients[0, state, state] = 0.0
    for row in range(outer_bits):
        old = coefficients[: row + 1]
        updated = np.full_like(coefficients, -math.inf)
        updated[: row + 1] = np.logaddexp(
            updated[: row + 1], log_matmul(old, log_zero)
        )
        updated[1 : row + 2] = np.logaddexp(
            updated[1 : row + 2], log_matmul(old, log_active)
        )
        coefficients = updated
    moments = np.full(outer_bits + 1, -math.inf)
    for weight in range(outer_bits + 1):
        moments[weight] = (
            float(logsumexp(coefficients[weight, 0, :]))
            - log_choose(outer_bits, weight)
        )
    return moments


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    outer_dimension = (
        args.outer_dimension
        if args.outer_dimension is not None
        else args.outer_bits // 2
    )
    if args.message_bits % outer_dimension:
        raise ValueError("outer dimension must divide message bits")
    outer_blocks = args.message_bits // outer_dimension
    output_bits = outer_blocks * args.outer_bits
    target_distance = math.floor(args.relative_distance * output_bits)
    region_bits = output_bits // args.outer_bits
    epochs = region_bits // args.step_bits
    spectrum = (
        (
            random_linear_extension_spectrum_logs(
                args.outer_spectrum,
                args.outer_bits,
                args.random_linear_extension_bits,
            )
            if args.random_linear_extension
            else load_outer_spectrum_logs(args.outer_spectrum, args.outer_bits)
        )
        if args.outer_spectrum is not None
        else modeled_even_floor_spectrum_logs(
            args.outer_bits, args.modeled_minimum_distance
        )
    )
    nonactivation = np.ones(args.step_bits + 1, dtype=np.float64)
    if args.ideal_random_activation:
        nonactivation[1:] = math.ldexp(1.0, -args.state_bits)
    else:
        activation_payload = json.loads(
            args.activation.read_text(encoding="utf-8")
        )
        for row in activation_payload["by_total_weight"]:
            weight = int(row["total_weight"])
            if weight > args.step_bits:
                continue
            key = "uniform_support_average_distinct_upper_bound"
            if key not in row:
                key = "uniform_256_support_average_distinct_upper_bound"
            nonactivation[weight] = float(row[key])

    best = np.full(args.outer_bits + 1, math.inf)
    best_tilt = np.full(args.outer_bits + 1, math.nan)
    live_spectrum = (
        load_nonzero_spectrum(
            args.live_spectrum, args.step_bits, args.state_bits
        )
        if args.live_spectrum is not None
        else None
    )
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        if args.three_state_refresh:
            if live_spectrum is None:
                raise ValueError("--three-state-refresh requires --live-spectrum")
            impulse = refreshed_three_state_impulse_matrices(
                z=z,
                step_bits=args.step_bits,
                state_bits=args.state_bits,
                live_spectrum=live_spectrum,
            )
        else:
            impulse = splitstate_impulse_matrices(
                z=z,
                step_bits=args.step_bits,
                state_bits=args.state_bits,
                constituent_distance=args.constituent_distance,
                live_moment_order=args.live_moment_order,
                live_model="support-averaged-preaddmul",
                nonactivation=nonactivation,
                live_spectrum=live_spectrum,
            )
        inactive_region, active_region = marked_region(
            impulse[0], impulse[1], epochs
        )
        if args.three_state_refresh:
            moments = weight_conditioned_log_moments_general(
                outer_bits=args.outer_bits,
                zero_row=inactive_region,
                active_row=active_region,
            )
        else:
            moments = weight_conditioned_log_moments(
                outer_bits=args.outer_bits,
                zero_row=inactive_region,
                active_row=active_region,
            )
        values = moments + target_distance * surprisal
        improved = values < best
        best[improved] = values[improved]
        best_tilt[improved] = log_surprisal

    rows = []
    logs = []
    for weight in range(1, args.outer_bits + 1):
        if not math.isfinite(float(spectrum[weight])):
            continue
        inner = min(0.0, float(best[weight]))
        total = math.log(outer_blocks) + float(spectrum[weight]) + inner
        logs.append(total)
        rows.append(
            {
                "outer_weight": weight,
                "best_log_surprisal": float(best_tilt[weight]),
                "inner_log2_upper": inner / LOG2,
                "pointwise_log2_upper": total / LOG2,
                "pointwise_margin_bits": -total / LOG2,
            }
        )
    aggregate = float(logsumexp(np.asarray(logs)))
    dominant = max(rows, key=lambda row: row["pointwise_log2_upper"])
    return {
        "schema": "riffle-splitstate-preaddmul-one-active-v1",
        "candidate": "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "outer_bits": args.outer_bits,
            "outer_dimension": outer_dimension,
            "outer_blocks": outer_blocks,
            "modeled_minimum_distance": args.modeled_minimum_distance,
            "outer_spectrum": (
                str(args.outer_spectrum) if args.outer_spectrum is not None else None
            ),
            "random_linear_extension": args.random_linear_extension,
            "random_linear_extension_bits": (
                args.random_linear_extension_bits
                if args.random_linear_extension
                else 0
            ),
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "activation_model": (
                "ideal-random-linear-map"
                if args.ideal_random_activation
                else str(args.activation)
            ),
            "live_spectrum": (
                str(args.live_spectrum)
                if args.live_spectrum is not None
                else None
            ),
            "epochs_per_region": epochs,
            "live_state_model": (
                "zero-uniform-live-punctured-live-refresh"
                if args.three_state_refresh
                else "zero-live-entrywise-envelope"
            ),
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
        },
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_row": dominant,
        "weight_rows": rows,
        "scope": (
            "Exact one-active placement calculation under the supplied outer "
            f"spectrum. The fixed state code is required to have distance "
            f"{args.constituent_distance}, full output support, and the "
            "factorial-moment profile supplied to the evaluator. Arithmetic "
            "is nearest binary64 and is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument(
        "--outer-dimension",
        type=int,
        help="outer message dimension (default: outer-bits // 2)",
    )
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--outer-spectrum", type=Path)
    parser.add_argument(
        "--random-linear-extension",
        action="store_true",
        help=(
            "treat outer-spectrum as a length-(outer-bits-1) source and "
            "append an independently sampled message-linear bit per block"
        ),
    )
    parser.add_argument(
        "--random-linear-extension-bits",
        type=int,
        default=1,
        help="number of independently sampled appended linear bits",
    )
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=32)
    parser.add_argument("--constituent-distance", type=int, default=20)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0),
    )
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--live-spectrum", type=Path)
    parser.add_argument("--ideal-random-activation", action="store_true")
    parser.add_argument("--three-state-refresh", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.random_linear_extension and args.outer_spectrum is None:
        parser.error("--random-linear-extension requires --outer-spectrum")
    if args.random_linear_extension_bits < 1:
        parser.error("--random-linear-extension-bits must be positive")
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"margin_bits,{payload['aggregate_margin_bits']:.9f}")
    print(
        "dominant_outer_weight,"
        f"{payload['dominant_row']['outer_weight']},"
        "pointwise_margin_bits,"
        f"{payload['dominant_row']['pointwise_margin_bits']:.9f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
