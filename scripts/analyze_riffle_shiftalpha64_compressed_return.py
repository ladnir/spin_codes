#!/usr/bin/env python3
"""Compress packet-accumulator path statistics by coefficient tilts.

For a path with H active packets and m zero prefix states, the exact-gap
certificate contains ``C(N-H+m,m)``.  For every theta in (0,1),

    C(N-H+m,m) <= theta^(-m) (1-theta)^(-(N-H+1)).

The factor ``theta^(-m)`` is applied when the 16-state walk returns to zero.
The dynamic program therefore retains only total input weight, packet support,
and the current four-bit state.  It no longer indexes the return count m.

With ``--compress-weight``, a second coefficient envelope removes the total
binary-weight coordinate:

    [x^W] F(x) <= x^(-W) F(x),  x > 0.

The script analyzes one equal-weight BCH-block profile.  It prints one JSON
receipt and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_bchblockperm_parallelacc_g4_goal02 as kernel  # noqa: E402


PART_BITS = 128
PACKET_BITS = 4
PACKETS_PER_PART = PART_BITS // PACKET_BITS
STATE_COUNT = 1 << PACKET_BITS
DEFAULT_PACKET_POSITIONS = 524_352


def configure_kernel(part_count: int, part_weight: int) -> None:
    kernel.PART_COUNT = part_count
    kernel.PART_WEIGHT = part_weight
    kernel.SUPERBLOCK_BITS = PART_BITS * part_count
    kernel.SUPERBLOCK_WEIGHT = part_weight * part_count
    kernel.PACKETS_PER_SUPERBLOCK = PACKETS_PER_PART * part_count


def evaluate_bound(
    *,
    distance: int,
    packet_positions: int,
    scaled_cost: float,
    zero_scale: float,
) -> dict[str, object]:
    """Evaluate one valid compressed positive-gap certificate."""
    if scaled_cost <= 0.0 or not 0.0 < zero_scale < packet_positions:
        return {"raw_log_bound": math.inf}

    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = kernel.weighted_path_log_counts(state_factor)
    normalization_log = kernel.log_binom(
        kernel.SUPERBLOCK_BITS, kernel.SUPERBLOCK_WEIGHT
    )
    terms: list[tuple[int, float]] = []
    zero_gap_log_factor = -math.log1p(-theta)
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            - normalization_log
            + (packet_positions - support + 1) * zero_gap_log_factor
            - kernel.log_binom(packet_positions, support)
        )
        terms.append((support, term))

    raw_log_bound = scaled_cost + float(logsumexp([term for _, term in terms]))
    unconditional_log_bound = min(0.0, raw_log_bound)
    conditioning_log = kernel.balance_log_probability()
    balanced_log_bound = min(
        0.0, unconditional_log_bound - conditioning_log
    )
    top_terms = sorted(terms, key=lambda item: item[1], reverse=True)[:12]
    return {
        "raw_log_bound": raw_log_bound,
        "scaled_positive_cost": scaled_cost,
        "positive_cost_parameter_u": u,
        "scaled_zero_return_parameter": zero_scale,
        "zero_return_tilt_theta": theta,
        "unconditional_log2_bound": unconditional_log_bound / math.log(2.0),
        "conditioning_penalty_bits": -conditioning_log / math.log(2.0),
        "balanced_profile_log2_bound": balanced_log_bound / math.log(2.0),
        "top_support_terms": [
            {
                "packet_support": support,
                "pre_chernoff_log2_contribution": term / math.log(2.0),
            }
            for support, term in top_terms
        ],
    }


def weighted_path_log_counts_with_weight_tilt(
    state_factor: np.ndarray, weight_tilt: float
) -> np.ndarray:
    """Sum paths by support and state after tilting total binary weight."""
    if weight_tilt <= 0.0:
        return np.full(kernel.PACKETS_PER_SUPERBLOCK + 1, math.inf)
    shape = (kernel.PACKETS_PER_SUPERBLOCK + 1, STATE_COUNT)
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0] = 1.0
    states = np.arange(STATE_COUNT)
    log_scale = 0.0

    for slot in range(kernel.PACKETS_PER_SUPERBLOCK):
        following = np.zeros_like(current)
        support_stop = slot + 1
        for value in range(STATE_COUNT):
            support_increment = int(value != 0)
            source = current[:support_stop, states ^ value]
            if value:
                source = source * state_factor[states]
            following[
                support_increment : support_increment + support_stop,
                states,
            ] += source * weight_tilt ** value.bit_count()
        scale = float(np.max(following))
        if not math.isfinite(scale) or scale <= 0.0:
            return np.full(kernel.PACKETS_PER_SUPERBLOCK + 1, math.inf)
        following /= scale
        log_scale += math.log(scale)
        current = following

    totals = np.sum(current, axis=1)
    result = np.full(kernel.PACKETS_PER_SUPERBLOCK + 1, -math.inf)
    positive = totals > 0.0
    result[positive] = np.log(totals[positive]) + log_scale
    return result


def evaluate_weight_compressed_bound(
    *,
    distance: int,
    packet_positions: int,
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
) -> dict[str, object]:
    """Also remove the exact total-weight dimension by coefficient tilting."""
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = weighted_path_log_counts_with_weight_tilt(
        state_factor, weight_tilt
    )
    normalization_log = kernel.log_binom(
        kernel.SUPERBLOCK_BITS, kernel.SUPERBLOCK_WEIGHT
    )
    weight_coefficient_cost = -kernel.SUPERBLOCK_WEIGHT * math.log(weight_tilt)
    terms: list[tuple[int, float]] = []
    zero_gap_log_factor = -math.log1p(-theta)
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            + weight_coefficient_cost
            - normalization_log
            + (packet_positions - support + 1) * zero_gap_log_factor
            - kernel.log_binom(packet_positions, support)
        )
        terms.append((support, term))

    raw_log_bound = scaled_cost + float(logsumexp([term for _, term in terms]))
    unconditional_log_bound = min(0.0, raw_log_bound)
    conditioning_log = kernel.balance_log_probability()
    balanced_log_bound = min(0.0, unconditional_log_bound - conditioning_log)
    top_terms = sorted(terms, key=lambda item: item[1], reverse=True)[:12]
    return {
        "raw_log_bound": raw_log_bound,
        "scaled_positive_cost": scaled_cost,
        "positive_cost_parameter_u": u,
        "scaled_zero_return_parameter": zero_scale,
        "zero_return_tilt_theta": theta,
        "binary_weight_coefficient_tilt": weight_tilt,
        "unconditional_log2_bound": unconditional_log_bound / math.log(2.0),
        "conditioning_penalty_bits": -conditioning_log / math.log(2.0),
        "balanced_profile_log2_bound": balanced_log_bound / math.log(2.0),
        "top_support_terms": [
            {
                "packet_support": support,
                "pre_chernoff_log2_contribution": term / math.log(2.0),
            }
            for support, term in top_terms
        ],
    }


def optimize_bound(
    *,
    distance: int,
    packet_positions: int,
    start_scaled_cost: float,
    start_zero_scale: float,
    maxiter: int,
) -> dict[str, object]:
    """Optimize in log coordinates; every evaluated point is a valid bound."""
    start = np.log(np.asarray([start_scaled_cost, start_zero_scale]))

    def objective(log_point: np.ndarray) -> float:
        scaled_cost, zero_scale = np.exp(log_point)
        return float(
            evaluate_bound(
                distance=distance,
                packet_positions=packet_positions,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
    )
    scaled_cost, zero_scale = np.exp(result.x)
    selected = evaluate_bound(
        distance=distance,
        packet_positions=packet_positions,
        scaled_cost=float(scaled_cost),
        zero_scale=float(zero_scale),
    )
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations": int(result.nfev),
        }
    )
    return selected


def optimize_weight_compressed_bound(
    *,
    distance: int,
    packet_positions: int,
    start_scaled_cost: float,
    start_zero_scale: float,
    start_weight_tilt: float,
    maxiter: int,
) -> dict[str, object]:
    """Optimize the return, positive-cost, and binary-weight tilts."""
    start = np.log(
        np.asarray([start_scaled_cost, start_zero_scale, start_weight_tilt])
    )

    def objective(log_point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(log_point)
        return float(
            evaluate_weight_compressed_bound(
                distance=distance,
                packet_positions=packet_positions,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
                weight_tilt=float(weight_tilt),
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
    )
    scaled_cost, zero_scale, weight_tilt = np.exp(result.x)
    selected = evaluate_weight_compressed_bound(
        distance=distance,
        packet_positions=packet_positions,
        scaled_cost=float(scaled_cost),
        zero_scale=float(zero_scale),
        weight_tilt=float(weight_tilt),
    )
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations": int(result.nfev),
        }
    )
    return selected


def validate_coefficient_envelope(theta: float) -> str:
    """Exhaustively verify the coefficient inequality on small integers."""
    for inactive in range(21):
        for zero_returns in range(11):
            exact = math.comb(inactive + zero_returns, zero_returns)
            envelope = (
                theta ** (-zero_returns)
                * (1.0 - theta) ** (-(inactive + 1))
            )
            if exact > envelope * (1.0 + 2e-13):
                raise RuntimeError("coefficient-envelope validation failed")
    return "PASS (231 small coefficient cases)"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part-count", type=int, required=True)
    parser.add_argument("--part-weight", type=int, required=True)
    parser.add_argument("--packet-positions", type=int, default=DEFAULT_PACKET_POSITIONS)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument("--scaled-cost", type=float, default=33.7)
    parser.add_argument("--zero-scale", type=float, default=18.0)
    parser.add_argument("--weight-tilt", type=float)
    parser.add_argument("--compress-weight", action="store_true")
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--optimizer-maxiter", type=int, default=30)
    args = parser.parse_args()

    if args.part_count < 1:
        raise ValueError("part count must be positive")
    if not 1 <= args.part_weight <= PART_BITS:
        raise ValueError("part weight lies outside one BCH block")
    if args.distance <= 0:
        raise ValueError("distance must be positive")
    if args.packet_positions < args.part_count * PACKETS_PER_PART:
        raise ValueError("packet count is smaller than the active superblock")
    if args.scaled_cost <= 0.0 or args.zero_scale <= 0.0:
        raise ValueError("tilt parameters must be positive")

    configure_kernel(args.part_count, args.part_weight)
    default_weight_tilt = args.part_weight / (PART_BITS - args.part_weight)
    weight_tilt = (
        args.weight_tilt if args.weight_tilt is not None else default_weight_tilt
    )
    if args.compress_weight:
        result = (
            optimize_weight_compressed_bound(
                distance=args.distance,
                packet_positions=args.packet_positions,
                start_scaled_cost=args.scaled_cost,
                start_zero_scale=args.zero_scale,
                start_weight_tilt=weight_tilt,
                maxiter=args.optimizer_maxiter,
            )
            if args.optimize
            else evaluate_weight_compressed_bound(
                distance=args.distance,
                packet_positions=args.packet_positions,
                scaled_cost=args.scaled_cost,
                zero_scale=args.zero_scale,
                weight_tilt=weight_tilt,
            )
        )
    else:
        result = (
            optimize_bound(
                distance=args.distance,
                packet_positions=args.packet_positions,
                start_scaled_cost=args.scaled_cost,
                start_zero_scale=args.zero_scale,
                maxiter=args.optimizer_maxiter,
            )
            if args.optimize
            else evaluate_bound(
                distance=args.distance,
                packet_positions=args.packet_positions,
                scaled_cost=args.scaled_cost,
                zero_scale=args.zero_scale,
            )
        )
    theta = float(result["zero_return_tilt_theta"])
    cells = (
        (kernel.PACKETS_PER_SUPERBLOCK + 1) * STATE_COUNT
        if args.compress_weight
        else (
            (kernel.SUPERBLOCK_WEIGHT + 1)
            * (kernel.PACKETS_PER_SUPERBLOCK + 1)
            * STATE_COUNT
        )
    )
    payload = {
        "schema": "riffle-shiftalpha64-compressed-return-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "profile": {
            "active_bch_blocks": args.part_count,
            "bch_block_weights": [args.part_weight] * args.part_count,
            "total_binary_input_weight": kernel.SUPERBLOCK_WEIGHT,
        },
        "parameters": {
            "global_packet_positions": args.packet_positions,
            "global_binary_length": PACKET_BITS * args.packet_positions,
            "failure_weight_inclusive": args.distance,
        },
        "compression": {
            "removed_dimensions": (
                ["zero prefix-state count m", "total superblock bit weight"]
                if args.compress_weight
                else ["zero prefix-state count m"]
            ),
            "retained_dimensions": [
                "active packet support H",
                "current state in GF(2)^4",
            ]
            + ([] if args.compress_weight else ["total superblock bit weight"]),
            "state_cells_per_array": cells,
            "state_bytes_for_two_float64_arrays": 2 * cells * 8,
            "coefficient_envelope": (
                "C(N-H+m,m) <= theta^(-m) "
                "(1-theta)^(-(N-H+1))"
            ),
            "binary_weight_coefficient_envelope": (
                "[x^W] F(x) <= x^(-W) F(x)"
                if args.compress_weight
                else "NOT_USED"
            ),
            "equal_block_normalization": (
                "After weight compression, pooled normalization and the "
                "balance factor cancel to C(128,h)^(-s)."
                if args.compress_weight
                else "The pooled-slice conditioning factor remains explicit."
            ),
        },
        "compressed_bound": result,
        "validation": {
            "coefficient_envelope": validate_coefficient_envelope(theta),
            "validity_under_optimizer_warning": (
                "PASS: the selected numerical parameters define a valid bound "
                "even when the optimizer stops early"
            ),
        },
        "scope": (
            "The result is a rigorous transfer for the declared equal-weight "
            "profile. It does not sum the outer spectrum."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
