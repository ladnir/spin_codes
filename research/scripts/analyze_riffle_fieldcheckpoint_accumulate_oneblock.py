#!/usr/bin/env python3
"""One-active diagnostic for FieldCheckpointAccumulate.

The outer geometry is BCHPerm-TransposeBitShuffle.  One active outer word
places one impulse in every selected region.  The checkpoint accumulator has
two post-checkpoint state classes: zero and uniform nonzero.  Exact epoch
moments give a two-by-two transfer, and the existing region-weight dynamic
program averages the selected regions without replacement.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
    weight_conditioned_log_moments,
)
from analyze_riffle_striped_random_outer import LOG2


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal01_one_active.json"
)


def geometric_sum(z: float, first: int, last: int) -> float:
    return math.fsum(z**exponent for exponent in range(first, last + 1))


def epoch_matrices(state_bits: int, epoch_bits: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    """Return exact zero-input and one-impulse epoch moment matrices."""
    if epoch_bits % state_bits:
        raise ValueError("state width must divide the epoch length")
    visits = epoch_bits // state_bits
    denominator = math.ldexp(1.0, state_bits) - 1.0
    lane_monomial = z**visits

    no_input_numerator = math.expm1(state_bits * math.log1p(lane_monomial))
    no_input_live = no_input_numerator / denominator
    zero_epoch = np.asarray(((1.0, 0.0), (0.0, no_input_live)), dtype=np.float64)

    zero_to_live = geometric_sum(z, 1, visits) / visits
    live_to_zero = geometric_sum(z, 0, visits - 1) / (visits * denominator)
    other_lane_numerator = math.expm1(
        (state_bits - 1) * math.log1p(lane_monomial)
    )
    endpoint_sum = math.fsum(
        z ** (visits - round_index) + z**round_index
        for round_index in range(visits)
    )
    live_to_live = other_lane_numerator * endpoint_sum / (visits * denominator)
    one_impulse = np.asarray(
        ((0.0, zero_to_live), (live_to_zero, live_to_live)),
        dtype=np.float64,
    )
    return zero_epoch, one_impulse


def _lane_prefix_polynomials(visits: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    """Sum the two starting-bit moments by the number of lane impulses."""
    start_zero = np.zeros(visits + 1, dtype=np.float64)
    start_one = np.zeros(visits + 1, dtype=np.float64)
    for support in range(1 << visits):
        parity = 0
        prefix_weight = 0
        for visit in range(visits):
            parity ^= (support >> visit) & 1
            prefix_weight += parity
        impulses = support.bit_count()
        start_zero[impulses] += z**prefix_weight
        start_one[impulses] += z ** (visits - prefix_weight)
    return start_zero, start_one


def _fixed_base_power(base: np.ndarray, exponent: int) -> np.ndarray:
    """Raise a short coefficient vector by repeated fixed-base convolution."""
    result = np.asarray((1.0,), dtype=np.float64)
    for _ in range(exponent):
        result = np.convolve(result, base)
    return result


def occupancy_epoch_matrices(
    state_bits: int, epoch_bits: int, z: float
) -> list[np.ndarray]:
    """Return the exact averaged epoch transfer for every impulse count.

    The support conditioned on count r is a uniform r-subset of the epoch.
    Coefficients are scaled by 2^-epoch_bits before polynomial powers so the
    target s=64, epoch_bits=1024 calculation stays in binary64 range.
    """
    if epoch_bits % state_bits:
        raise ValueError("state width must divide the epoch length")
    visits = epoch_bits // state_bits
    start_zero, start_one = _lane_prefix_polynomials(visits, z)
    even_start_zero = start_zero.copy()
    even_start_zero[1::2] = 0.0
    odd_start_one = start_one.copy()
    odd_start_one[0::2] = 0.0

    lane_scale = math.ldexp(1.0, -visits)
    start_zero_power = _fixed_base_power(start_zero * lane_scale, state_bits)
    even_zero_power = _fixed_base_power(
        even_start_zero * lane_scale, state_bits
    )
    endpoint_start_power = _fixed_base_power(
        (even_start_zero + odd_start_one) * lane_scale, state_bits
    )
    all_start_power = _fixed_base_power(
        (start_zero + start_one) * (lane_scale / 2.0), state_bits
    )

    live_states = math.ldexp(1.0, state_bits) - 1.0
    all_state_scale = math.ldexp(1.0, state_bits)
    matrices = []
    for impulses in range(epoch_bits + 1):
        support_probability = math.ldexp(
            float(math.comb(epoch_bits, impulses)), -epoch_bits
        )
        zero_to_zero = even_zero_power[impulses] / support_probability
        zero_total = start_zero_power[impulses] / support_probability
        live_to_zero = (
            endpoint_start_power[impulses] - even_zero_power[impulses]
        ) / (live_states * support_probability)
        live_total = (
            all_state_scale * all_start_power[impulses]
            - start_zero_power[impulses]
        ) / (live_states * support_probability)
        matrix = np.asarray(
            (
                (zero_to_zero, zero_total - zero_to_zero),
                (live_to_zero, live_total - live_to_zero),
            ),
            dtype=np.float64,
        )
        # Polynomial subtraction can leave ulp-scale negative entries.
        matrix[(matrix < 0.0) & (matrix > -2e-13)] = 0.0
        matrices.append(matrix)
    return matrices


def region_matrices(
    state_bits: int, epoch_bits: int, epochs_per_region: int, z: float
) -> tuple[np.ndarray, np.ndarray]:
    zero_epoch, one_impulse = epoch_matrices(state_bits, epoch_bits, z)
    zero_powers = [np.eye(2)]
    for _ in range(epochs_per_region):
        zero_powers.append(zero_powers[-1] @ zero_epoch)
    inactive = zero_powers[epochs_per_region]
    active = np.zeros((2, 2), dtype=np.float64)
    for epoch in range(epochs_per_region):
        active += zero_powers[epoch] @ one_impulse @ zero_powers[
            epochs_per_region - 1 - epoch
        ]
    active /= epochs_per_region
    return inactive, active


def brute_epoch_matrices(state_bits: int, epoch_bits: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    visits = epoch_bits // state_bits
    denominator = (1 << state_bits) - 1
    zero = np.zeros((2, 2), dtype=np.float64)
    one = np.zeros((2, 2), dtype=np.float64)
    zero[0, 0] = 1.0

    for state in range(1, 1 << state_bits):
        output_weight = visits * state.bit_count()
        zero[1, 1] += z**output_weight / denominator

    for position in range(epoch_bits):
        lane = position % state_bits
        zero_state = 0
        zero_weight = 0
        for index in range(epoch_bits):
            if index == position:
                zero_state ^= 1 << lane
            zero_weight += (zero_state >> (index % state_bits)) & 1
        one[0, 1] += z**zero_weight / epoch_bits

        for state in range(1, 1 << state_bits):
            current = state
            output_weight = 0
            for index in range(epoch_bits):
                if index == position:
                    current ^= 1 << lane
                output_weight += (current >> (index % state_bits)) & 1
            next_class = 0 if current == 0 else 1
            one[1, next_class] += z**output_weight / (epoch_bits * denominator)
    return zero, one


def brute_occupancy_epoch_matrix(
    state_bits: int, epoch_bits: int, impulses: int, z: float
) -> np.ndarray:
    """Exhaustive reference for a uniform fixed-size epoch support."""
    support_count = math.comb(epoch_bits, impulses)
    live_states = (1 << state_bits) - 1
    matrix = np.zeros((2, 2), dtype=np.float64)
    for positions in itertools.combinations(range(epoch_bits), impulses):
        support = set(positions)
        for start_class, states in (
            (0, (0,)),
            (1, range(1, 1 << state_bits)),
        ):
            denominator = support_count * (1 if start_class == 0 else live_states)
            for state in states:
                current = state
                output_weight = 0
                for index in range(epoch_bits):
                    lane = index % state_bits
                    if index in support:
                        current ^= 1 << lane
                    output_weight += (current >> lane) & 1
                end_class = 0 if current == 0 else 1
                matrix[start_class, end_class] += z**output_weight / denominator
    return matrix


def self_test() -> dict[str, float]:
    maximum_error = 0.0
    maximum_occupancy_error = 0.0
    for z in (0.41, 0.73, 1.0):
        exact = epoch_matrices(3, 6, z)
        brute = brute_epoch_matrices(3, 6, z)
        for exact_matrix, brute_matrix in zip(exact, brute):
            maximum_error = max(
                maximum_error,
                float(np.max(np.abs(exact_matrix - brute_matrix))),
            )
        occupancy = occupancy_epoch_matrices(3, 6, z)
        for impulses, exact_matrix in enumerate(occupancy):
            brute_matrix = brute_occupancy_epoch_matrix(3, 6, impulses, z)
            maximum_occupancy_error = max(
                maximum_occupancy_error,
                float(np.max(np.abs(exact_matrix - brute_matrix))),
            )
    if maximum_error > 2e-14:
        raise AssertionError("closed epoch transfer disagrees with brute force")
    if maximum_occupancy_error > 2e-14:
        raise AssertionError("occupancy transfer disagrees with brute force")

    zero, active = region_matrices(4, 8, 3, 1.0)
    stochastic_error = max(
        float(np.max(np.abs(zero.sum(axis=1) - 1.0))),
        float(np.max(np.abs(active.sum(axis=1) - 1.0))),
    )
    if stochastic_error > 2e-14:
        raise AssertionError("region transfer is not stochastic at z=1")
    return {
        "maximum_epoch_bruteforce_error": maximum_error,
        "maximum_occupancy_bruteforce_error": maximum_occupancy_error,
        "maximum_region_stochastic_error": stochastic_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    outer_dimension = args.outer_bits // 2
    if args.message_bits % outer_dimension:
        raise ValueError("outer dimension must divide the message length")
    outer_blocks = args.message_bits // outer_dimension
    output_bits = 2 * args.message_bits
    region_bits = output_bits // args.outer_bits
    epoch_bits = args.step_bits * args.checkpoint_steps
    if region_bits % epoch_bits:
        raise ValueError("checkpoint epochs must divide every transposed region")
    if epoch_bits % args.state_bits:
        raise ValueError("state width must divide the checkpoint epoch")
    epochs_per_region = region_bits // epoch_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum_logs = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )

    best_inner = np.full(args.outer_bits + 1, math.inf)
    best_surprisal = np.full(args.outer_bits + 1, math.nan)
    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        inactive, active = region_matrices(
            args.state_bits, epoch_bits, epochs_per_region, z
        )
        moments = weight_conditioned_log_moments(
            outer_bits=args.outer_bits,
            zero_row=inactive,
            active_row=active,
        )
        candidates = moments + distance * surprisal
        improved = candidates < best_inner
        best_inner[improved] = candidates[improved]
        best_surprisal[improved] = log_surprisal
        print(
            f"grid,{grid_index + 1},{grid_count},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    log_block_choices = math.log(outer_blocks)
    rows = []
    contributions = []
    for weight in range(1, args.outer_bits + 1):
        if not math.isfinite(float(spectrum_logs[weight])):
            continue
        inner_log = min(0.0, float(best_inner[weight]))
        contribution = log_block_choices + float(spectrum_logs[weight]) + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "outer_weight": weight,
                "spectrum_log2": float(spectrum_logs[weight]) / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "log_output_surprisal": float(best_surprisal[weight]),
            }
        )

    one_active_log = float(logsumexp(np.asarray(contributions)))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "riffle-fieldcheckpoint-accumulate-one-active-v1",
        "construction": "FieldCheckpointAccumulate",
        "probability_space": {
            "outer_spectrum": (
                "modeled real-valued complement-symmetric even "
                f"[{args.outer_bits},{outer_dimension},"
                f"{args.modeled_minimum_distance}]-shaped spectrum"
            ),
            "outer_coordinate_permutations": "independent by outer block",
            "region_bit_permutations": "independent by transposed region",
            "checkpoint_multipliers": (
                f"independent uniform nonzero elements of GF(2^{args.state_bits})"
            ),
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "region_bits": region_bits,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "checkpoint_steps": args.checkpoint_steps,
            "epoch_bits": epoch_bits,
            "epochs_per_region": epochs_per_region,
            "relative_distance": args.relative_distance,
            "distance": distance,
        },
        "self_test": self_test(),
        "one_active_log2_upper": one_active_log / LOG2,
        "one_active_lambda_bits_lower_float": -one_active_log / LOG2,
        "dominant_weight": dominant,
        "weight_rows": rows,
        "scope": (
            "One active outer block only. The outer spectrum is a real-valued "
            "model. Chernoff optimization uses an ordinary floating-point grid."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--step-bits", type=int, default=32)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--checkpoint-steps", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "one_active_lambda_bits,"
        f"{float(payload['one_active_lambda_bits_lower_float']):.6f}",
        flush=True,
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
