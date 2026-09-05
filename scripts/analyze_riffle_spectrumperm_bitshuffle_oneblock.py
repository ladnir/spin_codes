#!/usr/bin/env python3
"""One-active spectrum evaluator for SpectrumPerm-TransposeBitShuffle.

For outer weight w, the per-block coordinate permutation selects a uniform
w-subset of the B transposed rows.  A matrix-polynomial dynamic program
computes the exact averaged row-order transfer.  Chernoff tilts then bound the
inner low-weight event.  The default outer spectrum is the fractional ideal
random-code benchmark.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_striped_random_outer import (
    LOG2,
    log_choose,
    log_matmul_batch,
    log_matrix_entries,
    log_two_power_minus_one,
)
from analyze_riffle_transpose_packetshuffle_sparse_curve import one_marked_row


DEFAULT_OUTPUT = Path(
    "constructions/riffle_spectrumperm_transpose_bitshuffle_randomstepconv/"
    "receipts/goal01_random_spectrum.json"
)


def step_matrices(packet_bits: int, sigma: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    state_zero = math.ldexp(1.0, -sigma)
    output_moment = math.ldexp((1.0 + z) ** packet_bits, -packet_bits)
    off = state_zero * output_moment
    live = (1.0 - state_zero) * output_moment
    zero = np.asarray(((1.0, 0.0), (off, live)), dtype=np.float64)
    nonzero = np.asarray(((off, live), (off, live)), dtype=np.float64)
    return zero, nonzero


def row_matrices(
    packet_bits: int, sigma: int, positions: int, z: float
) -> tuple[np.ndarray, np.ndarray]:
    zero, nonzero = step_matrices(packet_bits, sigma, z)
    zero_row = np.linalg.matrix_power(zero, positions)
    active_row = one_marked_row(zero, nonzero, positions)
    return zero_row, active_row


def weight_conditioned_log_moments(
    *, outer_bits: int, zero_row: np.ndarray, active_row: np.ndarray
) -> np.ndarray:
    """Return log E[z^W | outer weight w] for every w."""
    log_zero = log_matrix_entries(tuple(zero_row.reshape(4)))
    log_active = log_matrix_entries(tuple(active_row.reshape(4)))
    coefficients = np.full((outer_bits + 1, 2, 2), -math.inf)
    coefficients[0, 0, 0] = 0.0
    coefficients[0, 1, 1] = 0.0
    for row in range(outer_bits):
        old = coefficients[: row + 1]
        updated = np.full_like(coefficients, -math.inf)
        zero_terms = log_matmul_batch(old, log_zero)
        active_terms = log_matmul_batch(old, log_active)
        updated[: row + 1] = np.logaddexp(updated[: row + 1], zero_terms)
        updated[1 : row + 2] = np.logaddexp(
            updated[1 : row + 2], active_terms
        )
        coefficients = updated
    moments = np.full(outer_bits + 1, -math.inf)
    for weight in range(outer_bits + 1):
        moments[weight] = (
            np.logaddexp(coefficients[weight, 0, 0], coefficients[weight, 0, 1])
            - log_choose(outer_bits, weight)
        )
    return moments


def brute_weight_average(
    zero_row: np.ndarray, active_row: np.ndarray, outer_bits: int, weight: int
) -> float:
    total = 0.0
    count = 0
    for active_positions in itertools.combinations(range(outer_bits), weight):
        active = set(active_positions)
        product = np.eye(2)
        for row in range(outer_bits):
            product = product @ (active_row if row in active else zero_row)
        total += product[0, 0] + product[0, 1]
        count += 1
    return total / count


def self_test() -> dict[str, float]:
    zero_row, active_row = row_matrices(2, 3, 3, 0.61)
    recovered = weight_conditioned_log_moments(
        outer_bits=6, zero_row=zero_row, active_row=active_row
    )
    maximum_error = 0.0
    for weight in range(7):
        brute = brute_weight_average(zero_row, active_row, 6, weight)
        maximum_error = max(maximum_error, abs(math.exp(recovered[weight]) - brute))
    if maximum_error > 4e-13:
        raise AssertionError("row-weight coefficient DP changed the exact average")
    ideal = ideal_random_spectrum_logs(32)
    modeled = modeled_even_floor_spectrum_logs(32, 6)
    target = log_two_power_minus_one(16)
    ideal_mass_error = abs(float(logsumexp(ideal)) - target)
    modeled_mass_error = abs(float(logsumexp(modeled)) - target)
    if max(ideal_mass_error, modeled_mass_error) > 3e-13:
        raise AssertionError("spectrum proxy mass changed")
    return {
        "maximum_row_weight_error": maximum_error,
        "ideal_spectrum_log_mass_error": ideal_mass_error,
        "modeled_spectrum_log_mass_error": modeled_mass_error,
    }


def ideal_random_spectrum_logs(outer_bits: int) -> np.ndarray:
    dimension = outer_bits // 2
    log_scale = log_two_power_minus_one(dimension) - log_two_power_minus_one(outer_bits)
    values = np.full(outer_bits + 1, -math.inf)
    for weight in range(1, outer_bits + 1):
        values[weight] = log_scale + log_choose(outer_bits, weight)
    return values


def modeled_even_floor_spectrum_logs(
    outer_bits: int, minimum_distance: int
) -> np.ndarray:
    """Return the real-valued complement-symmetric d-min proxy.

    The two endpoint words have multiplicity one.  The remaining mass is
    proportional to C(B,w) on even weights d <= w <= B-d.
    """
    dimension = outer_bits // 2
    if minimum_distance % 2 or not 2 <= minimum_distance < outer_bits // 2:
        raise ValueError("the modeled even minimum distance is invalid")
    weights = list(range(minimum_distance, outer_bits - minimum_distance + 1, 2))
    raw = np.asarray([log_choose(outer_bits, weight) for weight in weights])
    interior_log_mass = log_two_power_minus_one(dimension) + math.log1p(
        -math.exp(-log_two_power_minus_one(dimension))
    )
    # The expression above is log(2^dimension-2), evaluated without forming
    # the enormous integer as a float.
    normalizer = float(logsumexp(raw))
    values = np.full(outer_bits + 1, -math.inf)
    for weight, raw_log in zip(weights, raw):
        values[weight] = float(raw_log) + interior_log_mass - normalizer
    values[outer_bits] = 0.0
    return values


def load_spectrum_logs(path: Path, outer_bits: int) -> tuple[np.ndarray, str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    raw = payload["weight_counts"]
    counts = [0] * (outer_bits + 1)
    if isinstance(raw, list):
        if len(raw) != outer_bits + 1:
            raise ValueError("spectrum list has the wrong length")
        counts = [int(value) for value in raw]
    else:
        for weight, value in raw.items():
            counts[int(weight)] = int(value)
    if counts[0] != 1:
        raise ValueError("a linear-code spectrum must have A_0=1")
    if sum(counts) != 1 << (outer_bits // 2):
        raise ValueError("spectrum multiplicities do not sum to 2^(B/2)")
    values = np.full(outer_bits + 1, -math.inf)
    for weight in range(1, outer_bits + 1):
        if counts[weight]:
            values[weight] = math.log(counts[weight])
    return values, str(path)


def evaluate(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance: float,
    grid_min: float,
    grid_max: float,
    grid_step: float,
    spectrum_path: Path | None,
    modeled_minimum_distance: int | None,
) -> dict[str, object]:
    dimension = outer_bits // 2
    if message_bits % dimension:
        raise ValueError("outer dimension must divide the message length")
    blocks = message_bits // dimension
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    positions = blocks // packet_bits
    output_bits = outer_bits * positions * packet_bits
    distance = math.floor(relative_distance * output_bits)
    if spectrum_path is not None and modeled_minimum_distance is not None:
        raise ValueError("choose an exact spectrum or a modeled minimum distance")
    if spectrum_path is not None:
        spectrum_logs, spectrum_source = load_spectrum_logs(spectrum_path, outer_bits)
    elif modeled_minimum_distance is not None:
        spectrum_logs = modeled_even_floor_spectrum_logs(
            outer_bits, modeled_minimum_distance
        )
        spectrum_source = (
            "modeled_real_valued_even_random_like_spectrum_with_"
            f"minimum_distance_{modeled_minimum_distance}"
        )
    else:
        spectrum_logs = ideal_random_spectrum_logs(outer_bits)
        spectrum_source = "ideal_expected_random_linear_spectrum"

    best_inner = np.full(outer_bits + 1, math.inf)
    best_surprisal = np.full(outer_bits + 1, math.nan)
    grid_count = int(math.floor((grid_max - grid_min) / grid_step + 0.5)) + 1
    for grid_index in range(grid_count):
        log_surprisal = grid_min + grid_index * grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero_row, active_row = row_matrices(packet_bits, sigma, positions, z)
        moments = weight_conditioned_log_moments(
            outer_bits=outer_bits,
            zero_row=zero_row,
            active_row=active_row,
        )
        candidates = moments + distance * surprisal
        improved = candidates < best_inner
        best_inner[improved] = candidates[improved]
        best_surprisal[improved] = log_surprisal
        if grid_index % max(1, grid_count // 10) == 0:
            print(
                f"grid,{grid_index + 1},{grid_count},"
                f"log_surprisal,{log_surprisal:.6f}",
                flush=True,
            )

    rows = []
    pointwise_logs = []
    log_block_choices = math.log(blocks)
    for weight in range(1, outer_bits + 1):
        if not math.isfinite(float(spectrum_logs[weight])):
            continue
        inner_log = min(0.0, float(best_inner[weight]))
        contribution = log_block_choices + float(spectrum_logs[weight]) + inner_log
        pointwise_logs.append(contribution)
        rows.append(
            {
                "outer_weight": weight,
                "spectrum_log2": float(spectrum_logs[weight]) / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "log_output_surprisal": float(best_surprisal[weight]),
            }
        )
    total_log2 = float(logsumexp(pointwise_logs) / LOG2)
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "message_bits": message_bits,
        "outer_bits": outer_bits,
        "outer_dimension_bits": dimension,
        "outer_blocks": blocks,
        "packet_bits": packet_bits,
        "packet_positions_per_row": positions,
        "sigma": sigma,
        "output_bits": output_bits,
        "distance": distance,
        "spectrum_source": spectrum_source,
        "output_grid": {
            "minimum": grid_min,
            "maximum": grid_max,
            "step": grid_step,
            "count": grid_count,
        },
        "one_active_log2_expected_bad_upper": total_log2,
        "one_active_lambda_bits_lower_float": -total_log2,
        "dominant_weight": dominant,
        "weight_rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--g", type=int, default=4)
    parser.add_argument("--sigma", type=int, default=18)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--spectrum-json", type=Path)
    parser.add_argument("--modeled-minimum-distance", type=int)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    checks = self_test()
    if args.self_test_only:
        print(json.dumps(checks, indent=2, sort_keys=True))
        return
    result = evaluate(
        message_bits=args.message_bits,
        outer_bits=args.outer_bits,
        packet_bits=args.g,
        sigma=args.sigma,
        relative_distance=args.relative_distance,
        grid_min=args.grid_min,
        grid_max=args.grid_max,
        grid_step=args.grid_step,
        spectrum_path=args.spectrum_json,
        modeled_minimum_distance=args.modeled_minimum_distance,
    )
    payload = {
        "schema": "riffle-spectrumperm-bitshuffle-oneblock-v1",
        "candidate": "Riffle SpectrumPerm-TransposeBitShuffle-RandomStepConv g",
        "evidence_label": "EXACT_ROW_WEIGHT_AVERAGE_AND_NUMERICAL_CHERNOFF_GRID",
        "self_test": checks,
        "case": result,
        "scope": "One active outer block only. The default fractional random spectrum is an ideal benchmark, not a fixed constituent code.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "one_active_lambda_bits_lower_float": result[
                    "one_active_lambda_bits_lower_float"
                ],
                "dominant_weight": result["dominant_weight"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
