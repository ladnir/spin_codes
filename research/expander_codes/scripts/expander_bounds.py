#!/usr/bin/env python3
"""Enumerator-based first-moment diagnostics for EA and EC codes.

The script evaluates a Chernoff bound obtained from the exact output-weight
generating matrix.  It is intended for parameter exploration; the paper proof
should state the underlying coefficient enumerator before using this bound.
"""

from __future__ import annotations

import argparse
import math
from collections.abc import Iterable

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln


def log_binom(n: int, r: int) -> float:
    if r < 0 or r > n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(r + 1) - gammaln(n - r + 1))


def activation_probability(row_density: float, message_weight: int) -> float:
    """Probability that one Bernoulli-expander coordinate is nonzero."""
    # log1p avoids losing the small difference from one when row_density is tiny.
    parity_bias = math.exp(message_weight * math.log1p(-2.0 * row_density))
    return 0.5 * (1.0 - parity_bias)


def ea_matrix(q: float, z: float) -> np.ndarray:
    """Exact accumulator output-weight generating matrix."""
    return np.array(((1.0 - q, q * z), (q, (1.0 - q) * z)), dtype=float)


def ec_enumerator_matrix(
    input_marker: float,
    output_marker: float,
    memory: int,
    wrapping: bool,
) -> np.ndarray:
    """Exact bivariate input--output enumerator matrix for random convolution."""
    size = memory + 1
    matrix = np.zeros((size, size), dtype=float)
    uniform_states = memory - 1 if wrapping else memory
    for ell in range(uniform_states):
        matrix[ell, 0] = 0.5 * (1.0 + input_marker) * output_marker
        matrix[ell, ell + 1] = 0.5 * (1.0 + input_marker)
    if wrapping:
        # In state 10^(m-1), the fixed oldest tap flips the input bit.
        matrix[memory - 1, 0] = output_marker
        matrix[memory - 1, memory] = input_marker
    matrix[memory, 0] = input_marker * output_marker
    matrix[memory, memory] = 1.0
    return matrix


def ec_matrix(q: float, z: float, memory: int, wrapping: bool) -> np.ndarray:
    """Exact reduced EC generating matrix for uniform convolution taps.

    State ell<m records ell trailing zero outputs.  State m is the all-zero
    convolution state.  A transition producing output one receives weight z.
    """
    size = memory + 1
    matrix = np.zeros((size, size), dtype=float)
    uniform_states = memory - 1 if wrapping else memory
    for ell in range(uniform_states):
        matrix[ell, 0] += 0.5 * z
        matrix[ell, ell + 1] += 0.5
    if wrapping:
        # When only the oldest state bit is one, its tap is fixed to one.
        matrix[memory - 1, 0] = (1.0 - q) * z
        matrix[memory - 1, memory] = q
    matrix[memory, 0] = q * z
    matrix[memory, memory] = 1.0 - q
    return matrix


def normalized(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    scale = float(np.max(matrix))
    if scale == 0.0:
        return matrix, -math.inf
    return matrix / scale, math.log(scale)


def log_weight_mgf(matrix: np.ndarray, length: int, start_state: int) -> float:
    """Compute log(e_start^T matrix^length 1) with scaled squaring."""
    result = np.eye(matrix.shape[0], dtype=float)
    result_log_scale = 0.0
    power, power_log_scale = normalized(matrix)
    exponent = length

    while exponent:
        if exponent & 1:
            result, local_scale = normalized(result @ power)
            result_log_scale += power_log_scale + local_scale
        exponent >>= 1
        if exponent:
            power, local_scale = normalized(power @ power)
            power_log_scale = 2.0 * power_log_scale + local_scale

    endpoint_mass = float(np.sum(result[start_state, :]))
    return result_log_scale + math.log(endpoint_mass)


def log_tail_bound(
    *,
    code: str,
    q: float,
    length: int,
    cutoff: int,
    memory: int,
    wrapping: bool,
) -> tuple[float, float]:
    """Return (log bound, optimizing z) for Pr[output weight <= cutoff]."""
    if cutoff < 0:
        return -math.inf, 0.0
    if cutoff >= length:
        return 0.0, 1.0
    start_state = 0 if code == "ea" else memory

    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        matrix = (
            ea_matrix(q, z)
            if code == "ea"
            else ec_matrix(q, z, memory, wrapping)
        )
        return log_weight_mgf(matrix, length, start_state) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-40.0, 0.0),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 160},
    )
    # z=1 gives the trivial probability bound one and protects against
    # numerical optimizer noise near the boundary.
    best = min(0.0, float(optimum.fun))
    return best, math.exp(float(optimum.x))


def log_slice_tail_bound(
    *,
    input_weight: int,
    length: int,
    cutoff: int,
    memory: int,
    wrapping: bool,
) -> tuple[float, float, float]:
    """Chernoff-bound an EC output tail for a uniform fixed-weight input.

    Returns ``(log_bound, input_marker, output_marker)``.  The two markers
    extract the input shell and the low-output tail from the exact bivariate
    transfer matrix.
    """
    if not (0 <= input_weight <= length):
        raise ValueError("input_weight must lie in [0,length]")
    if input_weight == 0:
        return 0.0, 1.0, 1.0
    if cutoff < 0:
        return -math.inf, 1.0, 1.0
    if cutoff >= length:
        return 0.0, 1.0, 1.0

    def objective(markers: np.ndarray) -> float:
        log_x, log_z = float(markers[0]), float(markers[1])
        matrix = ec_enumerator_matrix(
            math.exp(log_x), math.exp(log_z), memory, wrapping
        )
        return (
            log_weight_mgf(matrix, length, memory)
            - input_weight * log_x
            - cutoff * log_z
            - log_binom(length, input_weight)
        )

    shell_saddle = math.log(input_weight / (length - input_weight)) if input_weight < length else 20.0
    starts = (
        (shell_saddle, math.log(max(cutoff, 1) / max(length - cutoff, 1))),
        (min(20.0, shell_saddle + 2.0), -max(1e-7, 4.0 * input_weight / length)),
        (shell_saddle, -1e-5),
    )
    best_value = 0.0
    best_markers = (1.0, 1.0)
    for start in starts:
        optimum = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="Nelder-Mead",
            options={"xatol": 1e-10, "fatol": 1e-9, "maxiter": 500},
        )
        if float(optimum.fun) < best_value:
            best_value = float(optimum.fun)
            best_markers = (
                math.exp(float(optimum.x[0])),
                math.exp(float(optimum.x[1])),
            )
    return best_value, best_markers[0], best_markers[1]


def scan_weights(k: int, dense_prefix: int = 128, geometric_points: int = 220) -> list[int]:
    values = set(range(1, min(k, dense_prefix) + 1))
    if k > dense_prefix:
        logs = np.linspace(math.log(dense_prefix + 1), math.log(k), geometric_points)
        values.update(int(round(math.exp(value))) for value in logs)
    values.update((max(1, k // 4), max(1, k // 2), max(1, 3 * k // 4), k))
    return sorted(value for value in values if 1 <= value <= k)


def evaluate(args: argparse.Namespace, weights: Iterable[int]) -> list[tuple[float, int, float, float]]:
    if args.systematic:
        length = int(round(args.k * (1.0 / args.rate - 1.0)))
        total_length = args.k + length
    else:
        length = int(round(args.k / args.rate))
        total_length = length
    row_density = args.row_weight / length
    results: list[tuple[float, int, float, float]] = []
    for message_weight in weights:
        cutoff = math.floor(args.delta * total_length) - (
            message_weight if args.systematic else 0
        )
        q = activation_probability(row_density, message_weight)
        log_tail, z = log_tail_bound(
            code=args.code,
            q=q,
            length=length,
            cutoff=cutoff,
            memory=args.memory,
            wrapping=args.wrapping,
        )
        log_first_moment_term = log_binom(args.k, message_weight) + log_tail
        results.append((log_first_moment_term, message_weight, q, z))
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--code", choices=("ea", "ec"), required=True)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--rate", type=float, required=True)
    parser.add_argument("--row-weight", type=float, required=True)
    parser.add_argument("--memory", type=int, default=25)
    parser.add_argument("--wrapping", action="store_true")
    parser.add_argument("--systematic", action="store_true")
    parser.add_argument("--delta", type=float, default=0.05)
    parser.add_argument("--dense-prefix", type=int, default=128)
    parser.add_argument("--geometric-points", type=int, default=220)
    parser.add_argument("--top", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not (0.0 < args.rate <= 1.0):
        raise SystemExit("--rate must lie in (0,1]")
    if not (0.0 < args.delta < 1.0):
        raise SystemExit("--delta must lie in (0,1)")
    length = (
        int(round(args.k * (1.0 / args.rate - 1.0)))
        if args.systematic
        else int(round(args.k / args.rate))
    )
    if length <= 0:
        raise SystemExit("the selected rate leaves no nonsystematic output coordinates")
    if not (0.0 < args.row_weight < 0.5 * length):
        raise SystemExit("--row-weight must induce a Bernoulli density in (0,1/2)")
    if args.code == "ec" and args.memory < 1:
        raise SystemExit("EC requires positive --memory")

    weights = scan_weights(args.k, args.dense_prefix, args.geometric_points)
    results = evaluate(args, weights)
    results.sort(reverse=True)

    print(
        f"code={args.code} k={args.k} n={length} rate={args.rate:g} "
        f"row_weight={args.row_weight:g} memory={args.memory} delta={args.delta:g} "
        f"wrapping={args.wrapping} systematic={args.systematic}"
    )
    print(f"scanned_message_weights={len(weights)}")
    print("log(first-moment term) / log2(term), message weight, activation q, Chernoff z")
    for log_term, message_weight, q, z in results[: args.top]:
        print(
            f"{log_term: .12g} / {log_term / math.log(2.0): .12g}, "
            f"r={message_weight}, q={q:.12g}, z={z:.12g}"
        )


if __name__ == "__main__":
    main()
