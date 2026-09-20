#!/usr/bin/env python3
"""Diagnostic landscape for parity-free random outer constituents.

The model fixes a rate-half random linear outer code independently at every
block, a global permutation of g-bit packets, and the one-lap RandomStepConv
inner.  It estimates the first moment over low-weight nonzero codewords for a
grid of outer block lengths and inner state sizes.

This is a landscape tool, not a certificate.  Outer coefficients use a
one-dimensional lattice saddle approximation.  The support sum uses adaptive
sampling followed by log-linear interpolation.  Boundary points selected for
a proof must later receive exact coefficients, a complete support cover, and
outward rounding.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp

from analyze_riffle_randomstepconv_g4_sigma20_goal01 import (
    optimize_tilts,
    transfer_coefficients,
)


LOG2 = math.log(2.0)
DEFAULT_MESSAGE_BITS = 1 << 20
DEFAULT_DISTANCE = 0.09
DEFAULT_OUTPUT = Path(
    "constructions/riffle_random_outer_landscape/receipts/"
    "g4_initial_slice.json"
)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


def log_expm1(value: float) -> float:
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def log_geometric_sum(first: float, delta: float, count: int) -> float:
    """Return log(sum(exp(first + j*delta), j=0..count-1))."""
    if count <= 0:
        return -math.inf
    if count == 1:
        return first
    if abs(delta) < 1e-10:
        return first + math.log(count)
    if delta < 0.0:
        return (
            first
            + math.log1p(-math.exp(count * delta))
            - math.log1p(-math.exp(delta))
        )
    return (
        first
        + (count - 1) * delta
        + math.log1p(-math.exp(-count * delta))
        - math.log1p(-math.exp(-delta))
    )


@dataclass(frozen=True)
class OuterPoint:
    support: int
    log_t: float
    coefficient_saddle_log2: float
    coefficient_chernoff_log2: float
    tilted_variance: float
    mean_active_blocks: float
    mean_support_per_active_block: float
    coefficient_method: str
    occupation_cutoff: int | None


class RandomOuterMoment:
    """Expected support enumerator for independent random rate-half blocks."""

    def __init__(self, message_bits: int, outer_bits: int, packet_bits: int) -> None:
        if outer_bits <= 0 or outer_bits % 2:
            raise ValueError("outer length must be positive and even")
        if packet_bits <= 0 or outer_bits % packet_bits:
            raise ValueError("packet size must divide the outer length")
        input_bits = outer_bits // 2
        if message_bits % input_bits:
            raise ValueError("outer dimension must divide the message length")
        self.message_bits = message_bits
        self.outer_bits = outer_bits
        self.input_bits = input_bits
        self.packet_bits = packet_bits
        self.packet_values = (1 << packet_bits) - 1
        self.packets_per_block = outer_bits // packet_bits
        self.blocks = message_bits // input_bits
        self.packet_positions = 2 * message_bits // packet_bits
        self.log_nonzero_ratio = (
            log_two_power_minus_one(input_bits)
            - log_two_power_minus_one(outer_bits)
        )

    def tilted_stats(
        self, log_t: float
    ) -> tuple[float, float, float, float, float]:
        """Return log Z, mean, variance, active probability, active mean."""
        log_at = math.log(self.packet_values) + log_t
        if log_at > 40.0:
            log_u = log_at + math.log1p(math.exp(-log_at))
            packet_probability = 1.0 / (1.0 + math.exp(-log_at))
        else:
            at = math.exp(log_at)
            log_u = math.log1p(at)
            packet_probability = at / (1.0 + at)

        exponent = self.packets_per_block * log_u
        log_active_mass = self.log_nonzero_ratio + log_expm1(exponent)
        log_z = float(np.logaddexp(0.0, log_active_mass))
        active_probability = math.exp(log_active_mass - log_z)

        nonzero_probability = -math.expm1(-exponent)
        active_mean = (
            self.packets_per_block
            * packet_probability
            / nonzero_probability
        )
        active_factorial_second = (
            self.packets_per_block
            * (self.packets_per_block - 1)
            * packet_probability
            * packet_probability
            / nonzero_probability
        )
        block_mean = active_probability * active_mean
        block_factorial_second = (
            active_probability * active_factorial_second
        )
        block_variance = max(
            np.finfo(float).tiny,
            block_factorial_second + block_mean - block_mean * block_mean,
        )
        return (
            log_z,
            self.blocks * block_mean,
            self.blocks * block_variance,
            active_probability,
            active_mean,
        )

    @lru_cache(maxsize=None)
    def coefficient_point(self, support: int) -> OuterPoint:
        if not 0 < support < self.packet_positions:
            raise ValueError("support must lie strictly inside packet range")

        def centered(log_t: float) -> float:
            return self.tilted_stats(log_t)[1] - support

        low = -80.0
        high = 40.0
        if centered(low) >= 0.0 or centered(high) <= 0.0:
            raise ArithmeticError("failed to bracket outer coefficient saddle")
        log_t = float(brentq(centered, low, high, xtol=2e-12, rtol=2e-14))
        log_z, _mean, variance, active_probability, active_mean = (
            self.tilted_stats(log_t)
        )
        chernoff = self.blocks * log_z - support * log_t
        saddle = chernoff - 0.5 * math.log(2.0 * math.pi * variance)
        method = "global_lattice_saddle"
        occupation_cutoff: int | None = None
        mean_active_blocks = self.blocks * active_probability
        # The global lattice saddle is unreliable when the tilted law is
        # mostly the all-zero message.  Once the mean occupation is a few
        # blocks, its error is small enough for a landscape diagnostic.  Keep
        # the exact kernel narrowly on the genuinely sparse regime.
        if mean_active_blocks < 2.5:
            occupation_cutoff = min(
                self.blocks,
                support,
                max(
                    32,
                    int(
                        math.ceil(
                            mean_active_blocks
                            + 16.0 * math.sqrt(mean_active_blocks + 1.0)
                        )
                    ),
                ),
            )
            minimum_occupation = max(
                1, math.ceil(support / self.packets_per_block)
            )
            contributions: list[float] = []
            weighted_occupations: list[int] = []
            for occupation in range(
                minimum_occupation, occupation_cutoff + 1
            ):
                subset_count = self.nonempty_subset_count(
                    occupation, support
                )
                if subset_count == 0:
                    continue
                contribution = (
                    math.lgamma(self.blocks + 1)
                    - math.lgamma(occupation + 1)
                    - math.lgamma(self.blocks - occupation + 1)
                    + occupation * self.log_nonzero_ratio
                    + math.log(subset_count)
                    + support * math.log(self.packet_values)
                )
                contributions.append(contribution)
                weighted_occupations.append(occupation)
            if contributions:
                saddle = float(logsumexp(contributions))
                probabilities = np.exp(
                    np.asarray(contributions) - saddle
                )
                mean_active_blocks = float(
                    probabilities @ np.asarray(weighted_occupations)
                )
                active_mean = support / mean_active_blocks
                method = "exact_occupation_sum_truncated"
        return OuterPoint(
            support=support,
            log_t=log_t,
            coefficient_saddle_log2=saddle / LOG2,
            coefficient_chernoff_log2=chernoff / LOG2,
            tilted_variance=variance,
            mean_active_blocks=mean_active_blocks,
            mean_support_per_active_block=active_mean,
            coefficient_method=method,
            occupation_cutoff=occupation_cutoff,
        )

    @lru_cache(maxsize=None)
    def nonempty_subset_count(self, occupation: int, support: int) -> int:
        """Count support subsets meeting every one of the active blocks."""
        total = 0
        block_packets = self.packets_per_block
        for omitted in range(occupation + 1):
            available = (occupation - omitted) * block_packets
            if available < support:
                continue
            term = math.comb(occupation, omitted) * math.comb(
                available, support
            )
            total = total - term if omitted & 1 else total + term
        if total < 0:
            raise ArithmeticError("nonempty subset inclusion-exclusion changed sign")
        return total


def initial_support_schedule(packet_positions: int, maximum: int) -> list[int]:
    stop = min(maximum, packet_positions - 1)
    supports = set(range(1, min(256, stop) + 1))
    supports.update(range(272, min(4096, stop) + 1, 16))
    current = 4352
    while current <= stop:
        supports.add(current)
        current = max(current + 256, int(math.ceil(current * 1.12 / 16.0)) * 16)
    supports.add(stop)
    return sorted(supports)


@lru_cache(maxsize=None)
def inner_point_episode(
    packet_positions: int,
    packet_bits: int,
    sigma: int,
    support: int,
    distance: int,
) -> tuple[float, int]:
    result = optimize_tilts(
        packet_positions, packet_bits, sigma, support, distance
    )
    return (
        float(result["diagnostic_log2_upper"]),
        int(result["dominant_termination_count"]),
    )


def log_matrix_power_moment(matrix: np.ndarray, exponent: int) -> float:
    """Return log([1,0] matrix^exponent [1,1]^T) with scaling."""
    power = np.asarray(matrix, dtype=np.float64)
    power_scale = 0.0
    vector = np.asarray([1.0, 0.0], dtype=np.float64)
    vector_scale = 0.0
    remaining = exponent
    while remaining:
        if remaining & 1:
            vector = vector @ power
            norm = float(np.max(vector))
            if not math.isfinite(norm) or norm <= 0.0:
                raise ArithmeticError("matrix moment vector lost positive scale")
            vector /= norm
            vector_scale += power_scale + math.log(norm)
        remaining >>= 1
        if not remaining:
            break
        power = power @ power
        norm = float(np.max(power))
        if not math.isfinite(norm) or norm <= 0.0:
            raise ArithmeticError("matrix moment power lost positive scale")
        power /= norm
        power_scale = 2.0 * power_scale + math.log(norm)
    return vector_scale + math.log(float(vector.sum()))


@lru_cache(maxsize=None)
def inner_point_matrix(
    packet_positions: int,
    packet_bits: int,
    sigma: int,
    support: int,
    distance: int,
) -> tuple[float, float, float]:
    """Joint support/output Chernoff bound from the exact two-state moment."""
    q = math.ldexp(1.0, -sigma)
    log_choose = (
        math.lgamma(packet_positions + 1)
        - math.lgamma(support + 1)
        - math.lgamma(packet_positions - support + 1)
    )

    def objective(parameters: np.ndarray) -> float:
        log_y = float(parameters[0])
        log_surprisal = float(parameters[1])
        y = math.exp(log_y)
        z = math.exp(-math.exp(log_surprisal))
        b = math.ldexp((1.0 + z) ** packet_bits, -packet_bits)
        qb = q * b
        live_b = (1.0 - q) * b
        matrix = np.asarray(
            [
                [1.0 + y * qb, y * live_b],
                [(1.0 + y) * qb, (1.0 + y) * live_b],
            ],
            dtype=np.float64,
        )
        log_moment = log_matrix_power_moment(matrix, packet_positions)
        return (
            log_moment
            - support * log_y
            - log_choose
            + distance * math.exp(log_surprisal)
        )

    base_log_y = math.log(support / (packet_positions - support))
    starts = (
        np.asarray([base_log_y, -4.0]),
        np.asarray([base_log_y, -1.0]),
        np.asarray([base_log_y, 1.0]),
    )
    results = [
        minimize(
            objective,
            start,
            method="Nelder-Mead",
            bounds=((-40.0, 40.0), (-18.0, 3.0)),
            options={"xatol": 2e-8, "fatol": 2e-8, "maxiter": 500},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: float(result.fun))
    value = min(0.0, float(best.fun) / LOG2)
    return value, float(best.x[0]), float(best.x[1])


def interpolated_logsum2(points: list[tuple[int, float]]) -> float:
    """Log2 sum over integer supports under log-linear interpolation."""
    ordered = sorted(points)
    natural_segments: list[float] = []
    for index in range(len(ordered) - 1):
        left_h, left_value = ordered[index]
        right_h, right_value = ordered[index + 1]
        width = right_h - left_h
        if width <= 0:
            continue
        first = left_value * LOG2
        delta = (right_value - left_value) * LOG2 / width
        natural_segments.append(log_geometric_sum(first, delta, width))
    natural_segments.append(ordered[-1][1] * LOG2)
    return float(logsumexp(natural_segments)) / LOG2


def local_maxima(rows: list[dict[str, object]]) -> list[int]:
    ordered = sorted(rows, key=lambda row: int(row["packet_support"]))
    maxima: list[int] = []
    for index, row in enumerate(ordered):
        value = float(row["pointwise_saddle_log2"])
        left = (
            float(ordered[index - 1]["pointwise_saddle_log2"])
            if index
            else -math.inf
        )
        right = (
            float(ordered[index + 1]["pointwise_saddle_log2"])
            if index + 1 < len(ordered)
            else -math.inf
        )
        if value >= left and value >= right:
            maxima.append(int(row["packet_support"]))
    maxima.sort(
        key=lambda support: float(
            next(
                row["pointwise_saddle_log2"]
                for row in ordered
                if int(row["packet_support"]) == support
            )
        ),
        reverse=True,
    )
    return maxima


def evaluate_sigma(
    outer: RandomOuterMoment,
    sigma: int,
    distance: int,
    coarse_supports: Iterable[int],
    refinement_radius: int,
    inner_method: str,
) -> dict[str, object]:
    rows_by_support: dict[int, dict[str, object]] = {}

    def evaluate(support: int) -> None:
        if support in rows_by_support:
            return
        outer_point = outer.coefficient_point(support)
        if inner_method == "matrix":
            inner_log2, inner_log_y, inner_log_surprisal = inner_point_matrix(
                outer.packet_positions,
                outer.packet_bits,
                sigma,
                support,
                distance,
            )
            terminations: int | None = None
        else:
            inner_log2, terminations = inner_point_episode(
                outer.packet_positions,
                outer.packet_bits,
                sigma,
                support,
                distance,
            )
            inner_log_y = None
            inner_log_surprisal = None
        rows_by_support[support] = {
            "packet_support": support,
            "outer_saddle_log2": outer_point.coefficient_saddle_log2,
            "outer_chernoff_log2": outer_point.coefficient_chernoff_log2,
            "outer_log_tilt": outer_point.log_t,
            "outer_tilted_variance": outer_point.tilted_variance,
            "inner_log2_upper": inner_log2,
            "inner_method": inner_method,
            "inner_log_support_tilt": inner_log_y,
            "inner_log_output_surprisal": inner_log_surprisal,
            "pointwise_saddle_log2": (
                outer_point.coefficient_saddle_log2 + inner_log2
            ),
            "pointwise_chernoff_log2": (
                outer_point.coefficient_chernoff_log2 + inner_log2
            ),
            "dominant_termination_count": terminations,
            "mean_active_outer_blocks": outer_point.mean_active_blocks,
            "mean_packet_support_per_active_block": (
                outer_point.mean_support_per_active_block
            ),
            "outer_coefficient_method": outer_point.coefficient_method,
            "outer_occupation_cutoff": outer_point.occupation_cutoff,
        }

    for support in coarse_supports:
        evaluate(support)

    coarse_rows = list(rows_by_support.values())
    peaks = local_maxima(coarse_rows)[:4]
    for peak in peaks:
        nearby = sorted(
            int(row["packet_support"])
            for row in coarse_rows
            if abs(int(row["packet_support"]) - peak) <= refinement_radius
        )
        left = max(1, peak - refinement_radius)
        right = min(outer.packet_positions - 1, peak + refinement_radius)
        if nearby:
            left = min(left, nearby[0])
            right = max(right, nearby[-1])
        for support in range(left, right + 1):
            evaluate(support)

    rows = sorted(rows_by_support.values(), key=lambda row: int(row["packet_support"]))
    dominant = max(rows, key=lambda row: float(row["pointwise_saddle_log2"]))
    if inner_method == "matrix":
        _episode_log2, episode_terminations = inner_point_episode(
            outer.packet_positions,
            outer.packet_bits,
            sigma,
            int(dominant["packet_support"]),
            distance,
        )
        dominant["episode_bound_dominant_termination_count"] = (
            episode_terminations
        )
    saddle_total = interpolated_logsum2(
        [
            (int(row["packet_support"]), float(row["pointwise_saddle_log2"]))
            for row in rows
        ]
    )
    chernoff_total = interpolated_logsum2(
        [
            (int(row["packet_support"]), float(row["pointwise_chernoff_log2"]))
            for row in rows
        ]
    )
    return {
        "sigma": sigma,
        "estimated_lambda_bits": -saddle_total,
        "interpolated_chernoff_margin_bits": -chernoff_total,
        "estimated_log2_first_moment": saddle_total,
        "interpolated_chernoff_log2_first_moment": chernoff_total,
        "dominant_profile": dominant,
        "evaluated_support_count": len(rows),
        "evaluated_support_min": int(rows[0]["packet_support"]),
        "evaluated_support_max": int(rows[-1]["packet_support"]),
        "support_rows": rows,
    }


def exact_outer_small_check() -> dict[str, object]:
    """Compare the exact generating polynomial with direct message counting."""
    message_bits = 8
    outer_bits = 8
    packet_bits = 2
    outer = RandomOuterMoment(message_bits, outer_bits, packet_bits)
    ratio = math.exp(outer.log_nonzero_ratio)
    block = np.zeros(outer.packets_per_block + 1)
    block[0] = 1.0
    for support in range(1, outer.packets_per_block + 1):
        block[support] = (
            ratio
            * math.comb(outer.packets_per_block, support)
            * outer.packet_values**support
        )
    polynomial = np.asarray([1.0])
    for _index in range(outer.blocks):
        polynomial = np.convolve(polynomial, block)
    expected_mass = float(1 << message_bits)
    error = abs(float(polynomial.sum()) - expected_mass)
    if error > 1e-10:
        raise AssertionError("random-outer generating polynomial mass changed")
    return {
        "message_bits": message_bits,
        "outer_bits": outer_bits,
        "packet_bits": packet_bits,
        "blocks": outer.blocks,
        "polynomial_mass": float(polynomial.sum()),
        "expected_message_count": expected_mass,
        "absolute_error": error,
        "status": "EXACT_FORMULA_MATCH",
    }


def exact_inner_small_checks() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for packet_bits, sigma, packet_positions, support, distance in (
        (1, 1, 3, 1, 0),
        (1, 2, 4, 2, 0),
        (2, 1, 3, 1, 0),
    ):
        coefficients = transfer_coefficients(
            packet_bits, sigma, packet_positions
        )
        exact = sum(
            mass
            for (row_support, weight), mass in coefficients.items()
            if row_support == support and weight <= distance
        ) / math.comb(packet_positions, support)
        bound_log2, log_y, log_surprisal = inner_point_matrix(
            packet_positions,
            packet_bits,
            sigma,
            support,
            distance,
        )
        bound = math.exp(bound_log2 * LOG2)
        exact_float = float(exact)
        if bound + 2e-10 < exact_float:
            raise AssertionError("matrix Chernoff value fell below exact probability")
        rows.append(
            {
                "packet_bits": packet_bits,
                "sigma": sigma,
                "packet_positions": packet_positions,
                "support": support,
                "distance": distance,
                "exact_probability": exact_float,
                "matrix_upper_probability": bound,
                "matrix_log_support_tilt": log_y,
                "matrix_log_output_surprisal": log_surprisal,
                "status": "UPPER_BOUND_CONFIRMED",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--g", type=int, choices=(1, 2, 4, 8), default=4)
    parser.add_argument("--outer-bits", nargs="+", type=int, default=[512])
    parser.add_argument("--sigmas", nargs="+", type=int, default=[16, 18, 20])
    parser.add_argument("--message-bits", type=int, default=DEFAULT_MESSAGE_BITS)
    parser.add_argument("--relative-distance", type=float, default=DEFAULT_DISTANCE)
    parser.add_argument("--max-support", type=int, default=32768)
    parser.add_argument("--refinement-radius", type=int, default=96)
    parser.add_argument(
        "--inner-method",
        choices=("matrix", "episode"),
        default="matrix",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    check = exact_outer_small_check()
    inner_checks = exact_inner_small_checks()
    if args.self_test_only:
        print(
            json.dumps(
                {"outer": check, "inner": inner_checks},
                indent=2,
                sort_keys=True,
            )
        )
        return

    cases: list[dict[str, object]] = []
    for outer_bits in args.outer_bits:
        outer = RandomOuterMoment(args.message_bits, outer_bits, args.g)
        distance = math.floor(
            args.relative_distance * args.g * outer.packet_positions
        )
        schedule = initial_support_schedule(
            outer.packet_positions, args.max_support
        )
        sigma_rows: list[dict[str, object]] = []
        for sigma in args.sigmas:
            print(
                f"case,g,{args.g},B,{outer_bits},sigma,{sigma},"
                f"coarse_supports,{len(schedule)}",
                flush=True,
            )
            result = evaluate_sigma(
                outer,
                sigma,
                distance,
                schedule,
                args.refinement_radius,
                args.inner_method,
            )
            sigma_rows.append(result)
            dominant = result["dominant_profile"]
            print(
                f"result,g,{args.g},B,{outer_bits},sigma,{sigma},"
                f"lambda,{result['estimated_lambda_bits']:.6f},"
                f"support,{dominant['packet_support']},"
                f"pointwise,{dominant['pointwise_saddle_log2']:.6f},"
                f"terminations,{dominant.get('episode_bound_dominant_termination_count', dominant.get('dominant_termination_count'))}",
                flush=True,
            )
        cases.append(
            {
                "outer_bits": outer_bits,
                "outer_dimension_bits": outer.input_bits,
                "outer_blocks": outer.blocks,
                "packets_per_outer_block": outer.packets_per_block,
                "packet_positions": outer.packet_positions,
                "distance": distance,
                "sigma_rows": sigma_rows,
            }
        )

    payload = {
        "schema": "riffle-random-outer-landscape-diagnostic-v1",
        "model": (
            "independent random rate-half outer blocks, no extra parity, "
            "global packet permutation, one-lap RandomStepConv"
        ),
        "message_bits": args.message_bits,
        "packet_bits": args.g,
        "relative_distance": args.relative_distance,
        "target_lambda_bits": 40,
        "inner_method": args.inner_method,
        "outer_small_exact_check": check,
        "inner_small_exact_checks": inner_checks,
        "cases": cases,
        "scope": (
            "Outer coefficients use a lattice saddle approximation. The "
            "support total uses adaptive sampling and log-linear interpolation. "
            "This receipt maps the landscape and is not a proof certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
