#!/usr/bin/env python3
"""Build a compressed all-support envelope for RandomStepConv Goal 02."""

from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from analyze_riffle_randomstepconv_g4_sigma20_goal01 import (
    TARGET_BITS,
    TARGET_G,
    TARGET_PACKETS,
    TARGET_SIGMA,
    decimal_context,
    log_binomial,
)


DEFAULT_DELTAS = (0.05, 0.09, 0.12)
INTERVAL_STARTS = (
    18,
    24,
    32,
    48,
    64,
    92,
    128,
    256,
    512,
    1024,
    2048,
    4096,
    8192,
    16384,
    32768,
    65536,
    131072,
    262144,
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_randomstepconv_g4_sigma20/receipts/"
    "goal02_all_support_envelope.json"
)
LOG2 = math.log(2.0)
# These bounds stay strictly inside the representable open interval (0,1)
# after the logistic transformations in decode_tilts.
OPTIMIZER_BOUNDS = ((-16.0, 4.0), (-35.0, 35.0))


def support_intervals() -> tuple[tuple[int, int], ...]:
    ends = tuple(start - 1 for start in INTERVAL_STARTS[1:]) + (TARGET_PACKETS,)
    intervals = tuple(zip(INTERVAL_STARTS, ends))
    if intervals[0][0] != 18 or intervals[-1][1] != TARGET_PACKETS:
        raise AssertionError("support intervals do not cover the required endpoints")
    for index, (left, right) in enumerate(intervals):
        if left > right:
            raise AssertionError("empty support interval")
        if index and intervals[index - 1][1] + 1 != left:
            raise AssertionError("support intervals are not contiguous")
    return intervals


def exact_discrete_convexity_check() -> dict[str, object]:
    # Consecutive slopes of -log(C(N,h)) differ by the logarithm of a
    # ratio whose cross-product difference is exactly N+1.
    difference = TARGET_PACKETS + 1
    return {
        "support_range": [0, TARGET_PACKETS - 2],
        "cross_product_difference": difference,
        "strictly_positive": difference > 0,
        "status": "EXACT_INTEGER_CHECK",
    }


def decode_tilts(parameters: np.ndarray | tuple[float, float]) -> tuple[float, float]:
    weight_tilt = math.exp(-math.exp(float(parameters[0])))
    radius = 1.0 / (1.0 + math.exp(-float(parameters[1])))
    return weight_tilt, radius


def log_common_bound(
    n: int,
    g: int,
    sigma: int,
    support: int,
    distance: int,
    weight_tilt: float,
    radius: float,
) -> float:
    if not 0 <= support <= n:
        raise ValueError("support must lie in [0,n]")
    if not 0.0 < weight_tilt < 1.0 or not 0.0 < radius < 1.0:
        raise ValueError("tilts must lie in (0,1)")
    q = math.ldexp(1.0, -sigma)
    b = math.ldexp((1.0 + weight_tilt) ** g, -g)
    d = (1.0 - q) * b
    gap_generating_value = (b - d * radius) / (
        (1.0 - radius) * (1.0 - d * radius)
    )
    return (
        -distance * math.log(weight_tilt)
        - log_binomial(n, support)
        - (n - support) * math.log(radius)
        - math.log1p(-radius)
        + support * math.log(gap_generating_value)
    )


def initial_parameters(support: int) -> np.ndarray:
    output_parameter = min(3.0, max(-18.0, math.log(max(1, support)) - 12.0))
    radius_parameter = math.log((TARGET_PACKETS - support + 1.0) / (support + 1.0)) + 1.5
    return np.asarray((output_parameter, min(35.0, max(-35.0, radius_parameter))))


def optimize_point(support: int, distance: int) -> dict[str, object]:
    def objective(parameters: np.ndarray) -> float:
        z, r = decode_tilts(parameters)
        return log_common_bound(
            TARGET_PACKETS, TARGET_G, TARGET_SIGMA, support, distance, z, r
        )

    starts = (
        initial_parameters(support),
        np.asarray((-8.0, 10.0)),
        np.asarray((-2.0, 0.0)),
    )
    results = [
        minimize(
            objective,
            start,
            method="Powell",
            bounds=OPTIMIZER_BOUNDS,
            options={"xtol": 1e-10, "ftol": 1e-12, "maxiter": 1000},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: float(result.fun))
    z, r = decode_tilts(best.x)
    return {
        "parameters": np.asarray(best.x),
        "weight_tilt": z,
        "coefficient_radius": r,
        "diagnostic_log2_upper": min(0.0, float(best.fun) / LOG2),
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
    }


def optimize_interval(left: int, right: int, distance: int) -> dict[str, object]:
    left_point = optimize_point(left, distance)
    right_point = optimize_point(right, distance)

    def endpoint_values(parameters: np.ndarray) -> tuple[float, float]:
        z, r = decode_tilts(parameters)
        return (
            log_common_bound(
                TARGET_PACKETS, TARGET_G, TARGET_SIGMA, left, distance, z, r
            ),
            log_common_bound(
                TARGET_PACKETS, TARGET_G, TARGET_SIGMA, right, distance, z, r
            ),
        )

    def objective(parameters: np.ndarray) -> float:
        return max(endpoint_values(parameters))

    candidates = [left_point["parameters"], right_point["parameters"]]
    midpoint = (np.asarray(candidates[0]) + np.asarray(candidates[1])) / 2.0
    candidates.append(midpoint)
    optimized = [
        minimize(
            objective,
            np.asarray(start),
            method="Powell",
            bounds=OPTIMIZER_BOUNDS,
            options={"xtol": 1e-10, "ftol": 1e-12, "maxiter": 1500},
        )
        for start in candidates
    ]
    all_candidates = [(np.asarray(start), None) for start in candidates]
    all_candidates.extend((np.asarray(result.x), result) for result in optimized)
    best_parameters, best_result = min(
        all_candidates, key=lambda item: objective(item[0])
    )
    z, r = decode_tilts(best_parameters)
    left_log, right_log = endpoint_values(best_parameters)
    middle = (left + right) // 2
    middle_log = log_common_bound(
        TARGET_PACKETS, TARGET_G, TARGET_SIGMA, middle, distance, z, r
    )
    point_middle = optimize_point(middle, distance)
    return {
        "weight_tilt": z,
        "coefficient_radius": r,
        "left_diagnostic_log2": min(0.0, left_log / LOG2),
        "middle_diagnostic_log2": min(0.0, middle_log / LOG2),
        "right_diagnostic_log2": min(0.0, right_log / LOG2),
        "middle_pointwise_diagnostic_log2": point_middle["diagnostic_log2_upper"],
        "middle_slack_bits": max(
            0.0,
            (middle_log / LOG2) - float(point_middle["diagnostic_log2_upper"]),
        ),
        "optimizer_success": best_result is None or bool(best_result.success),
        "optimizer_message": (
            "endpoint point certificate"
            if best_result is None
            else str(best_result.message)
        ),
    }


@lru_cache(maxsize=None)
def exact_binomial_decimal(n: int, support: int) -> Decimal:
    return Decimal(math.comb(n, support))


def outward_common_bound(
    support: int,
    distance: int,
    weight_tilt: float,
    radius: float,
) -> dict[str, str]:
    z = Decimal(format(weight_tilt, ".17g"))
    r = Decimal(format(radius, ".17g"))
    one = Decimal(1)
    q = one / Decimal(1 << TARGET_SIGMA)
    floor_context = decimal_context(ROUND_FLOOR)
    ceil_context = decimal_context(ROUND_CEILING)

    with localcontext(ceil_context):
        b_upper = (one + z) ** TARGET_G / Decimal(1 << TARGET_G)
        gap_numerator_upper = b_upper * (one - (one - q) * r)
        numerator_upper = gap_numerator_upper**support
    with localcontext(ceil_context):
        d_upper = (one - q) * b_upper
    with localcontext(floor_context):
        denominator_lower = (
            exact_binomial_decimal(TARGET_PACKETS, support)
            * (z**distance)
            * (r ** (TARGET_PACKETS - support))
            * ((one - r) ** (support + 1))
            * ((one - d_upper * r) ** support)
        )
    with localcontext(ceil_context):
        ratio_upper = numerator_upper / denominator_lower
        if ratio_upper >= one:
            ratio_upper = one
            log2_upper = Decimal(0)
        else:
            log2_upper = ratio_upper.ln(context=ceil_context) / Decimal(2).ln(
                context=ceil_context
            )
    return {
        "verified_weight_tilt_decimal": str(z),
        "verified_coefficient_radius_decimal": str(r),
        "outward_probability_upper": str(ratio_upper),
        "outward_log2_upper": str(log2_upper),
    }


def convexity_check(
    left: int, right: int, distance: int, weight_tilt: float, radius: float
) -> dict[str, object]:
    probes = sorted({left, min(right, left + 1), (left + right) // 2, max(left, right - 1), right})
    values = [
        log_common_bound(
            TARGET_PACKETS,
            TARGET_G,
            TARGET_SIGMA,
            support,
            distance,
            weight_tilt,
            radius,
        )
        for support in probes
    ]
    cap = max(values[0], values[-1])
    return {
        "probe_supports": probes,
        "probe_log2_bounds": [min(0.0, value / LOG2) for value in values],
        "all_probes_below_endpoint_cap": all(value <= cap + 1e-8 for value in values),
    }


def build_rows(deltas: tuple[float, ...]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for delta in deltas:
        distance = math.floor(delta * TARGET_BITS)
        for left, right in support_intervals():
            optimized = optimize_interval(left, right, distance)
            left_verified = outward_common_bound(
                left,
                distance,
                float(optimized["weight_tilt"]),
                float(optimized["coefficient_radius"]),
            )
            right_verified = outward_common_bound(
                right,
                distance,
                float(optimized["weight_tilt"]),
                float(optimized["coefficient_radius"]),
            )
            cap = max(
                Decimal(left_verified["outward_log2_upper"]),
                Decimal(right_verified["outward_log2_upper"]),
            )
            convexity = convexity_check(
                left,
                right,
                distance,
                float(optimized["weight_tilt"]),
                float(optimized["coefficient_radius"]),
            )
            if not convexity["all_probes_below_endpoint_cap"]:
                raise AssertionError("floating convexity probe failed")
            row: dict[str, object] = {
                "relative_binary_weight": delta,
                "distance": distance,
                "support_left": left,
                "support_right": right,
                "outward_interval_cap_log2": str(cap),
                "left_outward": left_verified,
                "right_outward": right_verified,
                "convexity_probe": convexity,
            }
            row.update(optimized)
            rows.append(row)
            print(
                f"delta,{delta:.6f},interval,{left}:{right},"
                f"cap_log2,{cap},mid_slack_bits,{optimized['middle_slack_bits']:.6f},"
                f"z,{left_verified['verified_weight_tilt_decimal']},"
                f"r,{left_verified['verified_coefficient_radius_decimal']}",
                flush=True,
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--relative-distances", nargs="+", type=float, default=list(DEFAULT_DELTAS)
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    deltas = tuple(args.relative_distances)
    rows = build_rows(deltas)
    intervals = support_intervals()
    payload = {
        "schema": "riffle-randomstepconv-goal02-v1",
        "construction": "Riffle RandomStepConv g=4 sigma=20",
        "evidence": {
            "fixed_support_identity": "EXACT_FROM_GOAL01",
            "interval_endpoint_reduction": "EXACT_BY_BINOMIAL_LOG_CONCAVITY",
            "optimizer": "FLOATING_DIAGNOSTIC",
            "endpoint_evaluation": "DECIMAL_OUTWARD_ROUNDED",
        },
        "parameters": {
            "g": TARGET_G,
            "sigma": TARGET_SIGMA,
            "packet_positions": TARGET_PACKETS,
            "binary_output_length": TARGET_BITS,
            "terminal_state": "discarded",
        },
        "coverage": {
            "minimum_support": intervals[0][0],
            "maximum_support": intervals[-1][1],
            "interval_count_per_distance": len(intervals),
            "contiguous": True,
        },
        "discrete_convexity_check": exact_discrete_convexity_check(),
        "rows": rows,
        "scope": (
            "The receipt bounds only the RandomStepConv inner for every fixed "
            "input-packet support from 18 through N. It does not consume the outer "
            "spectrum or certify an end-to-end distance."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
