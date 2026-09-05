#!/usr/bin/env python3
"""Exact checks and floating diagnostics for two-lap RandomStepConv Goal 01.

The first lap is an unretained burn-in. The second lap uses independent maps
and emits the retained output. The target calculation evaluates the exact gap
moment through a two-variable coefficient tilt; it does not expand any target-
size polynomial.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.stats import binom


TARGET_G = 4
TARGET_SIGMA = 20
TARGET_PACKETS = 524_352
TARGET_BITS = TARGET_G * TARGET_PACKETS
DEFAULT_SUPPORTS = (18, 24, 32, 64, 256, 1024, 16_384, 20_000, 24_000, 28_000, 32_767)
DEFAULT_DELTAS = (0.09,)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_randomstepconv_2lap_g4_sigma20/receipts/"
    "goal01_cyclic_live_coverage.json"
)


def log_binomial(n: int, h: int) -> float:
    return float(gammaln(n + 1) - gammaln(h + 1) - gammaln(n - h + 1))


def direct_support_moment(
    g: int, sigma: int, support_mask: int, packet_count: int, z: Fraction
) -> Fraction:
    """Evaluate both laps directly in the zero/live state quotient."""
    q = Fraction(1, 1 << sigma)
    a = 1 - q
    b = (1 + z) ** g / (1 << g)
    state = (Fraction(1), Fraction(0))
    for position in range(packet_count):
        active = bool((support_mask >> position) & 1)
        zero, live = state
        if active:
            mass = zero + live
            state = (mass * q, mass * a)
        else:
            state = (zero + live * q, live * a)

    for position in range(packet_count):
        active = bool((support_mask >> position) & 1)
        zero, live = state
        if active:
            mass = (zero + live) * b
            state = (mass * q, mass * a)
        else:
            state = (zero + live * b * q, live * b * a)
    return state[0] + state[1]


def gap_a(length: int, q: Fraction, u: Fraction) -> Fraction:
    a = 1 - q
    value = a**length
    for live_prefix in range(length):
        value += a**live_prefix * q * u ** (length - live_prefix)
    return value


def gap_b(length: int, q: Fraction, u: Fraction) -> Fraction:
    """Off-count moment for a zero run entered in a nonzero state."""
    if length == 0:
        return Fraction(1)
    a = 1 - q
    value = a ** (length - 1)
    for live_count in range(1, length):
        value += a ** (live_count - 1) * q * u ** (length - live_count)
    return value


def weak_compositions(total: int, parts: int) -> Iterable[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first,) + tail


def gap_formula_moment(
    g: int, sigma: int, packet_count: int, support: int, z: Fraction
) -> Fraction:
    if support == 0:
        return Fraction(1)
    q = Fraction(1, 1 << sigma)
    a = 1 - q
    b = (1 + z) ** g / (1 << g)
    u = 1 / b
    remaining = packet_count - support
    total = Fraction(0)
    # Composition order: prefix, internal gaps, suffix.
    for composition in weak_compositions(remaining, support + 1):
        prefix = composition[0]
        suffix = composition[-1]
        internal = composition[1:-1]
        prefix_moment = (
            a ** (suffix + 1) * gap_b(prefix, q, u)
            + (1 - a ** (suffix + 1)) * u**prefix
        )
        term = gap_a(suffix, q, u) * prefix_moment
        for length in internal:
            term *= gap_a(length, q, u)
        total += term
    return b**packet_count * total / math.comb(packet_count, support)


def exact_checks() -> list[dict[str, object]]:
    cases = ((1, 1, 4), (2, 1, 3), (1, 2, 4))
    tilts = (Fraction(1, 2), Fraction(2, 3))
    rows: list[dict[str, object]] = []
    for g, sigma, packet_count in cases:
        for z in tilts:
            for support in range(1, packet_count + 1):
                direct = sum(
                    direct_support_moment(g, sigma, mask, packet_count, z)
                    for mask in range(1 << packet_count)
                    if mask.bit_count() == support
                ) / math.comb(packet_count, support)
                formula = gap_formula_moment(g, sigma, packet_count, support, z)
                if direct != formula:
                    raise AssertionError(
                        "gap identity mismatch: "
                        f"g={g}, sigma={sigma}, N={packet_count}, h={support}, "
                        f"z={z}, direct={direct}, formula={formula}"
                    )
        rows.append(
            {
                "g": g,
                "sigma": sigma,
                "packet_positions": packet_count,
                "weight_tilts": [str(value) for value in tilts],
                "supports_checked": [1, packet_count],
                "status": "EXACT_MATCH",
            }
        )
    return rows


def log_a_function(radius: float, u: float, a: float) -> float:
    return (
        math.log1p(-a * u * radius)
        - math.log1p(-a * radius)
        - math.log1p(-u * radius)
    )


def log_b_function(radius: float, u: float, q: float, a: float) -> float:
    log_den_a = math.log1p(-a * radius)
    log_den_u = math.log1p(-u * radius)
    terms = [
        0.0,
        math.log(radius) - log_den_a,
        math.log(q) + math.log(u) + 2.0 * math.log(radius) - log_den_a - log_den_u,
    ]
    return float(logsumexp(terms))


def log_one_minus_exp(log_value: float) -> float:
    """Return log(1-exp(log_value)) for log_value <= 0."""
    if log_value > 0.0:
        if log_value < 1e-12:
            log_value = 0.0
        else:
            raise ArithmeticError(f"positive log ratio in subtraction: {log_value}")
    if log_value == 0.0:
        return -math.inf
    if log_value < -math.log(2.0):
        return math.log1p(-math.exp(log_value))
    return math.log(-math.expm1(log_value))


def log_h_function(radius: float, u: float, q: float, a: float) -> float:
    log_a_r = log_a_function(radius, u, a)
    log_a_ar = log_a_function(a * radius, u, a)
    log_b_r = log_b_function(radius, u, q, a)
    log_c_r = -math.log1p(-u * radius)

    first = math.log(a) + log_a_ar + log_b_r
    ratio = math.log(a) + log_a_ar - log_a_r
    second = log_a_r + log_one_minus_exp(ratio) + log_c_r
    return float(np.logaddexp(first, second))


def tilted_log_bound(
    n: int,
    g: int,
    sigma: int,
    support: int,
    distance: int,
    weight_tilt: float,
    radius_fraction: float,
) -> float:
    q = math.ldexp(1.0, -sigma)
    a = 1.0 - q
    b = math.ldexp((1.0 + weight_tilt) ** g, -g)
    u = 1.0 / b
    radius = b * radius_fraction
    log_a_r = log_a_function(radius, u, a)
    log_h_r = log_h_function(radius, u, q, a)
    return (
        -distance * math.log(weight_tilt)
        + n * math.log(b)
        - log_binomial(n, support)
        + log_h_r
        + (support - 1) * log_a_r
        - (n - support) * math.log(radius)
    )


def optimize_bound(n: int, g: int, sigma: int, support: int, distance: int) -> dict[str, object]:
    if not 1 <= support <= n:
        raise ValueError("support must lie in [1,n]")

    def optimize_radius(log_surprisal: float) -> tuple[float, float, float]:
        weight_tilt = math.exp(-math.exp(log_surprisal))

        def objective(radius_fraction: float) -> float:
            return tilted_log_bound(
                n, g, sigma, support, distance, weight_tilt, radius_fraction
            )

        result = minimize_scalar(
            objective,
            bounds=(1e-10, 1.0 - 1e-10),
            method="bounded",
            options={"xatol": 2e-12, "maxiter": 500},
        )
        return float(result.fun), weight_tilt, float(result.x)

    grid = np.linspace(-14.0, 2.0, 65)
    grid_rows = [optimize_radius(float(point)) for point in grid]
    index = min(range(len(grid_rows)), key=lambda i: grid_rows[i][0])
    low = float(grid[max(0, index - 1)])
    high = float(grid[min(len(grid) - 1, index + 1)])

    result = minimize_scalar(
        lambda value: optimize_radius(float(value))[0],
        bounds=(low, high),
        method="bounded",
        options={"xatol": 2e-9, "maxiter": 200},
    )
    log_bound, weight_tilt, radius_fraction = optimize_radius(float(result.x))
    b = math.ldexp((1.0 + weight_tilt) ** g, -g)
    return {
        "diagnostic_log2_upper": min(0.0, log_bound / math.log(2.0)),
        "weight_tilt": weight_tilt,
        "coefficient_radius": b * radius_fraction,
        "radius_fraction_of_singularity": radius_fraction,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
    }


def ideal_live_chernoff(n: int, g: int, distance: int) -> float:
    total_bits = n * g
    if distance >= total_bits / 2:
        return 0.0
    z = distance / (total_bits - distance)
    b = math.ldexp((1.0 + z) ** g, -g)
    return (-distance * math.log(z) + n * math.log(b)) / math.log(2.0)


def cyclic_cluster_witness(
    n: int,
    g: int,
    sigma: int,
    support: int,
    distance: int,
    reset_prefix_limit: int = 180_000,
) -> dict[str, object]:
    """Give a valid floating lower bound from one clustered-support event.

    A fixed cyclic interval of R positions crosses the retained-lap boundary.
    Both endpoints adjacent to the complementary gap are required active. A
    reset after k live zero packets in that gap leaves the rest of the gap off.
    Omitting reset times beyond reset_prefix_limit preserves a lower bound.
    """
    if support < 2:
        raise ValueError("the endpoint witness requires support at least two")
    q = math.ldexp(1.0, -sigma)
    a = 1.0 - q
    log_a = math.log(a)

    def score(radius: int) -> tuple[float, float, int]:
        gap = n - radius
        retained_reset_times = min(gap, reset_prefix_limit)
        resets = np.arange(retained_reset_times, dtype=np.float64)
        log_reset_weight = resets * log_a + math.log(q)
        log_low_weight = binom.logcdf(distance, g * (radius + resets), 0.5)
        log_inner_event = float(logsumexp(log_reset_weight + log_low_weight))
        log_support_event = (
            log_binomial(radius - 2, support - 2) - log_binomial(n, support)
        )
        return log_support_event + log_inner_event, log_inner_event, retained_reset_times

    low = max(support, 90_000)
    high = min(n - 1, 105_000)
    best = (-math.inf, low, -math.inf, 0)
    for radius in range(low, high + 1, 250):
        value, inner, count = score(radius)
        if value > best[0]:
            best = (value, radius, inner, count)
    center = best[1]
    for radius in range(max(low, center - 300), min(high, center + 300) + 1, 10):
        value, inner, count = score(radius)
        if value > best[0]:
            best = (value, radius, inner, count)
    return {
        "input_packet_support": support,
        "cyclic_interval_packets": best[1],
        "cyclic_interval_fraction": best[1] / n,
        "reset_times_included": best[3],
        "diagnostic_log2_probability_lower": best[0] / math.log(2.0),
        "diagnostic_log2_reset_and_weight_factor": best[2] / math.log(2.0),
        "status": "FLOATING_LOWER_BOUND_WITH_TRUNCATED_POSITIVE_SUM",
    }


def target_rows(supports: Iterable[int], deltas: Iterable[float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for delta in deltas:
        distance = math.floor(delta * TARGET_BITS)
        ideal = ideal_live_chernoff(TARGET_PACKETS, TARGET_G, distance)
        for support in supports:
            optimized = optimize_bound(
                TARGET_PACKETS, TARGET_G, TARGET_SIGMA, support, distance
            )
            row: dict[str, object] = {
                "relative_binary_weight": delta,
                "distance": distance,
                "input_packet_support": support,
                "ideal_all_live_chernoff_log2": ideal,
            }
            row.update(optimized)
            rows.append(row)
            print(
                f"delta,{delta:.6f},h,{support},D,{distance},"
                f"log2_upper,{optimized['diagnostic_log2_upper']:.9f},"
                f"ideal_live,{ideal:.9f},"
                f"z,{optimized['weight_tilt']:.17g},"
                f"r,{optimized['coefficient_radius']:.17g},"
                f"rho,{optimized['radius_fraction_of_singularity']:.17g}",
                flush=True,
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("check", "target", "all"), default="all")
    parser.add_argument("--supports", nargs="+", type=int, default=list(DEFAULT_SUPPORTS))
    parser.add_argument("--relative-distances", nargs="+", type=float, default=list(DEFAULT_DELTAS))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checks = exact_checks() if args.mode in ("check", "all") else []
    for row in checks:
        print(
            f"exact_check,g,{row['g']},sigma,{row['sigma']},"
            f"N,{row['packet_positions']},status,{row['status']}",
            flush=True,
        )
    rows = (
        target_rows(args.supports, args.relative_distances)
        if args.mode in ("target", "all")
        else []
    )
    cluster_witnesses = []
    if args.mode in ("target", "all") and 16_384 in args.supports:
        for delta in args.relative_distances:
            distance = math.floor(delta * TARGET_BITS)
            witness = cyclic_cluster_witness(
                TARGET_PACKETS, TARGET_G, TARGET_SIGMA, 16_384, distance
            )
            witness["relative_binary_weight"] = delta
            witness["distance"] = distance
            cluster_witnesses.append(witness)
            print(
                f"cluster_witness,delta,{delta:.6f},h,16384,"
                f"R,{witness['cyclic_interval_packets']},"
                f"log2_lower,{witness['diagnostic_log2_probability_lower']:.9f}",
                flush=True,
            )
    payload = {
        "schema": "riffle-randomstepconv-2lap-goal01-v1",
        "construction": "Riffle RandomStepConv-2Lap g=4 sigma=20",
        "evidence": {
            "gap_identity": "EXACT",
            "small_instance_checks": "EXHAUSTIVE_OVER_SUPPORTS_AND_EXACT_IN_TWO_STATE_QUOTIENT",
            "target_optimizer": "FLOATING_DIAGNOSTIC",
            "target_bound": "NOT_OUTWARD_ROUNDED",
        },
        "parameters": {
            "g": TARGET_G,
            "sigma": TARGET_SIGMA,
            "packet_positions": TARGET_PACKETS,
            "binary_output_length": TARGET_BITS,
            "laps_computed": 2,
            "laps_retained": 1,
        },
        "exact_checks": checks,
        "target_rows": rows,
        "cyclic_cluster_witnesses": cluster_witnesses,
        "scope": (
            "The receipt validates the exact retained-lap gap identity on small "
            "instances and records ordinary floating-point inner bounds at target size. "
            "It does not provide an outward-rounded or end-to-end distance certificate."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
