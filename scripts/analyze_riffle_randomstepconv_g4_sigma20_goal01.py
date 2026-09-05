#!/usr/bin/env python3
"""Exact checks and tilted bounds for RandomStepConv Goal 01.

The small-instance checker enumerates every setup matrix.  The target-size
path evaluates the two-state transfer matrix without expanding coefficients.
Floating-point optimization chooses tilts.  Decimal arithmetic then evaluates
an upper bound at the reported decimal tilts with directed rounding.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp


TARGET_G = 4
TARGET_SIGMA = 20
TARGET_PACKETS = 524_352
TARGET_BITS = TARGET_G * TARGET_PACKETS
DEFAULT_SUPPORTS = (18, 24, 32, 48, 64, 92, 128, 256, 512, 1024, 2048)
DEFAULT_DELTAS = (0.05, 0.09, 0.12)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_randomstepconv_g4_sigma20/receipts/"
    "goal01_inner_transfer_bound.json"
)


def log_binomial(n: int, h: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(h + 1) - math.lgamma(n - h + 1)


def episode_terms(
    n: int, g: int, sigma: int, h: int, weight_tilt: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return log bounds and their coefficient-saddle radii by termination count."""
    q = math.ldexp(1.0, -sigma)
    b = math.ldexp((1.0 + weight_tilt) ** g, -g)
    d = (1.0 - q) * b
    alpha = q * b / (1.0 - d)
    beta = d * (1.0 - b) / (1.0 - d)
    alpha = max(np.finfo(float).tiny, min(1.0, alpha))
    beta = max(np.finfo(float).tiny, beta)
    e = np.arange(h + 1, dtype=np.float64)
    k = h - e
    remaining = n - h
    p = remaining * (1.0 + d) + (e + 1.0) + k * d
    discriminant = np.maximum(0.0, p * p - 4.0 * d * (n + 1.0) * remaining)
    radii = 2.0 * remaining / (p + np.sqrt(discriminant))
    radii = np.clip(radii, 1e-300, 1.0 - 1e-15)
    log_choose_e = gammaln(h + 1.0) - gammaln(e + 1.0) - gammaln(k + 1.0)
    log_mixture = log_choose_e + e * math.log(alpha) + k * math.log(beta)
    log_coefficient = (
        -remaining * np.log(radii)
        - (e + 1.0) * np.log1p(-radii)
        - k * np.log1p(-d * radii)
    )
    return log_mixture + log_coefficient, radii


def optimize_tilts(n: int, g: int, sigma: int, h: int, distance: int) -> dict[str, float | bool | str | int]:
    if not 0 < h <= n:
        raise ValueError("support must lie in [1,n]")
    if not 0 <= distance < g * n:
        raise ValueError("distance must lie in [0,g*n)")
    log_choose = log_binomial(n, h)

    def objective(log_surprisal: float) -> float:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        terms, _radii = episode_terms(n, g, sigma, h, z)
        return float(logsumexp(terms)) - log_choose + distance * surprisal

    grid = np.linspace(-18.0, 3.0, 85)
    values = np.asarray([objective(float(point)) for point in grid])
    index = int(np.argmin(values))
    low = float(grid[max(0, index - 1)])
    high = float(grid[min(grid.size - 1, index + 1)])
    result = minimize_scalar(
        objective,
        bounds=(low, high),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 500},
    )
    best_log_surprisal = float(result.x)
    best_z = math.exp(-math.exp(best_log_surprisal))
    terms, radii = episode_terms(n, g, sigma, h, best_z)
    dominant = int(np.argmax(terms))
    diagnostic_log2 = min(0.0, float(result.fun) / math.log(2.0))
    return {
        "weight_tilt": best_z,
        "dominant_termination_count": dominant,
        "dominant_gap_tilt": float(radii[dominant]),
        "diagnostic_log2_upper": diagnostic_log2,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
    }


def decimal_context(rounding: str) -> Context:
    context = Context(prec=90, rounding=rounding, Emin=-9_999_999, Emax=9_999_999)
    for signal in context.traps:
        context.traps[signal] = False
    return context


def outward_bound(
    n: int,
    g: int,
    sigma: int,
    h: int,
    distance: int,
    weight_tilt: float,
) -> dict[str, str]:
    # The decimal string defines an exact rational output tilt.  Each
    # termination component uses its own exact decimal coefficient radius.
    z = Decimal(format(weight_tilt, ".17g"))
    one = Decimal(1)
    q = one / Decimal(1 << sigma)
    floor_context = decimal_context(ROUND_FLOOR)
    ceil_context = decimal_context(ROUND_CEILING)
    with localcontext(floor_context):
        b_lower = (one + z) ** g / Decimal(1 << g)
    with localcontext(ceil_context):
        b_upper = (one + z) ** g / Decimal(1 << g)
        d_upper = (one - q) * b_upper
        alpha_upper = q * b_upper / (one - d_upper)
        beta_upper = d_upper * (one - b_lower) / (one - d_upper)
    _terms, radii = episode_terms(n, g, sigma, h, float(z))
    sum_upper = Decimal(0)
    remaining = n - h
    for e in range(h + 1):
        k = h - e
        r = Decimal(format(float(radii[e]), ".17g"))
        with localcontext(ceil_context):
            term_numerator = (
                Decimal(math.comb(h, e))
                * (alpha_upper**e)
                * (beta_upper**k)
            )
        with localcontext(floor_context):
            term_denominator = (
                (r**remaining)
                * ((one - r) ** (e + 1))
                * ((one - d_upper * r) ** k)
            )
        with localcontext(ceil_context):
            sum_upper += term_numerator / term_denominator
    with localcontext(floor_context):
        denominator = Decimal(math.comb(n, h)) * (z**distance)
    with localcontext(ceil_context):
        ratio = sum_upper / denominator
        if ratio >= 1:
            log2_upper = Decimal(0)
            ratio = min(ratio, Decimal(1))
        else:
            numerator_log = ratio.ln(context=ceil_context)
            ln2_upper = Decimal(2).ln(context=ceil_context)
            log2_upper = numerator_log / ln2_upper
    return {
        "verified_weight_tilt_decimal": str(z),
        "outward_probability_upper": str(ratio),
        "outward_log2_upper": str(log2_upper),
    }


def transfer_coefficients(g: int, sigma: int, n: int) -> dict[tuple[int, int], Fraction]:
    """Exact coefficient sums over support patterns, with terminal state discarded."""
    q = Fraction(1, 1 << sigma)
    packet = tuple(Fraction(math.comb(g, weight), 1 << g) for weight in range(g + 1))
    current: dict[tuple[int, int, bool], Fraction] = {(0, 0, False): Fraction(1)}
    for _ in range(n):
        following: dict[tuple[int, int, bool], Fraction] = defaultdict(Fraction)
        for (support, output_weight, live), mass in current.items():
            for active in (False, True):
                if not active and not live:
                    following[(support, output_weight, False)] += mass
                    continue
                for weight, output_probability in enumerate(packet):
                    common = mass * output_probability
                    following[(support + int(active), output_weight + weight, False)] += common * q
                    following[(support + int(active), output_weight + weight, True)] += common * (1 - q)
        current = following
    result: dict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (support, output_weight, _live), mass in current.items():
        result[(support, output_weight)] += mass
    return dict(result)


def apply_binary_matrix(matrix_bits: int, dimension: int, vector: int) -> int:
    result = 0
    row_mask = (1 << dimension) - 1
    for row in range(dimension):
        row_bits = (matrix_bits >> (row * dimension)) & row_mask
        result |= ((row_bits & vector).bit_count() & 1) << row
    return result


def brute_map_coefficients(
    g: int, sigma: int, n: int, active_value: int
) -> dict[tuple[int, int], Fraction]:
    """Enumerate every matrix sequence and input support for tiny parameters."""
    if not 0 < active_value < (1 << g):
        raise ValueError("active packet value must be nonzero and fit in g bits")
    dimension = g + sigma
    map_count = 1 << (dimension * dimension)
    vector_count = 1 << dimension
    transitions = [
        tuple(apply_binary_matrix(matrix, dimension, vector) for vector in range(vector_count))
        for matrix in range(map_count)
    ]
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for support_mask in range(1 << n):
        support = support_mask.bit_count()
        for matrices in product(range(map_count), repeat=n):
            state = 0
            output_weight = 0
            for position, matrix in enumerate(matrices):
                packet = active_value if (support_mask >> position) & 1 else 0
                vector = packet | (state << g)
                image = transitions[matrix][vector]
                output_weight += (image & ((1 << g) - 1)).bit_count()
                state = image >> g
            counts[(support, output_weight)] += 1
    denominator = map_count**n
    return {key: Fraction(count, denominator) for key, count in counts.items()}


def exact_gap_moment(g: int, sigma: int, n: int, support: int, z: Fraction) -> Fraction:
    """Evaluate the rank-one gap formula by exhaustive gap composition."""
    if support == 0:
        return Fraction(1)
    q = Fraction(1, 1 << sigma)
    b = (1 + z) ** g / (1 << g)
    d = (1 - q) * b
    alpha = q * b / (1 - d)
    values = tuple(alpha + (1 - alpha) * d ** (gap + 1) for gap in range(n - support + 1))
    total = Fraction(0)
    for gaps in product(range(n - support + 1), repeat=support):
        if sum(gaps) > n - support:
            continue
        term = Fraction(1)
        for gap in gaps:
            term *= values[gap]
        total += term
    return total


def exact_checks() -> list[dict[str, object]]:
    cases = (
        (1, 1, 3, 1),
        (2, 1, 2, 1),
        (2, 1, 2, 3),
    )
    rows = []
    for g, sigma, n, active_value in cases:
        transfer = transfer_coefficients(g, sigma, n)
        brute = brute_map_coefficients(g, sigma, n, active_value)
        keys = set(transfer) | set(brute)
        mismatches = [
            (key, transfer.get(key, Fraction(0)), brute.get(key, Fraction(0)))
            for key in sorted(keys)
            if transfer.get(key, Fraction(0)) != brute.get(key, Fraction(0))
        ]
        if mismatches:
            raise AssertionError(f"exact map enumeration mismatch: {mismatches[:3]}")
        z = Fraction(1, 2)
        gap_mismatches = []
        for support in range(n + 1):
            transfer_moment = sum(
                mass * z**weight
                for (row_support, weight), mass in transfer.items()
                if row_support == support
            )
            gap_moment = exact_gap_moment(g, sigma, n, support, z)
            if transfer_moment != gap_moment:
                gap_mismatches.append((support, transfer_moment, gap_moment))
        if gap_mismatches:
            raise AssertionError(f"rank-one gap identity mismatch: {gap_mismatches[:3]}")
        rows.append(
            {
                "g": g,
                "sigma": sigma,
                "packet_positions": n,
                "active_packet_value": active_value,
                "maps_per_position": 1 << ((g + sigma) ** 2),
                "matrix_sequences": (1 << ((g + sigma) ** 2)) ** n,
                "coefficient_count": len(keys),
                "status": "EXACT_MATCH",
                "rank_one_gap_identity_at_z_half": "EXACT_MATCH",
            }
        )
    return rows


def target_rows(supports: Iterable[int], deltas: Iterable[float]) -> list[dict[str, object]]:
    rows = []
    for delta in deltas:
        distance = math.floor(delta * TARGET_BITS)
        for support in supports:
            optimized = optimize_tilts(
                TARGET_PACKETS,
                TARGET_G,
                TARGET_SIGMA,
                support,
                distance,
            )
            verified = outward_bound(
                TARGET_PACKETS,
                TARGET_G,
                TARGET_SIGMA,
                support,
                distance,
                float(optimized["weight_tilt"]),
            )
            row: dict[str, object] = {
                "relative_binary_weight": delta,
                "distance": distance,
                "input_packet_support": support,
            }
            row.update(optimized)
            row.update(verified)
            rows.append(row)
            print(
                f"delta,{delta:.6f},h,{support},D,{distance},"
                f"diagnostic_log2,{optimized['diagnostic_log2_upper']:.9f},"
                f"outward_log2,{verified['outward_log2_upper']},"
                f"e,{optimized['dominant_termination_count']},"
                f"r,{optimized['dominant_gap_tilt']:.17g},"
                f"z,{verified['verified_weight_tilt_decimal']}",
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
    if checks:
        for row in checks:
            print(
                f"exact_check,g,{row['g']},sigma,{row['sigma']},N,{row['packet_positions']},"
                f"active,{row['active_packet_value']},status,{row['status']}",
                flush=True,
            )
    rows = (
        target_rows(args.supports, args.relative_distances)
        if args.mode in ("target", "all")
        else []
    )
    payload = {
        "schema": "riffle-randomstepconv-goal01-v1",
        "construction": "Riffle RandomStepConv g=4 sigma=20",
        "evidence": {
            "transfer_identity": "EXACT",
            "rank_one_gap_identity": "EXACT",
            "small_map_checks": "EXHAUSTIVE",
            "target_optimizer": "FLOATING_DIAGNOSTIC",
            "target_evaluation_at_reported_tilts": "DECIMAL_OUTWARD_ROUNDED",
        },
        "parameters": {
            "g": TARGET_G,
            "sigma": TARGET_SIGMA,
            "packet_positions": TARGET_PACKETS,
            "binary_output_length": TARGET_BITS,
            "terminal_state": "discarded",
        },
        "exact_checks": checks,
        "target_rows": rows,
        "scope": (
            "The receipt bounds only the RandomStepConv inner for fixed input-packet "
            "support. It does not consume the outer spectrum or certify an end-to-end distance."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
