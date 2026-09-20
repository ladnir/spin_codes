#!/usr/bin/env python3
"""Interval certificate for the bulk of a full-length designed BA-3 code.

The construction is

    direct_sum(C), uniform interleaver, accumulator,
    uniform interleaver, accumulator.

For normalized weights x, y, and z, every individual expected-spectrum term
is bounded by

    (N+1)^4 2^(N F(x,y,z)),

where

    F(x,y,z) = a_C(x) + p_acc(x,y) + p_acc(y,z).

Here a_C is a coefficient Chernoff bound for W_C^m and p_acc is the entropy
bound for an accumulator transition.  There are fewer than (N+1)^3 triples,
so the complete bulk is at most

    (N+1)^7 2^(N max F).

For z < 1/2, p_acc(y,z) increases with z.  The maximum through distance D is
therefore attained at z=D/N.  The remaining triangular (x,y) domain is
parameterized by

    y = x/2 + t (2 delta - x/2),  0 <= t <= 1.

This script covers that rectangle with boxes and evaluates every box with
mpmath interval arithmetic.  The coefficient bound uses one fixed positive r
per x box.  It remains valid even if the floating-point saddle used to choose
r is inaccurate, because every positive r gives a Chernoff upper bound.
"""

from __future__ import annotations

import argparse
from collections import deque
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp
from scipy.optimize import brentq


@dataclass(frozen=True)
class Spectrum:
    length: int
    dimension: int
    counts: tuple[int, ...]


def read_spectrum(path: Path, local_length: int, local_dimension: int) -> Spectrum:
    counts = [0] * (local_length + 1)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            counts[int(row["weight"])] += int(row["count"])
    if counts[0] != 1 or sum(counts) != 1 << local_dimension:
        raise ValueError("invalid constituent spectrum")
    return Spectrum(local_length, local_dimension, tuple(counts))


def mean_weight(spectrum: Spectrum, log_r: float) -> float:
    terms = [
        count * math.exp(weight * log_r)
        for weight, count in enumerate(spectrum.counts)
        if count
    ]
    weighted = [
        weight * count * math.exp(weight * log_r)
        for weight, count in enumerate(spectrum.counts)
        if count
    ]
    return math.fsum(weighted) / math.fsum(terms)


def choose_r(spectrum: Spectrum, normalized_weight: float) -> float:
    target = spectrum.length * normalized_weight

    def equation(log_r: float) -> float:
        return mean_weight(spectrum, log_r) - target

    return math.exp(brentq(equation, -60.0, 0.0))


def endpoint(value: object, upper: bool) -> mp.mpf:
    interval = value
    return mp.mpf(interval.b if upper else interval.a)


def point_phi(value: mp.mpf):
    if value == 0:
        return mp.iv.mpf([0, 0])
    q = mp.iv.mpf([value, value])
    return q * mp.iv.log(q) / mp.iv.log(2)


def phi_range(value):
    """Range enclosure of q log2(q) for an interval in [0,1]."""

    lo = endpoint(value, False)
    hi = endpoint(value, True)
    tolerance = mp.mpf("1e-45")
    if lo < 0 and lo > -tolerance:
        lo = mp.mpf(0)
    if hi > 1 and hi < 1 + tolerance:
        hi = mp.mpf(1)
    if lo < 0 or hi > 1 or lo > hi:
        raise ValueError(f"phi domain escaped [0,1]: [{lo}, {hi}]")

    lo_value = point_phi(lo)
    hi_value = point_phi(hi)
    lower = min(endpoint(lo_value, False), endpoint(hi_value, False))
    upper = max(endpoint(lo_value, True), endpoint(hi_value, True))
    stationary = 1 / mp.e
    if lo <= stationary <= hi:
        stationary_value = point_phi(stationary)
        lower = min(lower, endpoint(stationary_value, False))
        upper = max(upper, endpoint(stationary_value, True))
    return mp.iv.mpf([lower, upper])


def accumulator_exponent_box(input_weight, output_weight, one_minus_output,
                             first_gap, second_gap):
    """Interval extension of p_acc in expanded q log(q) form."""

    return (
        input_weight
        + phi_range(one_minus_output)
        + phi_range(output_weight)
        - phi_range(first_gap)
        - phi_range(second_gap)
        + phi_range(1 - input_weight)
    )


def outer_exponent_box(spectrum: Spectrum, x, r_float: float):
    """Chernoff upper bound valid for every x in the supplied interval."""

    r_text = format(r_float, ".17g")
    r = mp.iv.mpf([r_text, r_text])
    polynomial = mp.iv.mpf([0, 0])
    for weight, count in enumerate(spectrum.counts):
        if count:
            polynomial += count * r**weight
    return (
        mp.iv.log(polynomial) / (spectrum.length * mp.iv.log(2))
        - x * mp.iv.log(r) / mp.iv.log(2)
    )


def box_upper(spectrum: Spectrum, delta, x_lo: float, x_hi: float,
              t_lo: float, t_hi: float) -> float:
    x = mp.iv.mpf([format(x_lo, ".17g"), format(x_hi, ".17g")])
    t = mp.iv.mpf([format(t_lo, ".17g"), format(t_hi, ".17g")])

    two_delta_minus_half_x = 2 * delta - x / 2
    y = x / 2 + t * two_delta_minus_half_x
    one_minus_y = 1 - x / 2 - t * two_delta_minus_half_x
    first_gap = 1 - x - t * two_delta_minus_half_x
    second_gap = t * two_delta_minus_half_x

    p_first = accumulator_exponent_box(
        x, y, one_minus_y, first_gap, second_gap
    )

    half_y = x / 4 + t * (delta - x / 4)
    final_first_gap = 1 - delta - half_y
    final_second_gap = (1 - t) * (delta - x / 4)
    p_second = accumulator_exponent_box(
        y,
        delta,
        1 - delta,
        final_first_gap,
        final_second_gap,
    )

    r = choose_r(spectrum, (x_lo + x_hi) / 2)
    outer = outer_exponent_box(spectrum, x, r)
    total = outer + p_first + p_second
    # Preserve the interval upper endpoint when converting back to binary64.
    return math.nextafter(float(total.b), math.inf)


def certify(spectrum: Spectrum, length: int, distance: int, sparse_cap: int,
            x_boxes: int, t_boxes: int, target_exponent: float,
            max_depth: int, max_box_evaluations: int) -> dict[str, object]:
    mp.mp.dps = 60
    mp.iv.dps = 60
    delta_text = f"{distance}/{length}"
    delta = mp.iv.mpf(distance) / length
    x_min = (sparse_cap + 1) / length
    x_max = 4 * distance / length

    queue: deque[tuple[float, float, float, float, int]] = deque()
    for x_index in range(x_boxes):
        x_lo = x_min + (x_max - x_min) * x_index / x_boxes
        x_hi = x_min + (x_max - x_min) * (x_index + 1) / x_boxes
        for t_index in range(t_boxes):
            t_lo = t_index / t_boxes
            t_hi = (t_index + 1) / t_boxes
            queue.append((x_lo, x_hi, t_lo, t_hi, 0))

    evaluations = 0
    accepted = 0
    subdivided = 0
    deepest = 0
    maximum_accepted = -math.inf
    worst_accepted: dict[str, float] | None = None
    unresolved: dict[str, float] | None = None
    while queue:
        if evaluations >= max_box_evaluations:
            unresolved = {"reason": "box evaluation limit reached"}
            break
        x_lo, x_hi, t_lo, t_hi, depth = queue.popleft()
        upper = box_upper(spectrum, delta, x_lo, x_hi, t_lo, t_hi)
        evaluations += 1
        deepest = max(deepest, depth)
        if upper <= target_exponent:
            accepted += 1
            if upper > maximum_accepted:
                maximum_accepted = upper
                worst_accepted = {
                    "x_lo": x_lo,
                    "x_hi": x_hi,
                    "t_lo": t_lo,
                    "t_hi": t_hi,
                    "depth": depth,
                    "exponent_upper": upper,
                }
            continue
        if depth >= max_depth:
            unresolved = {
                "reason": "maximum subdivision depth reached",
                "x_lo": x_lo,
                "x_hi": x_hi,
                "t_lo": t_lo,
                "t_hi": t_hi,
                "depth": depth,
                "exponent_upper": upper,
            }
            break
        x_mid = (x_lo + x_hi) / 2
        t_mid = (t_lo + t_hi) / 2
        next_depth = depth + 1
        queue.extend(
            (
                (x_lo, x_mid, t_lo, t_mid, next_depth),
                (x_lo, x_mid, t_mid, t_hi, next_depth),
                (x_mid, x_hi, t_lo, t_mid, next_depth),
                (x_mid, x_hi, t_mid, t_hi, next_depth),
            )
        )
        subdivided += 1

    complete = unresolved is None and not queue
    maximum = target_exponent if complete else math.inf
    type_overhead = 7 * math.log2(length + 1)
    bulk_log2_upper = length * maximum + type_overhead
    return {
        "schema": "full-designed-ba3-bulk-interval-certificate-v1",
        "construction": {
            "length": length,
            "dimension": length // 2,
            "local_length": spectrum.length,
            "local_dimension": spectrum.dimension,
        },
        "target": {
            "distance": distance,
            "relative_distance": distance / length,
            "relative_distance_exact": delta_text,
            "sparse_outer_weight_cap": sparse_cap,
        },
        "domain": {
            "outer_relative_weight_min": x_min,
            "outer_relative_weight_max": x_max,
            "parameterization": "y=x/2+t*(2*delta-x/2), t in [0,1]",
            "final_weight": "z=delta by monotonicity for z<1/2",
        },
        "interval_cover": {
            "initial_x_boxes": x_boxes,
            "initial_t_boxes": t_boxes,
            "initial_box_count": x_boxes * t_boxes,
            "adaptive_acceptance_threshold": target_exponent,
            "box_evaluations": evaluations,
            "accepted_leaf_boxes": accepted,
            "subdivided_boxes": subdivided,
            "maximum_depth_used": deepest,
            "maximum_depth_allowed": max_depth,
            "mpmath_interval_decimal_precision": mp.iv.dps,
            "worst_accepted_box": worst_accepted,
            "unresolved_box": unresolved,
        },
        "finite_union_bound": {
            "transition_prefactor": "(N+1)^4",
            "weight_triple_count": "less than (N+1)^3",
            "total_prefactor": "(N+1)^7",
            "log2_prefactor_bits": type_overhead,
            "certified_max_exponent_per_output_bit_upper": maximum,
            "largest_accepted_box_upper": maximum_accepted,
            "scaled_exponent_upper_bits": length * maximum,
            "bulk_log2_expected_bad_mass_upper": bulk_log2_upper,
            "bulk_margin_bits_lower": -bulk_log2_upper,
        },
        "status": (
            "certificate passes"
            if complete and bulk_log2_upper < 0
            else "certificate incomplete"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--local-length", type=int, required=True)
    parser.add_argument("--local-dimension", type=int, required=True)
    parser.add_argument("--length", type=int, default=1 << 21)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument("--sparse-cap", type=int, default=256)
    parser.add_argument("--x-boxes", type=int, default=4)
    parser.add_argument("--t-boxes", type=int, default=4)
    parser.add_argument("--target-exponent", type=float, default=-1e-4)
    parser.add_argument("--max-depth", type=int, default=16)
    parser.add_argument("--max-box-evaluations", type=int, default=500000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    spectrum = read_spectrum(
        args.spectrum, args.local_length, args.local_dimension
    )
    result = certify(
        spectrum,
        args.length,
        args.distance,
        args.sparse_cap,
        args.x_boxes,
        args.t_boxes,
        args.target_exponent,
        args.max_depth,
        args.max_box_evaluations,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
