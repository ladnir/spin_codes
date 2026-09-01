#!/usr/bin/env python3
"""Outward certificate for a concave Golay--BA-3 spectrum majorant.

The limiting expected spectrum exponent is

    sup_{a,b} g(a) + p(a,b) + p(b,w).

The verifier covers each weight interval by one affine support line.  All
box endpoints and support-line coefficients are rational.  Transcendental
operations use mpmath interval arithmetic, and conversion to binary64 is
expanded outwards by one representable number.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from mpmath import iv


iv.dps = 70
ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)
LEFT_ENDPOINT = Fraction(13, 125)
RIGHT_ENDPOINT = ONE - LEFT_ENDPOINT
LN2 = iv.log(iv.mpf(2))

# Each pair is (slope, intercept) on the left half.  The zero-slope line is
# exact: the spectrum contains at most 2^(B/2) words.  The other intercepts
# include diagnostic slack; this program proves their validity.
LEFT_LINES = [
    (Fraction(1666, 1000), Fraction("-0.173244")),
    (Fraction(166, 100), Fraction("-0.172329198693950")),
    (Fraction(165, 100), Fraction("-0.170746000276704")),
    (Fraction(163, 100), Fraction("-0.167527934405350")),
    (Fraction(16, 10), Fraction("-0.162588357378692")),
    (Fraction(155, 100), Fraction("-0.154045203238254")),
    (Fraction(14, 10), Fraction("-0.126133170725692")),
    (Fraction(12, 10), Fraction("-0.0832718284462267")),
    (Fraction(1), Fraction("-0.033294157620498")),
    (Fraction(8, 10), Fraction("0.0245470741467448")),
    (Fraction(6, 10), Fraction("0.0909341020092481")),
    (Fraction(4, 10), Fraction("0.166460678214840")),
    (Fraction(2, 10), Fraction("0.251585223433616")),
    (Fraction(1, 10), Fraction("0.297841147503783")),
    (Fraction(5, 100), Fraction("0.321905562625601")),
]


def as_interval(value: Fraction):
    return iv.mpf(value.numerator) / value.denominator


def interval(lo: Fraction, hi: Fraction):
    return iv.mpf([as_interval(lo).a, as_interval(hi).b])


def lower_float(value) -> float:
    return math.nextafter(float(value.a), -math.inf)


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def add_upper(*values: float) -> float:
    return math.nextafter(math.fsum(values), math.inf)


def xlogx(value):
    return value * iv.log(value)


@lru_cache(maxsize=None)
def xlogx_point(value: Fraction) -> tuple[float, float]:
    if value == 0:
        return 0.0, 0.0
    result = xlogx(as_interval(value))
    return lower_float(result), upper_float(result)


def xlogx_range(lo: Fraction, hi: Fraction) -> tuple[float, float]:
    if lo < 0 or hi > 1 or lo > hi:
        raise ValueError("invalid x log x interval")
    lo_bounds = xlogx_point(lo)
    hi_bounds = xlogx_point(hi)
    upper = math.nextafter(max(lo_bounds[1], hi_bounds[1]), math.inf)
    inv_e = iv.mpf(1) / iv.e
    inv_e_lo = Fraction(str(lower_float(inv_e)))
    inv_e_hi = Fraction(str(upper_float(inv_e)))
    if lo <= inv_e_lo and inv_e_hi <= hi:
        lower = lower_float(-iv.mpf(1) / iv.e)
    else:
        lower = math.nextafter(min(lo_bounds[0], hi_bounds[0]), -math.inf)
    return lower, upper


def p_point(a: Fraction, b: Fraction):
    aa = as_interval(a)
    bb = as_interval(b)
    return (
        aa * LN2
        + xlogx(bb)
        + xlogx(1 - bb)
        - xlogx(bb - aa / 2)
        - xlogx(1 - bb - aa / 2)
        + xlogx(1 - aa)
    )


def p_gradient(a, b):
    left = b - a / 2
    right = 1 - b - a / 2
    return (
        iv.log(2 * iv.sqrt(left * right) / (1 - a)),
        iv.log(b * right / ((1 - b) * left)),
    )


@lru_cache(maxsize=None)
def p_taylor_upper(
    a_lo: Fraction,
    a_hi: Fraction,
    b_lo: Fraction,
    b_hi: Fraction,
) -> float:
    if (
        b_lo - a_hi * HALF <= 0
        or ONE - b_hi - a_hi * HALF <= 0
        or a_hi >= 1
        or b_lo <= 0
        or b_hi >= 1
    ):
        return math.inf
    a_mid = (a_lo + a_hi) * HALF
    b_mid = (b_lo + b_hi) * HALF
    if not a_mid * HALF < b_mid < ONE - a_mid * HALF:
        return math.inf
    a = interval(a_lo, a_hi)
    b = interval(b_lo, b_hi)
    derivative_a, derivative_b = p_gradient(a, b)
    enclosure = (
        p_point(a_mid, b_mid)
        + derivative_a * interval(a_lo - a_mid, a_hi - a_mid)
        + derivative_b * interval(b_lo - b_mid, b_hi - b_mid)
    )
    return upper_float(enclosure)


@lru_cache(maxsize=None)
def p_upper(
    a_lo: Fraction,
    a_hi: Fraction,
    b_lo: Fraction,
    b_hi: Fraction,
) -> float:
    left_lo = max(ZERO, b_lo - a_hi * HALF)
    left_hi = min(ONE, b_hi - a_lo * HALF)
    right_lo = max(ZERO, ONE - b_hi - a_hi * HALF)
    right_hi = min(ONE, ONE - b_lo - a_lo * HALF)
    if left_lo > left_hi or right_lo > right_hi:
        return -math.inf
    a_term = interval(a_lo, a_hi) * LN2
    _, b_upper = xlogx_range(b_lo, b_hi)
    _, one_b_upper = xlogx_range(ONE - b_hi, ONE - b_lo)
    left_lower, _ = xlogx_range(left_lo, left_hi)
    right_lower, _ = xlogx_range(right_lo, right_hi)
    _, one_a_upper = xlogx_range(ONE - a_hi, ONE - a_lo)
    range_upper = add_upper(
        upper_float(a_term),
        b_upper,
        one_b_upper,
        -left_lower,
        -right_lower,
        one_a_upper,
    )
    return min(
        range_upper,
        p_taylor_upper(a_lo, a_hi, b_lo, b_hi),
    )


def golay_mean_weight(log_r: float) -> float:
    terms = {
        0: 0.0,
        8: math.log(759.0) + 8 * log_r,
        12: math.log(2576.0) + 12 * log_r,
        16: math.log(759.0) + 16 * log_r,
        24: 24 * log_r,
    }
    scale = max(terms.values())
    scaled = {weight: math.exp(value - scale) for weight, value in terms.items()}
    denominator = sum(scaled.values())
    return sum(weight * value for weight, value in scaled.items()) / (
        24.0 * denominator
    )


@lru_cache(maxsize=None)
def choose_r(alpha: Fraction) -> Fraction:
    target = float(alpha)
    lo = -80.0
    hi = 80.0
    for _ in range(160):
        mid = (lo + hi) / 2
        if golay_mean_weight(mid) < target:
            lo = mid
        else:
            hi = mid
    return Fraction(format(math.exp((lo + hi) / 2), ".17g"))


@lru_cache(maxsize=None)
def g_upper(alpha_lo: Fraction, alpha_hi: Fraction) -> float:
    midpoint = (alpha_lo + alpha_hi) * HALF
    r = as_interval(choose_r(midpoint))
    polynomial = 1 + 759 * r**8 + 2576 * r**12 + 759 * r**16 + r**24
    objective = iv.log(polynomial) / 24 - interval(alpha_lo, alpha_hi) * iv.log(r)
    return upper_float(objective)


@dataclass(frozen=True)
class Segment:
    index: int
    omega_lo: Fraction
    omega_hi: Fraction
    slope: Fraction
    intercept: Fraction


def left_segments() -> list[Segment]:
    lines = LEFT_LINES + [(ZERO, Fraction(0))]
    # Replace the placeholder intercept on the last line by ln(2)/2 only
    # after its exact intersection with the preceding rational line has been
    # enclosed.  A rational upper endpoint is used for the interval search;
    # the central line itself is proved by the code-size bound.
    ln2_half_upper = Fraction(format(upper_float(LN2 / 2), ".17g"))
    lines[-1] = (ZERO, ln2_half_upper)
    boundaries = [LEFT_ENDPOINT]
    for (s0, b0), (s1, b1) in zip(lines, lines[1:], strict=False):
        boundary = (b1 - b0) / (s0 - s1)
        boundaries.append(boundary)
    if boundaries != sorted(boundaries) or boundaries[-1] >= HALF:
        raise AssertionError("support-line intersections are not ordered")
    result = []
    for index, ((slope, intercept), lo, hi) in enumerate(
        zip(lines[:-1], boundaries[:-1], boundaries[1:], strict=True)
    ):
        result.append(Segment(index, lo, hi, slope, intercept))
    return result


@dataclass(frozen=True)
class Box:
    segment: Segment
    alpha_lo: Fraction
    alpha_hi: Fraction
    beta_lo: Fraction
    beta_hi: Fraction
    omega_lo: Fraction
    omega_hi: Fraction

    def feasible(self) -> bool:
        if self.beta_hi < self.alpha_lo * HALF:
            return False
        if self.beta_lo > ONE - self.alpha_lo * HALF:
            return False
        if self.omega_hi < self.beta_lo * HALF:
            return False
        if self.omega_lo > ONE - self.beta_lo * HALF:
            return False
        return True

    def line_lower(self) -> float:
        slope = as_interval(self.segment.slope)
        intercept = as_interval(self.segment.intercept)
        omega = interval(self.omega_lo, self.omega_hi)
        return lower_float(slope * omega + intercept)

    def objective_upper(self) -> float:
        if not self.feasible():
            return -math.inf
        component = add_upper(
            g_upper(self.alpha_lo, self.alpha_hi),
            p_upper(self.alpha_lo, self.alpha_hi, self.beta_lo, self.beta_hi),
            p_upper(self.beta_lo, self.beta_hi, self.omega_lo, self.omega_hi),
            -self.line_lower(),
        )
        return min(component, self.taylor_upper())

    def taylor_upper(self) -> float:
        if not (
            self.beta_lo > self.alpha_hi * HALF
            and self.beta_hi < ONE - self.alpha_hi * HALF
            and self.omega_lo > self.beta_hi * HALF
            and self.omega_hi < ONE - self.beta_hi * HALF
            and self.alpha_hi < ONE
            and self.beta_lo > ZERO
            and self.beta_hi < ONE
        ):
            return math.inf
        aa = (self.alpha_lo + self.alpha_hi) * HALF
        bb = (self.beta_lo + self.beta_hi) * HALF
        ww = (self.omega_lo + self.omega_hi) * HALF
        r = as_interval(choose_r(aa))
        polynomial = 1 + 759 * r**8 + 2576 * r**12 + 759 * r**16 + r**24
        center = (
            iv.log(polynomial) / 24
            - as_interval(aa) * iv.log(r)
            + p_point(aa, bb)
            + p_point(bb, ww)
            - as_interval(self.segment.slope) * as_interval(ww)
            - as_interval(self.segment.intercept)
        )
        alpha = interval(self.alpha_lo, self.alpha_hi)
        beta = interval(self.beta_lo, self.beta_hi)
        omega = interval(self.omega_lo, self.omega_hi)
        p1a, p1b = p_gradient(alpha, beta)
        p2b, p2w = p_gradient(beta, omega)
        gradient_alpha = -iv.log(r) + p1a
        gradient_beta = p1b + p2b
        gradient_omega = p2w - as_interval(self.segment.slope)
        enclosure = (
            center
            + gradient_alpha
            * interval(self.alpha_lo - aa, self.alpha_hi - aa)
            + gradient_beta * interval(self.beta_lo - bb, self.beta_hi - bb)
            + gradient_omega
            * interval(self.omega_lo - ww, self.omega_hi - ww)
        )
        return upper_float(enclosure)

    def split(self) -> tuple["Box", "Box"]:
        widths = [
            self.alpha_hi - self.alpha_lo,
            self.beta_hi - self.beta_lo,
            2 * (self.omega_hi - self.omega_lo),
        ]
        dimension = max(range(3), key=lambda index: widths[index])
        bounds = [
            [self.alpha_lo, self.alpha_hi],
            [self.beta_lo, self.beta_hi],
            [self.omega_lo, self.omega_hi],
        ]
        midpoint = sum(bounds[dimension], ZERO) / 2
        left = [pair[:] for pair in bounds]
        right = [pair[:] for pair in bounds]
        left[dimension][1] = midpoint
        right[dimension][0] = midpoint
        return (
            Box(self.segment, *left[0], *left[1], *left[2]),
            Box(self.segment, *right[0], *right[1], *right[2]),
        )


def reflected(segment: Segment) -> Segment:
    # s(1-w)+b = -s w + (s+b).
    return Segment(
        -segment.index - 1,
        ONE - segment.omega_hi,
        ONE - segment.omega_lo,
        -segment.slope,
        segment.slope + segment.intercept,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-boxes", type=int, default=4_000_000)
    parser.add_argument("--max-depth", type=int, default=90)
    parser.add_argument("--progress-every", type=int, default=100_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_concave_majorant.json"
        ),
    )
    args = parser.parse_args()

    left = left_segments()
    segments = left + [reflected(segment) for segment in reversed(left)]
    heap = []
    serial = 0
    for segment in segments:
        box = Box(segment, ZERO, ONE, ZERO, ONE, segment.omega_lo, segment.omega_hi)
        heapq.heappush(heap, (-box.objective_upper(), 0, serial, box))
        serial += 1

    processed = 0
    accepted = 0
    deepest = 0
    max_accepted = -math.inf
    unresolved = None
    while heap and processed < args.max_boxes:
        negative_upper, depth, _, box = heapq.heappop(heap)
        upper = -negative_upper
        processed += 1
        deepest = max(deepest, depth)
        if upper < 0:
            accepted += 1
            max_accepted = max(max_accepted, upper)
            continue
        if depth >= args.max_depth:
            unresolved = (upper, depth, box)
            break
        for child in box.split():
            child_upper = child.objective_upper()
            if child_upper < 0:
                accepted += 1
                max_accepted = max(max_accepted, child_upper)
            else:
                heapq.heappush(heap, (-child_upper, depth + 1, serial, child))
                serial += 1
        if args.progress_every and processed % args.progress_every == 0:
            print(
                json.dumps(
                    {
                        "event": "progress",
                        "processed": processed,
                        "queued": len(heap),
                        "accepted": accepted,
                        "next_upper": -heap[0][0] if heap else None,
                        "deepest": deepest,
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    if unresolved is None and heap:
        upper, depth, box = -heap[0][0], heap[0][1], heap[0][3]
        unresolved = (upper, depth, box)

    ln2_half_upper = Fraction(format(upper_float(LN2 / 2), ".17g"))
    payload = {
        "schema": "golay-ba3-concave-majorant-v1",
        "status": "proved" if unresolved is None else "unresolved",
        "claim": (
            "The listed lower envelope of affine lines majorizes the limiting "
            "expected Golay--BA-3 spectrum exponent on [0.104,0.896]."
        ),
        "interval_dps": iv.dps,
        "weight_interval": [str(LEFT_ENDPOINT), str(RIGHT_ENDPOINT)],
        "central_constant_upper": str(ln2_half_upper),
        "central_constant_reason": "At most 2^(B/2) codewords exist.",
        "left_segments": [
            {
                "weight": [str(item.omega_lo), str(item.omega_hi)],
                "slope": str(item.slope),
                "intercept": str(item.intercept),
            }
            for item in left
        ],
        "right_segments_are_reflections": True,
        "processed_boxes": processed,
        "accepted_boxes": accepted,
        "deepest_box": deepest,
        "largest_accepted_upper_natural": max_accepted,
        "unresolved": None
        if unresolved is None
        else {
            "upper": unresolved[0],
            "depth": unresolved[1],
            "segment": unresolved[2].segment.index,
            "alpha": [str(unresolved[2].alpha_lo), str(unresolved[2].alpha_hi)],
            "beta": [str(unresolved[2].beta_lo), str(unresolved[2].beta_hi)],
            "omega": [str(unresolved[2].omega_lo), str(unresolved[2].omega_hi)],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "left_segments"}, indent=2, sort_keys=True))
    if unresolved is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
