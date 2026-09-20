"""Outward-rounded branch verifier for the Golay--BA-3 outer exponent.

The verifier proves the central likelihood inequality (17) in
LINEAR_OUTER_CERTIFICATE.md.  Every search box has rational endpoints.  The
transcendental point evaluations use mpmath interval arithmetic, and every
conversion to binary64 is expanded by one representable number.

This file does not verify the two endpoint sums in equation (16).
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from mpmath import iv


iv.dps = 60
ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)
LN2_INTERVAL = iv.log(iv.mpf(2))


def as_interval(value: Fraction):
    return iv.mpf(value.numerator) / value.denominator


def interval(lo: Fraction, hi: Fraction):
    return iv.mpf([as_interval(lo).a, as_interval(hi).b])


def lower_float(value) -> float:
    return math.nextafter(float(value.a), -math.inf)


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


@lru_cache(maxsize=None)
def xlogx_point(value: Fraction) -> tuple[float, float]:
    if value == 0:
        return 0.0, 0.0
    x = as_interval(value)
    result = x * iv.log(x)
    return lower_float(result), upper_float(result)


def xlogx_range(lo: Fraction, hi: Fraction) -> tuple[float, float]:
    """Return an enclosure for x log x on [lo,hi] within [0,1]."""
    if lo < 0 or hi > 1 or lo > hi:
        raise ValueError(f"invalid x log x interval [{lo}, {hi}]")
    lo_bounds = xlogx_point(lo)
    hi_bounds = xlogx_point(hi)
    upper = math.nextafter(max(lo_bounds[1], hi_bounds[1]), math.inf)

    inv_e = iv.mpf(1) / iv.e
    inv_e_lo = Fraction(str(lower_float(inv_e)))
    inv_e_hi = Fraction(str(upper_float(inv_e)))
    if lo <= inv_e_lo and inv_e_hi <= hi:
        minimum = -iv.mpf(1) / iv.e
        lower = lower_float(minimum)
    else:
        lower = math.nextafter(min(lo_bounds[0], hi_bounds[0]), -math.inf)
    return lower, upper


def add_upper(*values: float) -> float:
    return math.nextafter(math.fsum(values), math.inf)


def xlogx_interval(value):
    return value * iv.log(value)


def p_point(a: Fraction, b: Fraction):
    a_iv = as_interval(a)
    b_iv = as_interval(b)
    return (
        a_iv * LN2_INTERVAL
        + xlogx_interval(b_iv)
        + xlogx_interval(1 - b_iv)
        - xlogx_interval(b_iv - a_iv / 2)
        - xlogx_interval(1 - b_iv - a_iv / 2)
        + xlogx_interval(1 - a_iv)
    )


def p_taylor_upper(
    a_lo: Fraction,
    a_hi: Fraction,
    b_lo: Fraction,
    b_hi: Fraction,
) -> float:
    """Mean-value enclosure that retains the cancellations in p(a,b)."""
    x_lo = b_lo - a_hi * HALF
    y_lo = ONE - b_hi - a_hi * HALF
    if (
        x_lo <= 0
        or y_lo <= 0
        or a_hi >= 1
        or b_lo <= 0
        or b_hi >= 1
    ):
        return math.inf

    a_mid = (a_lo + a_hi) * HALF
    b_mid = (b_lo + b_hi) * HALF
    if not (a_mid * HALF < b_mid < ONE - a_mid * HALF):
        return math.inf

    a = interval(a_lo, a_hi)
    b = interval(b_lo, b_hi)
    x = b - a / 2
    y = 1 - b - a / 2
    z = 1 - a
    derivative_a = iv.log(2 * iv.sqrt(x * y) / z)
    derivative_b = iv.log(b * y / ((1 - b) * x))
    delta_a = interval(a_lo - a_mid, a_hi - a_mid)
    delta_b = interval(b_lo - b_mid, b_hi - b_mid)
    enclosure = (
        p_point(a_mid, b_mid)
        + derivative_a * delta_a
        + derivative_b * delta_b
    )
    return upper_float(enclosure)


def p_gradient(a, b):
    x = b - a / 2
    y = 1 - b - a / 2
    derivative_a = iv.log(2 * iv.sqrt(x * y) / (1 - a))
    derivative_b = iv.log(b * y / ((1 - b) * x))
    return derivative_a, derivative_b


def p_upper(
    a_lo: Fraction,
    a_hi: Fraction,
    b_lo: Fraction,
    b_hi: Fraction,
) -> float:
    """Upper-bound the accumulator exponent p(a,b) on a feasible box."""
    if not (ZERO <= a_lo <= a_hi <= ONE and ZERO <= b_lo <= b_hi <= ONE):
        raise ValueError("accumulator box leaves [0,1]^2")

    b_minus_a_lo = max(ZERO, b_lo - a_hi * HALF)
    b_minus_a_hi = min(ONE, b_hi - a_lo * HALF)
    other_lo = max(ZERO, ONE - b_hi - a_hi * HALF)
    other_hi = min(ONE, ONE - b_lo - a_lo * HALF)
    if b_minus_a_lo > b_minus_a_hi or other_lo > other_hi:
        return -math.inf

    a_term = interval(a_lo, a_hi) * LN2_INTERVAL
    _, f_b_upper = xlogx_range(b_lo, b_hi)
    _, f_one_b_upper = xlogx_range(ONE - b_hi, ONE - b_lo)
    f_ba_lower, _ = xlogx_range(b_minus_a_lo, b_minus_a_hi)
    f_other_lower, _ = xlogx_range(other_lo, other_hi)
    _, f_one_a_upper = xlogx_range(ONE - a_hi, ONE - a_lo)
    range_upper = add_upper(
        upper_float(a_term),
        f_b_upper,
        f_one_b_upper,
        -f_ba_lower,
        -f_other_lower,
        f_one_a_upper,
    )
    return min(
        range_upper,
        p_taylor_upper(a_lo, a_hi, b_lo, b_hi),
    )


def negative_entropy_upper(lo: Fraction, hi: Fraction) -> float:
    """Upper-bound -h(x), a convex function, by its endpoint values."""
    values = []
    for endpoint in (lo, hi):
        x = as_interval(endpoint)
        if endpoint == 0 or endpoint == 1:
            values.append(iv.mpf(0))
        else:
            values.append(xlogx_interval(x) + xlogx_interval(1 - x))
    return math.nextafter(max(upper_float(value) for value in values), math.inf)


def golay_mean_weight_from_log(log_r: float) -> float:
    powers = {
        0: 1.0,
        8: 759.0,
        12: 2576.0,
        16: 759.0,
        24: 1.0,
    }
    logs = {
        weight: math.log(count) + weight * log_r
        for weight, count in powers.items()
    }
    scale = max(logs.values())
    scaled = {
        weight: math.exp(value - scale) for weight, value in logs.items()
    }
    denominator = sum(scaled.values())
    numerator = sum(weight * value for weight, value in scaled.items())
    return numerator / (24.0 * denominator)


def choose_r(alpha: Fraction) -> Fraction:
    target = float(alpha)
    low_log = -80.0
    high_log = 80.0
    for _ in range(160):
        middle_log = 0.5 * (low_log + high_log)
        if golay_mean_weight_from_log(middle_log) < target:
            low_log = middle_log
        else:
            high_log = middle_log
    chosen = math.exp(0.5 * (low_log + high_log))
    return Fraction(format(chosen, ".17g"))


def g_upper(alpha_lo: Fraction, alpha_hi: Fraction) -> float:
    """Upper-bound g(alpha) by evaluating its infimum objective at one r."""
    midpoint = (alpha_lo + alpha_hi) * HALF
    r_fraction = choose_r(midpoint)
    r = as_interval(r_fraction)
    polynomial = (
        1
        + 759 * r**8
        + 2576 * r**12
        + 759 * r**16
        + r**24
    )
    alpha = interval(alpha_lo, alpha_hi)
    objective = iv.log(polynomial) / 24 - alpha * iv.log(r)
    return upper_float(objective)


@dataclass(frozen=True)
class Box:
    alpha_lo: Fraction
    alpha_hi: Fraction
    beta_lo: Fraction
    beta_hi: Fraction
    omega_lo: Fraction
    omega_hi: Fraction

    def intersects_feasible_region(self) -> bool:
        if self.beta_hi < self.alpha_lo * HALF:
            return False
        if self.beta_lo > ONE - self.alpha_lo * HALF:
            return False
        if self.omega_hi < self.beta_lo * HALF:
            return False
        if self.omega_lo > ONE - self.beta_lo * HALF:
            return False
        return True

    def objective_upper(self, include_likelihood: bool) -> float:
        if not self.intersects_feasible_region():
            return -math.inf
        first = g_upper(self.alpha_lo, self.alpha_hi)
        second = p_upper(
            self.alpha_lo, self.alpha_hi, self.beta_lo, self.beta_hi
        )
        third = p_upper(
            self.beta_lo, self.beta_hi, self.omega_lo, self.omega_hi
        )
        terms = [
            first,
            second,
            third,
        ]
        if include_likelihood:
            terms.extend(
                [
                    negative_entropy_upper(self.omega_lo, self.omega_hi),
                    upper_float(LN2_INTERVAL),
                ]
            )
        component_upper = add_upper(*terms)
        return min(
            component_upper,
            self.objective_taylor_upper(include_likelihood),
        )

    def objective_taylor_upper(self, include_likelihood: bool) -> float:
        """Enclose the complete objective so derivative cancellations remain."""
        if not (
            self.beta_lo > self.alpha_hi * HALF
            and self.beta_hi < ONE - self.alpha_hi * HALF
            and self.omega_lo > self.beta_hi * HALF
            and self.omega_hi < ONE - self.beta_hi * HALF
            and self.alpha_hi < ONE
            and self.beta_lo > ZERO
            and self.beta_hi < ONE
            and self.omega_lo > ZERO
            and self.omega_hi < ONE
        ):
            return math.inf

        alpha_mid = (self.alpha_lo + self.alpha_hi) * HALF
        beta_mid = (self.beta_lo + self.beta_hi) * HALF
        omega_mid = (self.omega_lo + self.omega_hi) * HALF
        r_fraction = choose_r(alpha_mid)
        r = as_interval(r_fraction)

        polynomial = (
            1
            + 759 * r**8
            + 2576 * r**12
            + 759 * r**16
            + r**24
        )
        alpha_point = as_interval(alpha_mid)
        omega_point = as_interval(omega_mid)
        center = (
            iv.log(polynomial) / 24
            - alpha_point * iv.log(r)
            + p_point(alpha_mid, beta_mid)
            + p_point(beta_mid, omega_mid)
        )
        if include_likelihood:
            center += (
                xlogx_interval(omega_point)
                + xlogx_interval(1 - omega_point)
                + LN2_INTERVAL
            )

        alpha = interval(self.alpha_lo, self.alpha_hi)
        beta = interval(self.beta_lo, self.beta_hi)
        omega = interval(self.omega_lo, self.omega_hi)
        p1_alpha, p1_beta = p_gradient(alpha, beta)
        p2_beta, p2_omega = p_gradient(beta, omega)
        gradient_alpha = -iv.log(r) + p1_alpha
        gradient_beta = p1_beta + p2_beta
        gradient_omega = p2_omega
        if include_likelihood:
            gradient_omega += iv.log(omega / (1 - omega))

        alpha_radius = (self.alpha_hi - self.alpha_lo) * HALF
        beta_radius = (self.beta_hi - self.beta_lo) * HALF
        omega_radius = (self.omega_hi - self.omega_lo) * HALF
        delta_alpha = interval(-alpha_radius, alpha_radius)
        delta_beta = interval(-beta_radius, beta_radius)
        delta_omega = interval(-omega_radius, omega_radius)
        enclosure = (
            center
            + gradient_alpha * delta_alpha
            + gradient_beta * delta_beta
            + gradient_omega * delta_omega
        )
        return upper_float(enclosure)

    def normalized_widths(self) -> tuple[Fraction, Fraction, Fraction]:
        return (
            self.alpha_hi - self.alpha_lo,
            self.beta_hi - self.beta_lo,
            (self.omega_hi - self.omega_lo) / Fraction(99, 125),
        )

    def split(self) -> tuple["Box", "Box"]:
        widths = self.normalized_widths()
        dimension = max(range(3), key=lambda index: widths[index])
        if dimension == 0:
            midpoint = (self.alpha_lo + self.alpha_hi) * HALF
            return (
                Box(
                    self.alpha_lo,
                    midpoint,
                    self.beta_lo,
                    self.beta_hi,
                    self.omega_lo,
                    self.omega_hi,
                ),
                Box(
                    midpoint,
                    self.alpha_hi,
                    self.beta_lo,
                    self.beta_hi,
                    self.omega_lo,
                    self.omega_hi,
                ),
            )
        if dimension == 1:
            midpoint = (self.beta_lo + self.beta_hi) * HALF
            return (
                Box(
                    self.alpha_lo,
                    self.alpha_hi,
                    self.beta_lo,
                    midpoint,
                    self.omega_lo,
                    self.omega_hi,
                ),
                Box(
                    self.alpha_lo,
                    self.alpha_hi,
                    midpoint,
                    self.beta_hi,
                    self.omega_lo,
                    self.omega_hi,
                ),
            )
        midpoint = (self.omega_lo + self.omega_hi) * HALF
        return (
            Box(
                self.alpha_lo,
                self.alpha_hi,
                self.beta_lo,
                self.beta_hi,
                self.omega_lo,
                midpoint,
            ),
            Box(
                self.alpha_lo,
                self.alpha_hi,
                self.beta_lo,
                self.beta_hi,
                midpoint,
                self.omega_hi,
            ),
        )

    def as_json(self) -> dict[str, list[str]]:
        return {
            "alpha": [str(self.alpha_lo), str(self.alpha_hi)],
            "beta": [str(self.beta_lo), str(self.beta_hi)],
            "omega": [str(self.omega_lo), str(self.omega_hi)],
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["central", "low-tail"], default="central")
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--max-boxes", type=int, default=2_000_000)
    parser.add_argument("--max-depth", type=int, default=80)
    parser.add_argument("--progress-every", type=int, default=100_000)
    args = parser.parse_args()

    if args.mode == "central":
        initial = Box(
            ZERO, ONE, ZERO, ONE, Fraction(13, 125), Fraction(112, 125)
        )
        threshold = 0.36 if args.threshold is None else args.threshold
        include_likelihood = True
    else:
        initial = Box(
            Fraction(1, 125),
            ONE,
            ZERO,
            ONE,
            Fraction(1, 500),
            Fraction(13, 125),
        )
        threshold = 0.0 if args.threshold is None else args.threshold
        include_likelihood = False
    initial_upper = initial.objective_upper(include_likelihood)
    heap: list[tuple[float, int, int, Box]] = [(-initial_upper, 0, 0, initial)]
    serial = 1
    processed = 0
    accepted = 0
    max_accepted_upper = -math.inf
    deepest = 0
    unresolved: tuple[float, int, Box] | None = None

    while heap and processed < args.max_boxes:
        negative_upper, depth, _, box = heapq.heappop(heap)
        upper = -negative_upper
        processed += 1
        deepest = max(deepest, depth)
        if upper < threshold:
            accepted += 1
            max_accepted_upper = max(max_accepted_upper, upper)
            continue
        if depth >= args.max_depth:
            unresolved = (upper, depth, box)
            break
        for child in box.split():
            child_upper = child.objective_upper(include_likelihood)
            if child_upper < threshold:
                accepted += 1
                max_accepted_upper = max(max_accepted_upper, child_upper)
                continue
            heapq.heappush(heap, (-child_upper, depth + 1, serial, child))
            serial += 1
        if args.progress_every and processed % args.progress_every == 0:
            next_upper = -heap[0][0] if heap else -math.inf
            print(
                json.dumps(
                    {
                        "event": "progress",
                        "processed": processed,
                        "queued": len(heap),
                        "accepted": accepted,
                        "max_accepted_upper": max_accepted_upper,
                        "next_upper": next_upper,
                        "deepest": deepest,
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    if unresolved is None and heap:
        upper, depth, box = -heap[0][0], heap[0][1], heap[0][3]
        unresolved = (upper, depth, box)

    if unresolved is None:
        result = {
            "status": "proved",
            "claim": (
                "central Golay--BA-3 likelihood exponent is below threshold"
                if args.mode == "central"
                else "Golay--BA-3 low-tail expected-spectrum exponent is below threshold"
            ),
            "mode": args.mode,
            "threshold": threshold,
            "processed": processed,
            "accepted": accepted,
            "max_accepted_upper": max_accepted_upper,
            "deepest": deepest,
            "remaining_boxes": 0,
            "interval_dps": iv.dps,
        }
    else:
        upper, depth, box = unresolved
        result = {
            "status": "unresolved",
            "mode": args.mode,
            "threshold": threshold,
            "processed": processed,
            "accepted": accepted,
            "max_accepted_upper": max_accepted_upper,
            "deepest": deepest,
            "remaining_boxes": len(heap),
            "worst_upper": upper,
            "worst_depth": depth,
            "worst_box": box.as_json(),
            "interval_dps": iv.dps,
        }
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "proved":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
