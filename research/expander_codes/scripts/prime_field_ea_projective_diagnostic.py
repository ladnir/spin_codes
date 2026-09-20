#!/usr/bin/env python3
"""Saddle diagnostic for the Bernoulli EA projective gap enumerator.

The proved finite enumerator separates the structural event with fewer than
``r`` hit zero-gaps from the part suppressed by surplus edge-label equations.
This script applies positive-coefficient Chernoff bounds to those two pieces.
It is a floating-point diagnostic, not an interval certificate.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from expander_bounds import log_binom
from prime_field_ea_diagnostic import gv_distance
from regular_ec_diagnostic import diagnostic_weights


@dataclass(frozen=True)
class ProjectiveGapTerm:
    message_weight: int
    structural_log2: float
    field_log2: float
    total_log2: float
    structural_z: float
    structural_v: float
    field_z: float


def log_gap_factor(theta: float, weight: int, log_z: float, log_v: float) -> float:
    """Return log F_r(z,v) without forming small powers directly."""
    z = math.exp(log_z)
    v = math.exp(log_v)
    log_a = weight * math.log1p(-theta)
    a_z = math.exp(log_a + log_z)
    log_hit = log_v - math.log1p(-z)
    if v == 1.0:
        return log_hit
    log_miss = math.log1p(-v) + log_a - math.log1p(-a_z)
    return float(np.logaddexp(log_hit, log_miss))


def structural_log_bound(
    *, n: int, cutoff: int, theta: float, message_weight: int
) -> tuple[float, float, float]:
    """Bound the summed mass of containers with hit-gap count below r."""
    zero_count = n - cutoff
    threshold = message_weight - 1

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return (
            -cutoff * log_z
            - threshold * log_v
            + zero_count * log_gap_factor(
                theta, message_weight, log_z, log_v
            )
            - math.log1p(-math.exp(log_z))
        )

    z_start = max(1e-8, min(0.95, cutoff / (n + 1.0)))
    expected_edges = theta * message_weight * n
    v_start = max(1e-8, min(0.8, message_weight / max(1.0, expected_edges)))
    optimum = minimize(
        objective,
        np.log(np.array([z_start, v_start], dtype=np.float64)),
        method="L-BFGS-B",
        bounds=((-50.0, -1e-10), (-80.0, -1e-10)),
        options={"ftol": 1e-12, "gtol": 1e-9, "maxiter": 500},
    )
    return float(optimum.fun), math.exp(float(optimum.x[0])), math.exp(float(optimum.x[1]))


def field_log_bound(
    *, prime: int, n: int, cutoff: int, theta: float, message_weight: int
) -> tuple[float, float]:
    """Bound the gap enumerator weighted by (p-1)^(r-1-c)."""
    zero_count = n - cutoff
    log_s = math.log(prime - 1)
    log_v = -log_s

    def objective(log_z: float) -> float:
        return (
            (message_weight - 1) * log_s
            - cutoff * log_z
            + zero_count * log_gap_factor(
                theta, message_weight, log_z, log_v
            )
            - math.log1p(-math.exp(log_z))
        )

    optimum = minimize_scalar(
        objective,
        bounds=(-50.0, -1e-10),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 300},
    )
    return float(optimum.fun), math.exp(float(optimum.x))


def logaddexp(left: float, right: float) -> float:
    return float(np.logaddexp(left, right))


def projective_gap_logterm(
    *, prime: int, k: int, n: int, cutoff: int, degree: float,
    message_weight: int,
) -> ProjectiveGapTerm:
    theta = degree / n
    structural, structural_z, structural_v = structural_log_bound(
        n=n, cutoff=cutoff, theta=theta, message_weight=message_weight
    )
    field, field_z = field_log_bound(
        prime=prime,
        n=n,
        cutoff=cutoff,
        theta=theta,
        message_weight=message_weight,
    )
    support_count = log_binom(k, message_weight)
    structural += support_count
    field += support_count
    total = logaddexp(structural, field)
    inverse_log_two = 1.0 / math.log(2.0)
    return ProjectiveGapTerm(
        message_weight=message_weight,
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
        structural_z=structural_z,
        structural_v=structural_v,
        field_z=field_z,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime", type=int, default=170141183460469231731687303715884105727
    )
    parser.add_argument("--n", type=int, default=2_097_100)
    parser.add_argument("--k", type=int)
    parser.add_argument("--degree", type=float, default=100.0)
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--gv-fraction", type=float, default=0.995)
    parser.add_argument("--weights", type=int, nargs="+")
    args = parser.parse_args()
    k = args.k if args.k is not None else args.n // 2
    delta_gv = gv_distance(args.prime, k / args.n)
    cutoff = (
        args.cutoff
        if args.cutoff is not None
        else math.floor(args.gv_fraction * delta_gv * args.n)
    )
    weights = args.weights if args.weights else diagnostic_weights(k)
    print(
        f"p={args.prime} k={k} n={args.n} rate={k/args.n:.9f} "
        f"degree={args.degree:.6f} cutoff={cutoff} "
        f"relative_cutoff={cutoff/args.n:.12f} gv={delta_gv:.12f}"
    )
    terms = []
    for weight in weights:
        term = projective_gap_logterm(
            prime=args.prime,
            k=k,
            n=args.n,
            cutoff=cutoff,
            degree=args.degree,
            message_weight=weight,
        )
        terms.append(term)
        print(
            f"r={weight}\tstruct={term.structural_log2:.6f}\t"
            f"field={term.field_log2:.6f}\ttotal={term.total_log2:.6f}\t"
            f"z0={term.structural_z:.7g}\tv={term.structural_v:.7g}\t"
            f"z1={term.field_z:.7g}"
        )
    worst = max(terms, key=lambda item: item.total_log2)
    print(
        f"sampled worst: r={worst.message_weight}, "
        f"log2_bound={worst.total_log2:.6f}"
    )


if __name__ == "__main__":
    main()
