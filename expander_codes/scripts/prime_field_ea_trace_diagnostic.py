#!/usr/bin/env python3
"""Trace-placement diagnostic for labeled prime-field Bernoulli EA.

The transfer counts occupied expander coordinates and exact zero/nonzero
accumulator traces.  A marker records independent zero-return equations.  The
projective union pays only equations beyond the message-support dimension.
This is a floating-point saddle diagnostic, not an interval certificate.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from expander_bounds import log_binom, log_weight_mgf
from prime_field_ea_diagnostic import gv_distance
from regular_ec_diagnostic import diagnostic_weights


@dataclass(frozen=True)
class ProjectiveTraceTerm:
    message_weight: int
    occupancy_probability: float
    structural_log2: float
    field_log2: float
    total_log2: float
    structural_z: float
    structural_v: float
    field_z: float
    field_v: float


def occupancy_probability(theta: float, message_weight: int) -> float:
    return -math.expm1(message_weight * math.log1p(-theta))


def empty_probability(theta: float, message_weight: int) -> float:
    return math.exp(message_weight * math.log1p(-theta))


def trace_matrix(occupied_probability: float, z: float, v: float) -> np.ndarray:
    """Transfer for output weight and zero-state occupied coordinates."""
    q = occupied_probability
    return np.array(
        ((1.0 - q + q * v, q * z), (q * v, z)),
        dtype=np.float64,
    )


def trace_matrix_from_empty(empty: float, z: float, v: float) -> np.ndarray:
    """Equivalent transfer that preserves a tiny empty-coordinate mass."""
    occupied = 1.0 - empty
    return np.array(
        ((empty + occupied * v, occupied * z), (occupied * v, z)),
        dtype=np.float64,
    )


def structural_trace_log_bound(
    *, n: int, cutoff: int, empty: float, message_weight: int
) -> tuple[float, float, float]:
    """Bound traces with at most r-1 zero-return equations."""
    threshold = message_weight - 1
    log_scale = min(
        700.0,
        max(120.0, (-math.log(empty) if empty else 700.0) + 100.0),
    )

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        z = math.exp(log_z)
        v = math.exp(log_v)
        return (
            log_weight_mgf(trace_matrix_from_empty(empty, z, v), n, 0)
            - cutoff * log_z
            - threshold * log_v
        )

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    expected_events = (1.0 - empty) * n
    v_ratio = message_weight / max(1.0, expected_events)
    starts = (
        (z_start, max(1e-12, min(0.8, v_ratio))),
        (z_start, 0.5),
        (z_start, 1e-6),
    )
    candidates = []
    for start_z, start_v in starts:
        candidates.append(minimize(
            objective,
            np.log(np.array([start_z, start_v], dtype=np.float64)),
            method="L-BFGS-B",
            bounds=((-log_scale, -1e-11), (-log_scale, -1e-11)),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
        ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    return (
        float(optimum.fun),
        math.exp(float(optimum.x[0])),
        math.exp(float(optimum.x[1])),
    )


def field_trace_log_bound(
    *, prime: int, n: int, cutoff: int, empty: float,
    message_weight: int,
) -> tuple[float, float, float]:
    """Bound traces with at least r equations and their surplus field weight."""
    log_s = math.log(prime - 1)
    log_scale = min(
        700.0,
        max(120.0, (-math.log(empty) if empty else 700.0) + 100.0),
    )

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        z = math.exp(log_z)
        v = math.exp(log_v)
        return (
            -log_s
            - message_weight * log_v
            + log_weight_mgf(trace_matrix_from_empty(empty, z, v), n, 0)
            - cutoff * log_z
        )

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    expected_events = (1.0 - empty) * n
    threshold_ratio = message_weight / max(1.0, expected_events)
    minimum_v = math.exp(-log_s)
    starts = (
        (z_start, max(minimum_v, min(0.8, threshold_ratio))),
        (z_start, max(minimum_v, 0.5)),
        (z_start, minimum_v),
    )
    candidates = []
    for start_z, start_v in starts:
        candidates.append(minimize(
            objective,
            np.log(np.array([start_z, start_v], dtype=np.float64)),
            method="L-BFGS-B",
            bounds=((-log_scale, -1e-11), (-log_s, -1e-11)),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
        ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    return (
        float(optimum.fun),
        math.exp(float(optimum.x[0])),
        math.exp(float(optimum.x[1])),
    )


def projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, degree: float,
    message_weight: int,
) -> ProjectiveTraceTerm:
    theta = degree / n
    empty = empty_probability(theta, message_weight)
    q = -math.expm1(math.log(empty)) if empty else 1.0
    structural, structural_z, structural_v = structural_trace_log_bound(
        n=n,
        cutoff=cutoff,
        empty=empty,
        message_weight=message_weight,
    )
    field, field_z, field_v = field_trace_log_bound(
        prime=prime,
        n=n,
        cutoff=cutoff,
        empty=empty,
        message_weight=message_weight,
    )
    support_count = log_binom(k, message_weight)
    structural += support_count
    field += support_count
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return ProjectiveTraceTerm(
        message_weight=message_weight,
        occupancy_probability=q,
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
        structural_z=structural_z,
        structural_v=structural_v,
        field_z=field_z,
        field_v=field_v,
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
    parser.add_argument("--summary-only", action="store_true")
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
        term = projective_trace_logterm(
            prime=args.prime,
            k=k,
            n=args.n,
            cutoff=cutoff,
            degree=args.degree,
            message_weight=weight,
        )
        terms.append(term)
        if not args.summary_only:
            print(
                f"r={weight}\tq={term.occupancy_probability:.8g}\t"
                f"struct={term.structural_log2:.6f}\t"
                f"field={term.field_log2:.6f}\t"
                f"total={term.total_log2:.6f}\t"
                f"z0={term.structural_z:.7g}\tv={term.structural_v:.7g}\t"
                f"z1={term.field_z:.7g}\tu={term.field_v:.7g}"
            )
    worst = max(terms, key=lambda item: item.total_log2)
    print(
        f"sampled worst: r={worst.message_weight}, "
        f"log2_bound={worst.total_log2:.6f}"
    )


if __name__ == "__main__":
    main()
