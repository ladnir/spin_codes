#!/usr/bin/env python3
"""Projective trace diagnostic for labeled left-regular prime-field EA.

The dense branch dominates each exact regional empty set by the probability
that all of its coordinates are missed.  This avoids the Poisson conditioning
penalty.  The result is a floating-point diagnostic, not a certificate.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from expander_bounds import log_binom, log_weight_mgf
from prime_field_ea_diagnostic import gv_distance
from prime_field_ea_trace_diagnostic import (
    field_trace_log_bound as poisson_field_trace_log_bound,
    structural_trace_log_bound as poisson_structural_trace_log_bound,
)
from regular_ec_diagnostic import (
    combine_uniform_slice_transfers,
    diagnostic_weights,
    log_poisson_conditioning_penalty,
)


@dataclass(frozen=True)
class RegularProjectiveTraceTerm:
    message_weight: int
    empty_set_weight: float
    structural_log2: float
    field_log2: float
    total_log2: float


def empty_set_weight(message_weight: int, region_length: int) -> float:
    """Per-coordinate weight in the factorized regular empty-set bound."""
    return math.exp(message_weight * math.log1p(-1.0 / region_length))


def dominated_trace_matrix(beta: float, z: float, v: float) -> np.ndarray:
    """Return beta times the empty transfer plus the occupied trace transfer."""
    return np.array(
        ((beta + v, z), (v, (1.0 + beta) * z)),
        dtype=np.float64,
    )


def _optimize_two_markers(
    *, objective, z_start: float, v_start: float, lower_z: float,
    lower_v: float,
) -> tuple[float, float, float]:
    starts = (
        (z_start, v_start),
        (z_start, 0.5),
        (z_start, max(math.exp(lower_v), 1e-12)),
    )
    candidates = []
    for start_z, start_v in starts:
        candidates.append(minimize(
            objective,
            np.log(np.array([start_z, start_v], dtype=np.float64)),
            method="L-BFGS-B",
            bounds=((lower_z, -1e-11), (lower_v, -1e-11)),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
        ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    return (
        float(optimum.fun),
        math.exp(float(optimum.x[0])),
        math.exp(float(optimum.x[1])),
    )


def structural_log_bound(
    *, n: int, cutoff: int, beta: float, message_weight: int
) -> tuple[float, float, float]:
    threshold = message_weight - 1
    scale = min(700.0, max(120.0, -math.log(beta) + 100.0))

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        z = math.exp(log_z)
        v = math.exp(log_v)
        return (
            log_weight_mgf(dominated_trace_matrix(beta, z, v), n, 0)
            - cutoff * log_z
            - threshold * log_v
        )

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    v_start = max(1e-12, min(0.8, message_weight / max(1.0, n)))
    return _optimize_two_markers(
        objective=objective,
        z_start=z_start,
        v_start=v_start,
        lower_z=-scale,
        lower_v=-scale,
    )


def field_log_bound(
    *, prime: int, n: int, cutoff: int, beta: float, message_weight: int
) -> tuple[float, float, float]:
    log_s = math.log(prime - 1)
    scale = min(700.0, max(120.0, -math.log(beta) + 100.0))

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        z = math.exp(log_z)
        v = math.exp(log_v)
        return (
            -log_s
            - message_weight * log_v
            + log_weight_mgf(dominated_trace_matrix(beta, z, v), n, 0)
            - cutoff * log_z
        )

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    v_start = max(math.exp(-log_s), min(0.8, message_weight / max(1.0, n)))
    return _optimize_two_markers(
        objective=objective,
        z_start=z_start,
        v_start=v_start,
        lower_z=-scale,
        lower_v=-log_s,
    )


def regular_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    beta = empty_set_weight(message_weight, region_length)
    structural, _, _ = structural_log_bound(
        n=n, cutoff=cutoff, beta=beta, message_weight=message_weight
    )
    field, _, _ = field_log_bound(
        prime=prime,
        n=n,
        cutoff=cutoff,
        beta=beta,
        message_weight=message_weight,
    )
    support_count = log_binom(k, message_weight)
    structural += support_count
    field += support_count
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return RegularProjectiveTraceTerm(
        message_weight=message_weight,
        empty_set_weight=beta,
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
    )


def poisson_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    """Poissonized regional trace bound, including its conditioning cost."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    empty = math.exp(-message_weight / region_length)
    structural, _, _ = poisson_structural_trace_log_bound(
        n=n, cutoff=cutoff, empty=empty, message_weight=message_weight
    )
    field, _, _ = poisson_field_trace_log_bound(
        prime=prime,
        n=n,
        cutoff=cutoff,
        empty=empty,
        message_weight=message_weight,
    )
    common = (
        log_binom(k, message_weight)
        + region_count * log_poisson_conditioning_penalty(message_weight)
    )
    structural += common
    field += common
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return RegularProjectiveTraceTerm(
        message_weight=message_weight,
        empty_set_weight=empty,
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
    )


def saddle_trace_matrix(
    endpoint_marker: float, z: float, v: float
) -> np.ndarray:
    """Return Q0 + (exp(t)-1) Q1 for the regional endpoint EGF."""
    empty, occupied = trace_input_matrices(z, v)
    return empty + math.expm1(endpoint_marker) * occupied


def _optimize_three_markers(
    *, objective, z_start: float, v_start: float, t_start: float,
    lower_z: float, lower_v: float, lower_t: float = -30.0,
    upper_t: float = math.log(700.0),
) -> tuple[float, float, float, float]:
    t_start = max(math.exp(lower_t), min(math.exp(upper_t), t_start))
    starts = (
        (z_start, v_start, t_start),
        (z_start, v_start, max(1e-12, 0.8 * t_start)),
        (z_start, v_start, min(math.exp(upper_t), 1.2 * t_start)),
    )
    candidates = []
    for start_z, start_v, start_t in starts:
        candidates.append(minimize(
            objective,
            np.log(np.array([start_z, start_v, start_t], dtype=np.float64)),
            method="L-BFGS-B",
            bounds=(
                (lower_z, -1e-11),
                (lower_v, -1e-11),
                (lower_t, upper_t),
            ),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 700},
        ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    return (
        float(optimum.fun),
        math.exp(float(optimum.x[0])),
        math.exp(float(optimum.x[1])),
        math.exp(float(optimum.x[2])),
    )


def saddle_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    """Chernoff-bound the exact regional endpoint EGF with a free saddle."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    r = message_weight
    regional_constant = math.lgamma(r + 1) - r * math.log(region_length)

    def base_objective(log_z: float, log_v: float, log_t: float) -> float:
        z = math.exp(log_z)
        v = math.exp(log_v)
        t = math.exp(log_t)
        coefficient_factor = region_count * (
            regional_constant - r * log_t
        )
        return (
            coefficient_factor
            + log_weight_mgf(saddle_trace_matrix(t, z, v), n, 0)
            - cutoff * log_z
        )

    def structural_objective(point: np.ndarray) -> float:
        log_z, log_v, log_t = map(float, point)
        return base_objective(log_z, log_v, log_t) - (r - 1) * log_v

    log_s = math.log(prime - 1)

    def field_objective(point: np.ndarray) -> float:
        log_z, log_v, log_t = map(float, point)
        return -log_s + base_objective(log_z, log_v, log_t) - r * log_v

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    v_start = max(1e-12, min(0.8, r / max(1.0, n)))
    t_start = r / region_length
    lower = -max(120.0, log_s + 12.0)
    structural, _, _, _ = _optimize_three_markers(
        objective=structural_objective,
        z_start=z_start,
        v_start=v_start,
        t_start=t_start,
        lower_z=lower,
        lower_v=lower,
    )
    field, _, _, _ = _optimize_three_markers(
        objective=field_objective,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        t_start=t_start,
        lower_z=lower,
        lower_v=-log_s,
    )
    support_count = log_binom(k, r)
    structural += support_count
    field += support_count
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return RegularProjectiveTraceTerm(
        message_weight=r,
        empty_set_weight=math.exp(-r / region_length),
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
    )


def balanced_slot_trace_matrix(
    group_size: int, slot_marker: float, z: float, v: float
) -> np.ndarray:
    """Trace matrix for one balanced output group of ``group_size`` slots."""
    empty, occupied = trace_input_matrices(z, v)
    log_factor = group_size * math.log1p(slot_marker)
    if log_factor > 700.0:
        raise OverflowError("balanced slot marker is outside the stable range")
    return empty + math.expm1(log_factor) * occupied


def balanced_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    """Coefficient-saddle bound for a two-sided regular regional expander."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k in the balanced ensemble")
    group_size = k // region_length
    if group_size < 1:
        raise ValueError("balanced output groups must be nonempty")
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    log_s = math.log(prime - 1)

    def finish(base_objective, *, has_slot_marker: bool) -> RegularProjectiveTraceTerm:
        def structural_objective(point: np.ndarray) -> float:
            if has_slot_marker:
                log_z, log_v, log_x = map(float, point)
                base = base_objective(log_z, log_v, log_x)
            else:
                log_z, log_v = map(float, point)
                base = base_objective(log_z, log_v)
            return base - (r - 1) * log_v

        def field_objective(point: np.ndarray) -> float:
            if has_slot_marker:
                log_z, log_v, log_x = map(float, point)
                base = base_objective(log_z, log_v, log_x)
            else:
                log_z, log_v = map(float, point)
                base = base_objective(log_z, log_v)
            return -log_s + base - r * log_v

        z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
        v_start = max(1e-12, min(0.8, r / max(1.0, n)))
        lower = -max(120.0, log_s + 12.0)
        if has_slot_marker:
            x_start = r / max(1.0, k - r)
            upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
            structural, _, _, _ = _optimize_three_markers(
                objective=structural_objective,
                z_start=z_start,
                v_start=v_start,
                t_start=x_start,
                lower_z=lower,
                lower_v=lower,
                upper_t=upper_x,
            )
            field, _, _, _ = _optimize_three_markers(
                objective=field_objective,
                z_start=z_start,
                v_start=max(math.exp(-log_s), v_start),
                t_start=x_start,
                lower_z=lower,
                lower_v=-log_s,
                upper_t=upper_x,
            )
        else:
            structural, _, _ = _optimize_two_markers(
                objective=structural_objective,
                z_start=z_start,
                v_start=v_start,
                lower_z=lower,
                lower_v=lower,
            )
            field, _, _ = _optimize_two_markers(
                objective=field_objective,
                z_start=z_start,
                v_start=max(math.exp(-log_s), v_start),
                lower_z=lower,
                lower_v=-log_s,
            )
        support_count = log_binom(k, r)
        structural += support_count
        field += support_count
        total = float(np.logaddexp(structural, field))
        inverse_log_two = 1.0 / math.log(2.0)
        empty_probability = (
            0.0
            if k - r < group_size
            else math.exp(
                log_binom(k - r, group_size) - log_binom(k, group_size)
            )
        )
        return RegularProjectiveTraceTerm(
            message_weight=r,
            empty_set_weight=empty_probability,
            structural_log2=structural * inverse_log_two,
            field_log2=field * inverse_log_two,
            total_log2=total * inverse_log_two,
        )

    if r == k:
        def full_base(log_z: float, log_v: float) -> float:
            _, occupied = trace_input_matrices(
                math.exp(log_z), math.exp(log_v)
            )
            return log_weight_mgf(occupied, n, 0) - cutoff * log_z

        return finish(full_base, has_slot_marker=False)

    regional_normalization = region_count * log_binom(k, r)

    def coefficient_base(log_z: float, log_v: float, log_x: float) -> float:
        z = math.exp(log_z)
        v = math.exp(log_v)
        x = math.exp(log_x)
        return (
            -regional_normalization
            - region_count * r * log_x
            + log_weight_mgf(
                balanced_slot_trace_matrix(group_size, x, z, v), n, 0
            )
            - cutoff * log_z
        )

    return finish(coefficient_base, has_slot_marker=True)


def occupancy_size_distribution(region_length: int, draws: int) -> np.ndarray:
    """Return the number of occupied bins after ``draws`` uniform draws."""
    if region_length < 1 or draws < 0:
        raise ValueError("invalid occupancy parameters")
    distribution = np.array([1.0], dtype=np.float64)
    for _ in range(draws):
        following = np.zeros(len(distribution) + 1, dtype=np.float64)
        for occupied, probability in enumerate(distribution):
            following[occupied] += probability * occupied / region_length
            following[occupied + 1] += (
                probability * (region_length - occupied) / region_length
            )
        distribution = following
    return distribution


def balanced_occupancy_size_distribution(
    region_length: int, group_size: int, draws: int
) -> np.ndarray:
    """Occupied groups after drawing slots without replacement."""
    total_slots = region_length * group_size
    if region_length < 1 or group_size < 1 or not 0 <= draws <= total_slots:
        raise ValueError("invalid balanced occupancy parameters")
    distribution = np.array([1.0], dtype=np.float64)
    for used in range(draws):
        following = np.zeros(len(distribution) + 1, dtype=np.float64)
        remaining = total_slots - used
        for occupied, probability in enumerate(distribution):
            following[occupied] += (
                probability * (group_size * occupied - used) / remaining
            )
            following[occupied + 1] += (
                probability * group_size * (region_length - occupied) / remaining
            )
        distribution = following
    return distribution


def trace_input_matrices(z: float, v: float) -> tuple[np.ndarray, np.ndarray]:
    """Return trace transfers for an empty and an occupied coordinate."""
    empty = np.array(((1.0, 0.0), (0.0, z)), dtype=np.float64)
    occupied = np.array(((v, z), (v, z)), dtype=np.float64)
    return empty, occupied


def uniform_occupancy_trace_transfers(
    *, region_length: int, max_occupied: int, z: float, v: float,
) -> list[np.ndarray]:
    """Return exact uniform-occupied-set transfers through one region."""
    if not 0 <= max_occupied <= region_length:
        raise ValueError("max_occupied must lie in [0, region_length]")
    empty, occupied = trace_input_matrices(z, v)
    result: tuple[int, list[np.ndarray]] = (0, [np.eye(2, dtype=np.float64)])
    power: tuple[int, list[np.ndarray]] = (1, [empty, occupied])
    remaining = region_length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers(
                result, power, max_occupied
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers(
                power, power, max_occupied
            )
    return result[1]


def exact_region_trace_matrix(
    *, region_length: int, message_weight: int, z: float, v: float,
) -> np.ndarray:
    """Average the trace transfer over the exact regional occupancy law."""
    occupancy = occupancy_size_distribution(region_length, message_weight)
    slices = uniform_occupancy_trace_transfers(
        region_length=region_length,
        max_occupied=min(region_length, message_weight),
        z=z,
        v=v,
    )
    region = np.zeros((2, 2), dtype=np.float64)
    for occupied, probability in enumerate(occupancy):
        if probability:
            region += probability * slices[occupied]
    return region


def exact_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    """Exact regular-incidence trace bound for one message support size."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count

    def base_objective(log_z: float, log_v: float) -> float:
        region = exact_region_trace_matrix(
            region_length=region_length,
            message_weight=message_weight,
            z=math.exp(log_z),
            v=math.exp(log_v),
        )
        return log_weight_mgf(region, region_count, 0) - cutoff * log_z

    def structural_objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return base_objective(log_z, log_v) - (message_weight - 1) * log_v

    log_s = math.log(prime - 1)

    def field_objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return -log_s + base_objective(log_z, log_v) - message_weight * log_v

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    v_start = max(1e-12, min(0.8, message_weight / max(1.0, n)))
    lower = -max(120.0, log_s + 12.0)
    structural, _, _ = _optimize_two_markers(
        objective=structural_objective,
        z_start=z_start,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    field, _, _ = _optimize_two_markers(
        objective=field_objective,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        lower_z=lower,
        lower_v=-log_s,
    )
    support_count = log_binom(k, message_weight)
    structural += support_count
    field += support_count
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return RegularProjectiveTraceTerm(
        message_weight=message_weight,
        empty_set_weight=math.exp(-message_weight / region_length),
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
    )


def exact_balanced_projective_trace_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> RegularProjectiveTraceTerm:
    """Exact small-support trace bound for the two-sided regular ensemble."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k in the balanced ensemble")
    group_size = k // region_length
    occupancy = balanced_occupancy_size_distribution(
        region_length, group_size, message_weight
    )

    def region_matrix(z: float, v: float) -> np.ndarray:
        slices = uniform_occupancy_trace_transfers(
            region_length=region_length,
            max_occupied=min(region_length, message_weight),
            z=z,
            v=v,
        )
        region = np.zeros((2, 2), dtype=np.float64)
        for occupied, probability in enumerate(occupancy):
            if probability:
                region += probability * slices[occupied]
        return region

    def base_objective(log_z: float, log_v: float) -> float:
        return (
            log_weight_mgf(
                region_matrix(math.exp(log_z), math.exp(log_v)),
                region_count,
                0,
            )
            - cutoff * log_z
        )

    def structural_objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return base_objective(log_z, log_v) - (message_weight - 1) * log_v

    log_s = math.log(prime - 1)

    def field_objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return -log_s + base_objective(log_z, log_v) - message_weight * log_v

    z_start = max(1e-9, min(0.95, cutoff / max(1.0, n)))
    v_start = max(1e-12, min(0.8, message_weight / max(1.0, n)))
    lower = -max(120.0, log_s + 12.0)
    structural, _, _ = _optimize_two_markers(
        objective=structural_objective,
        z_start=z_start,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    field, _, _ = _optimize_two_markers(
        objective=field_objective,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        lower_z=lower,
        lower_v=-log_s,
    )
    support_count = log_binom(k, message_weight)
    structural += support_count
    field += support_count
    total = float(np.logaddexp(structural, field))
    inverse_log_two = 1.0 / math.log(2.0)
    return RegularProjectiveTraceTerm(
        message_weight=message_weight,
        empty_set_weight=(k - message_weight) / k,
        structural_log2=structural * inverse_log_two,
        field_log2=field * inverse_log_two,
        total_log2=total * inverse_log_two,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime", type=int, default=170141183460469231731687303715884105727
    )
    parser.add_argument("--n", type=int, default=2_097_140)
    parser.add_argument("--k", type=int)
    parser.add_argument("--regions", type=int, default=194)
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--gv-fraction", type=float, default=1.0)
    parser.add_argument("--weights", type=int, nargs="+")
    parser.add_argument(
        "--exact-limit", type=int, default=5,
        help="use the exact regional occupancy transfer through this weight",
    )
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
        f"regions={args.regions} region_length={args.n//args.regions} "
        f"cutoff={cutoff} relative_cutoff={cutoff/args.n:.12f} "
        f"gv={delta_gv:.12f}"
    )
    terms = []
    for weight in weights:
        dense = regular_projective_trace_logterm(
            prime=args.prime,
            k=k,
            n=args.n,
            cutoff=cutoff,
            region_count=args.regions,
            message_weight=weight,
        )
        poisson = poisson_projective_trace_logterm(
            prime=args.prime,
            k=k,
            n=args.n,
            cutoff=cutoff,
            region_count=args.regions,
            message_weight=weight,
        )
        saddle = saddle_projective_trace_logterm(
            prime=args.prime,
            k=k,
            n=args.n,
            cutoff=cutoff,
            region_count=args.regions,
            message_weight=weight,
        )
        exact = (
            exact_projective_trace_logterm(
                prime=args.prime,
                k=k,
                n=args.n,
                cutoff=cutoff,
                region_count=args.regions,
                message_weight=weight,
            )
            if weight <= args.exact_limit
            else None
        )
        candidates = [dense, poisson, saddle]
        if exact is not None:
            candidates.append(exact)
        term = min(candidates, key=lambda item: item.total_log2)
        terms.append(term)
        if not args.summary_only:
            exact_text = (
                f"exact={exact.total_log2:.6f}\t" if exact is not None else ""
            )
            print(
                f"r={weight}\tbeta={term.empty_set_weight:.8g}\t"
                f"dense={dense.total_log2:.6f}\t"
                f"poisson={poisson.total_log2:.6f}\t"
                f"saddle={saddle.total_log2:.6f}\t"
                f"{exact_text}"
                f"best={term.total_log2:.6f}"
            )
    worst = max(terms, key=lambda item: item.total_log2)
    print(
        f"sampled combined worst: r={worst.message_weight}, "
        f"log2_bound={worst.total_log2:.6f}"
    )


if __name__ == "__main__":
    main()
