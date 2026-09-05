#!/usr/bin/env python3
"""Floating-point diagnostics for two-sided regular prime-field EC.

The convolution is nonwrapping and time varying.  At time ``t`` it samples
an independent feedback vector uniformly from ``F_p^memory``.  The certified
diagnostic counts every fresh zero equation before applying the projective
cap, then uses a positive coefficient saddle for the balanced regional
expander.  Probability-transfer helpers remain for uncapped experiments only.
This script selects candidate parameters; it is not an interval certificate.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.fft import fft2, ifft2, next_fast_len
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

from expander_bounds import log_binom, log_weight_mgf, normalized
from prime_field_ea_diagnostic import gv_distance
from prime_field_regular_ea_diagnostic import log_qary_hamming_ball_geometric
from prime_field_regular_ea_trace_diagnostic import (
    _optimize_three_markers,
    _optimize_two_markers,
    balanced_occupancy_size_distribution,
    diagnostic_weights,
)
from regular_ec_diagnostic import combine_uniform_slice_transfers


@dataclass(frozen=True)
class BiregularECTerm:
    message_weight: int
    structural_log2: float
    field_log2: float
    total_log2: float


@dataclass(frozen=True)
class BiregularECPrefixBound:
    support_start: int
    support_limit: int
    output_marker: float
    equation_marker: float
    total_log2: float


@dataclass(frozen=True)
class ConstraintECBandBound:
    """A support-band bound and the two projective-cap saddles."""

    support_start: int
    support_limit: int
    structural_log2: float
    field_log2: float
    total_log2: float
    structural_markers: tuple[float, ...]
    field_markers: tuple[float, ...]


@dataclass(frozen=True)
class ShellStackECTerm:
    """Coarse projective cap after summing one fixed-line convolution tail."""

    message_weight: int
    line_tail_log2: float
    projective_union_log2: float
    total_log2: float
    output_marker: float
    cap_saturated: bool


def optimize_output_marker_modes(
    *, objective, z_starts: tuple[float, ...], lower_z: float
) -> tuple[float, float]:
    candidates = [
        minimize(
            lambda point: objective(float(point[0])),
            np.array([math.log(z_start)], dtype=np.float64),
            method="L-BFGS-B",
            bounds=((lower_z, -1e-11),),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
        )
        for z_start in z_starts
    ]
    optimum = min(candidates, key=lambda result: float(result.fun))
    return float(optimum.fun), math.exp(float(optimum.x[0]))


def optimize_output_slot_marker_modes(
    *, objective, z_starts: tuple[float, ...], x_start: float,
    lower_z: float, upper_x: float,
) -> tuple[float, float, float]:
    lower_x = -30.0
    x_start = max(math.exp(lower_x), min(math.exp(upper_x), x_start))
    candidates = []
    for z_start in z_starts:
        for scale in (0.8, 1.0, 1.2):
            start_x = max(
                math.exp(lower_x), min(math.exp(upper_x), scale * x_start)
            )
            candidates.append(minimize(
                objective,
                np.log(np.array([z_start, start_x], dtype=np.float64)),
                method="L-BFGS-B",
                bounds=((lower_z, -1e-11), (lower_x, upper_x)),
                options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 700},
            ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    return (
        float(optimum.fun),
        math.exp(float(optimum.x[0])),
        math.exp(float(optimum.x[1])),
    )


def output_marker_starts(
    *, prime: int, n: int, cutoff: int, memory: int, lower_log: float
) -> tuple[float, ...]:
    """Return starts for the active-output and rare-termination saddle modes."""
    active = cutoff / ((prime - 1) * (n - cutoff))
    termination = math.exp(-memory * math.log(prime) / max(1, cutoff))
    floor = math.exp(lower_log)
    values = {
        max(floor, min(0.999999999, active)),
        max(floor, min(0.999999999, termination)),
        max(floor, min(0.999999999, math.sqrt(active * termination))),
        0.5,
    }
    return tuple(sorted(values))


def optimize_two_marker_modes(
    *, objective, z_starts: tuple[float, ...], v_start: float,
    lower_z: float, lower_v: float,
) -> tuple[float, float, float]:
    return min(
        (
            _optimize_two_markers(
                objective=objective,
                z_start=z_start,
                v_start=v_start,
                lower_z=lower_z,
                lower_v=lower_v,
            )
            for z_start in z_starts
        ),
        key=lambda result: result[0],
    )


def optimize_substochastic_two_markers(
    *, objective, z_starts: tuple[float, ...], v_start: float,
    lower_log: float, v_floor: float,
) -> tuple[float, float, float]:
    """Optimize positive markers with ``z+v<=1`` by construction."""
    if not 0.0 <= v_floor < 1.0:
        raise ValueError("v_floor must lie in [0,1)")
    scale = 1.0 - v_floor
    candidates = []
    for z_candidate in z_starts:
        z = max(
            math.exp(lower_log) * scale,
            min(scale * (1.0 - 1e-11), z_candidate),
        )
        residual = max(math.exp(lower_log), 1.0 - v_floor - z)
        w = max(
            math.exp(lower_log),
            min(1.0, (v_start - v_floor) / residual),
        )

        def transformed(point: np.ndarray) -> float:
            log_z_scale, log_w = map(float, point)
            local_z = scale * math.exp(log_z_scale)
            local_v = v_floor + (1.0 - v_floor - local_z) * math.exp(log_w)
            if local_v <= 0.0:
                return 1e300
            return objective(np.array(
                [math.log(local_z), math.log(local_v)], dtype=np.float64
            ))

        candidates.append(minimize(
            transformed,
            np.array([math.log(z / scale), math.log(w)], dtype=np.float64),
            method="L-BFGS-B",
            bounds=((lower_log, -1e-11), (lower_log, 0.0)),
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
        ))
    optimum = min(candidates, key=lambda result: float(result.fun))
    log_z_scale, log_w = map(float, optimum.x)
    z = scale * math.exp(log_z_scale)
    v = v_floor + (1.0 - v_floor - z) * math.exp(log_w)
    return float(optimum.fun), z, v


def optimize_three_marker_modes(
    *, objective, z_starts: tuple[float, ...], v_start: float,
    x_start: float, lower_z: float, lower_v: float, upper_x: float,
) -> tuple[float, float, float, float]:
    return min(
        (
            _optimize_three_markers(
                objective=objective,
                z_start=z_start,
                v_start=v_start,
                t_start=x_start,
                lower_z=lower_z,
                lower_v=lower_v,
                upper_t=upper_x,
            )
            for z_start in z_starts
        ),
        key=lambda result: result[0],
    )


def _log_expm1_positive(value: float) -> float:
    """Return log(exp(value)-1) without overflow or cancellation."""
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def balanced_point_mass_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> float:
    """Bound one projective support family using convolution invertibility.

    Conditional on the balanced incidence pattern, independent edge labels
    give maximum mass at most ``(p-1)^(-u)`` when ``u`` output coordinates
    are occupied.  A coefficient saddle averages this bound over incidence.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    log_s = math.log(prime - 1)

    if r == k:
        regional_log_mass = -region_length * log_s
    else:
        def objective(log_x: float) -> float:
            log_one_plus_x = float(np.logaddexp(0.0, log_x))
            log_occupied_factor = _log_expm1_positive(
                group_size * log_one_plus_x
            )
            log_polynomial = float(
                np.logaddexp(0.0, log_occupied_factor - log_s)
            )
            return region_length * log_polynomial - r * log_x

        optimum = minimize_scalar(
            objective,
            bounds=(-50.0, log_s + 50.0),
            method="bounded",
            options={"xatol": 1e-11, "maxiter": 300},
        )
        regional_log_mass = float(optimum.fun) - log_binom(k, r)

    return (
        log_binom(k, r)
        + (r - 1) * log_s
        + log_qary_hamming_ball_geometric(prime, n, cutoff)
        + region_count * regional_log_mass
    ) / math.log(2.0)


def nonwrapping_trace_matrices(
    *, prime: int, memory: int, output_marker: float, equation_marker: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return probability transfers for uncapped fixed-message experiments."""
    if prime < 2 or memory < 1:
        raise ValueError("invalid prime-field convolution parameters")
    if output_marker <= 0.0 or equation_marker <= 0.0:
        raise ValueError("markers must be positive")

    size = memory + 1
    empty = np.zeros((size, size), dtype=np.float64)
    occupied = np.zeros((size, size), dtype=np.float64)
    zero_probability = 1.0 / prime
    nonzero_probability = (prime - 1.0) / prime

    # In a nonzero convolution state, a fresh uniform feedback vector makes
    # the output uniform in F_p.  The expander input does not affect the law.
    for state in range(memory):
        empty[state, 0] = occupied[state, 0] = (
            nonzero_probability * output_marker
        )
        empty[state, state + 1] = occupied[state, state + 1] = zero_probability

    # In the zero convolution state, the output equals the expander input.
    empty[memory, memory] = 1.0
    occupied[memory, 0] = output_marker
    occupied[memory, memory] = equation_marker
    return empty, occupied


def nonwrapping_constraint_trace_matrices(
    *, memory: int, output_marker: float, equation_marker: float
) -> tuple[np.ndarray, np.ndarray]:
    """Trace output weight and every fresh zero constraint.

    A zero produced from a nonzero convolution state is one affine equation in
    the fresh feedback vector.  A zero expanded symbol in the zero state is
    one equation in fresh edge labels.  Nonzero outcomes are assigned weight
    one, which is an upper bound on their probability.
    """
    if memory < 1 or output_marker <= 0.0 or equation_marker <= 0.0:
        raise ValueError("invalid constraint-trace parameters")
    size = memory + 1
    empty = np.zeros((size, size), dtype=np.float64)
    occupied = np.zeros((size, size), dtype=np.float64)
    for state in range(memory):
        empty[state, 0] = occupied[state, 0] = output_marker
        empty[state, state + 1] = occupied[state, state + 1] = equation_marker
    empty[memory, memory] = 1.0
    occupied[memory, 0] = output_marker
    occupied[memory, memory] = equation_marker
    return empty, occupied


def nonwrapping_singleton_constraint_trace_matrices(
    *, memory: int, output_marker: float, equation_marker: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Trace empty, singleton, and collision expander coordinates.

    A singleton coordinate is a nonzero message symbol times a nonzero edge
    label, so it cannot vanish while the convolution state is zero.  Lumping
    it with collision coordinates loses substantial information at sparse
    message support.
    """
    empty, collision = nonwrapping_constraint_trace_matrices(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    singleton = collision.copy()
    singleton[memory, memory] = 0.0
    return empty, singleton, collision


def balanced_constraint_slot_trace_matrix(
    *, memory: int, group_size: int, slot_marker: float,
    output_marker: float, equation_marker: float,
) -> np.ndarray:
    empty, occupied = nonwrapping_constraint_trace_matrices(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    log_factor = group_size * math.log1p(slot_marker)
    if log_factor > 700.0:
        raise OverflowError("slot marker is outside the stable range")
    return empty + math.expm1(log_factor) * occupied


def balanced_singleton_constraint_slot_trace_matrix(
    *, memory: int, group_size: int, slot_marker: float,
    output_marker: float, equation_marker: float,
) -> np.ndarray:
    """Positive saddle matrix retaining singleton occupied groups."""
    empty, singleton, collision = (
        nonwrapping_singleton_constraint_trace_matrices(
            memory=memory,
            output_marker=output_marker,
            equation_marker=equation_marker,
        )
    )
    log_factor = group_size * math.log1p(slot_marker)
    if log_factor > 700.0:
        raise OverflowError("slot marker is outside the stable range")
    singleton_factor = group_size * slot_marker
    collision_factor = math.expm1(log_factor) - singleton_factor
    return empty + singleton_factor * singleton + collision_factor * collision


def balanced_slot_trace_matrix(
    *, prime: int, memory: int, group_size: int, slot_marker: float,
    output_marker: float, equation_marker: float,
) -> np.ndarray:
    """Return Q0 + ((1+x)^q-1) Q1 for one balanced output group."""
    empty, occupied = nonwrapping_trace_matrices(
        prime=prime,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    log_factor = group_size * math.log1p(slot_marker)
    if log_factor > 700.0:
        raise OverflowError("slot marker is outside the stable range")
    return empty + math.expm1(log_factor) * occupied


def uniform_occupancy_trace_transfers(
    *, length: int, max_occupied: int, prime: int, memory: int,
    output_marker: float, equation_marker: float,
) -> list[np.ndarray]:
    """Return normalized transfers for uniform occupied subsets."""
    empty, occupied = nonwrapping_trace_matrices(
        prime=prime,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    identity = np.eye(memory + 1, dtype=np.float64)
    result: tuple[int, list[np.ndarray]] = (0, [identity])
    power: tuple[int, list[np.ndarray]] = (1, [empty, occupied])
    remaining = length
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


def _float_uniform_subset_transfers(
    *, empty: np.ndarray, occupied: np.ndarray, length: int,
    max_occupied: int,
) -> list[np.ndarray]:
    """Float-only uniform-subset convolution without enormous binomials."""
    size = empty.shape[0]
    identity = np.eye(size, dtype=np.float64)
    result: tuple[int, list[np.ndarray]] = (0, [identity])
    power: tuple[int, list[np.ndarray]] = (1, [empty, occupied])

    def combine(
        left: tuple[int, list[np.ndarray]],
        right: tuple[int, list[np.ndarray]],
    ) -> tuple[int, list[np.ndarray]]:
        left_length, left_slices = left
        right_length, right_slices = right
        total_length = left_length + right_length
        combined: list[np.ndarray] = []
        for weight in range(min(max_occupied, total_length) + 1):
            lo = max(0, weight - right_length)
            hi = min(weight, left_length, len(left_slices) - 1)
            left_weights = np.arange(lo, hi + 1, dtype=np.int64)
            right_weights = weight - left_weights
            valid = right_weights < len(right_slices)
            left_weights = left_weights[valid]
            right_weights = right_weights[valid]
            matrix = np.zeros((size, size), dtype=np.float64)
            if not len(left_weights):
                combined.append(matrix)
                continue

            first = int(left_weights[0])
            log_probability = (
                math.lgamma(left_length + 1)
                - math.lgamma(first + 1)
                - math.lgamma(left_length - first + 1)
                + math.lgamma(right_length + 1)
                - math.lgamma(weight - first + 1)
                - math.lgamma(right_length - weight + first + 1)
                - math.lgamma(total_length + 1)
                + math.lgamma(weight + 1)
                + math.lgamma(total_length - weight + 1)
            )
            probabilities = np.empty(len(left_weights), dtype=np.float64)
            probabilities[0] = math.exp(log_probability)
            for index in range(1, len(left_weights)):
                previous = int(left_weights[index - 1])
                probabilities[index] = probabilities[index - 1] * (
                    (left_length - previous)
                    * (weight - previous)
                    / ((previous + 1)
                       * (right_length - weight + previous + 1))
                )
            products = np.matmul(
                np.stack([left_slices[int(value)] for value in left_weights]),
                np.stack([right_slices[int(value)] for value in right_weights]),
            )
            matrix = np.tensordot(probabilities, products, axes=(0, 0))
            combined.append(matrix)
        return total_length, combined

    remaining = length
    while remaining:
        if remaining & 1:
            result = combine(result, power)
        remaining >>= 1
        if remaining:
            power = combine(power, power)
    return result[1]


def _float_uniform_group_power(
    *, group_size: int, group_slices: list[np.ndarray], group_count: int,
    max_weight: int,
) -> list[np.ndarray]:
    """Normalized coefficients of a matrix group polynomial power.

    ``group_slices[t]`` is the transfer conditional on selecting exactly
    ``t`` of one group's ``group_size`` support slots.  Returned slice ``r``
    is normalized by ``binom(group_size * group_count, r)``.
    """
    if len(group_slices) != group_size + 1:
        raise ValueError("one conditional slice is required per group weight")
    size = group_slices[0].shape[0]
    identity = np.eye(size, dtype=np.float64)

    def combine(
        left: tuple[int, list[np.ndarray]],
        right: tuple[int, list[np.ndarray]],
    ) -> tuple[int, list[np.ndarray]]:
        left_slots, left_slices = left
        right_slots, right_slices = right
        total_slots = left_slots + right_slots
        combined: list[np.ndarray] = []
        for weight in range(min(max_weight, total_slots) + 1):
            lo = max(0, weight - right_slots)
            hi = min(weight, left_slots, len(left_slices) - 1)
            left_weights = np.arange(lo, hi + 1, dtype=np.int64)
            right_weights = weight - left_weights
            valid = right_weights < len(right_slices)
            left_weights = left_weights[valid]
            right_weights = right_weights[valid]
            if not len(left_weights):
                combined.append(np.zeros((size, size), dtype=np.float64))
                continue
            first = int(left_weights[0])
            log_probability = (
                math.lgamma(left_slots + 1)
                - math.lgamma(first + 1)
                - math.lgamma(left_slots - first + 1)
                + math.lgamma(right_slots + 1)
                - math.lgamma(weight - first + 1)
                - math.lgamma(right_slots - weight + first + 1)
                - math.lgamma(total_slots + 1)
                + math.lgamma(weight + 1)
                + math.lgamma(total_slots - weight + 1)
            )
            probabilities = np.empty(len(left_weights), dtype=np.float64)
            probabilities[0] = math.exp(log_probability)
            for index in range(1, len(left_weights)):
                previous = int(left_weights[index - 1])
                probabilities[index] = probabilities[index - 1] * (
                    (left_slots - previous)
                    * (weight - previous)
                    / ((previous + 1)
                       * (right_slots - weight + previous + 1))
                )
            products = np.matmul(
                np.stack([left_slices[int(value)] for value in left_weights]),
                np.stack([right_slices[int(value)] for value in right_weights]),
            )
            combined.append(np.tensordot(probabilities, products, axes=(0, 0)))
        return total_slots, combined

    result: tuple[int, list[np.ndarray]] = (0, [identity])
    power = (group_size, group_slices)
    remaining = group_count
    while remaining:
        if remaining & 1:
            result = combine(result, power)
        remaining >>= 1
        if remaining:
            power = combine(power, power)
    return result[1]


def singleton_constraint_uniform_support_trace_transfers(
    *, region_length: int, group_size: int, max_weight: int, memory: int,
    output_marker: float, equation_marker: float,
) -> list[np.ndarray]:
    """Regional support coefficients distinguishing singleton groups."""
    empty, singleton, collision = (
        nonwrapping_singleton_constraint_trace_matrices(
            memory=memory,
            output_marker=output_marker,
            equation_marker=equation_marker,
        )
    )
    group_slices = [empty, singleton]
    group_slices.extend(collision for _ in range(2, group_size + 1))
    return _float_uniform_group_power(
        group_size=group_size,
        group_slices=group_slices,
        group_count=region_length,
        max_weight=max_weight,
    )


def singleton_constraint_uniform_support_trace_transfers_scaled(
    *, region_length: int, group_size: int, max_weight: int, memory: int,
    output_marker: float, equation_marker: float,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Underflow-safe singleton regional slices and natural-log scales."""
    empty, singleton, collision = (
        nonwrapping_singleton_constraint_trace_matrices(
            memory=memory,
            output_marker=output_marker,
            equation_marker=equation_marker,
        )
    )
    raw_group_slices = [empty, singleton]
    raw_group_slices.extend(collision for _ in range(2, group_size + 1))
    group_slices: list[np.ndarray] = []
    group_scales: list[float] = []
    for matrix in raw_group_slices:
        scaled, log_scale = normalized(matrix)
        group_slices.append(scaled)
        group_scales.append(log_scale)
    size = memory + 1

    def combine(
        left: tuple[int, list[np.ndarray], np.ndarray],
        right: tuple[int, list[np.ndarray], np.ndarray],
    ) -> tuple[int, list[np.ndarray], np.ndarray]:
        left_slots, left_slices, left_scales = left
        right_slots, right_slices, right_scales = right
        total_slots = left_slots + right_slots
        combined: list[np.ndarray] = []
        combined_scales: list[float] = []
        for weight in range(min(max_weight, total_slots) + 1):
            lo = max(0, weight - right_slots)
            hi = min(weight, left_slots, len(left_slices) - 1)
            left_weights = np.arange(lo, hi + 1, dtype=np.int64)
            right_weights = weight - left_weights
            valid = right_weights < len(right_slices)
            left_weights = left_weights[valid]
            right_weights = right_weights[valid]
            if not len(left_weights):
                combined.append(np.zeros((size, size), dtype=np.float64))
                combined_scales.append(-math.inf)
                continue
            log_probabilities = (
                gammaln(left_slots + 1)
                - gammaln(left_weights + 1)
                - gammaln(left_slots - left_weights + 1)
                + gammaln(right_slots + 1)
                - gammaln(right_weights + 1)
                - gammaln(right_slots - right_weights + 1)
                - gammaln(total_slots + 1)
                + gammaln(weight + 1)
                + gammaln(total_slots - weight + 1)
            )
            term_scales = (
                log_probabilities
                + left_scales[left_weights]
                + right_scales[right_weights]
            )
            finite = np.isfinite(term_scales)
            left_weights = left_weights[finite]
            right_weights = right_weights[finite]
            term_scales = term_scales[finite]
            if not len(term_scales):
                combined.append(np.zeros((size, size), dtype=np.float64))
                combined_scales.append(-math.inf)
                continue
            maximum = float(np.max(term_scales))
            products = np.matmul(
                np.stack([left_slices[int(value)] for value in left_weights]),
                np.stack([right_slices[int(value)] for value in right_weights]),
            )
            matrix = np.tensordot(
                np.exp(term_scales - maximum), products, axes=(0, 0)
            )
            matrix, local_scale = normalized(matrix)
            combined.append(matrix)
            combined_scales.append(maximum + local_scale)
        return total_slots, combined, np.asarray(combined_scales)

    result = (0, [np.eye(size, dtype=np.float64)], np.asarray([0.0]))
    power = (group_size, group_slices, np.asarray(group_scales))
    remaining = region_length
    while remaining:
        if remaining & 1:
            result = combine(result, power)
        remaining >>= 1
        if remaining:
            power = combine(power, power)
    return result[1], result[2]


def _normalized_nonnegative_tensor(
    tensor: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Normalize a nonnegative coefficient tensor by its largest entry."""
    scale = float(np.max(tensor))
    if scale == 0.0:
        return tensor, -math.inf
    return tensor / scale, math.log(scale)


def _truncated_bivariate_matrix_product(
    left: np.ndarray, right: np.ndarray, *, max_x_degree: int,
    max_equation_degree: int,
) -> np.ndarray:
    """Multiply matrix polynomials in ``x`` and ``v``, then truncate.

    The first two axes are the polynomial degrees and the final two axes are
    matrix indices.  FFT convolution keeps the sparse-support diagnostic
    practical at convolution memories near ``log2(n)``.
    """
    if left.ndim != 4 or right.ndim != 4:
        raise ValueError("matrix-polynomial tensors must have four axes")
    if left.shape[3] != right.shape[2]:
        raise ValueError("matrix-polynomial dimensions do not conform")
    full_x = left.shape[0] + right.shape[0] - 1
    full_a = left.shape[1] + right.shape[1] - 1
    fft_x = next_fast_len(full_x)
    fft_a = next_fast_len(full_a)
    left_fft = fft2(left, s=(fft_x, fft_a), axes=(0, 1))
    right_fft = fft2(right, s=(fft_x, fft_a), axes=(0, 1))
    product_fft = np.matmul(left_fft, right_fft)
    product = ifft2(product_fft, axes=(0, 1)).real
    product = product[
        :min(max_x_degree + 1, full_x),
        :min(max_equation_degree + 1, full_a),
    ]
    # Exact coefficients are nonnegative.  The inverse FFT can leave tiny
    # negative roundoff in coefficients that should be zero.
    largest = float(np.max(np.abs(product)))
    if largest and float(np.min(product)) < -2e-10 * largest:
        raise ArithmeticError("FFT polynomial convolution lost positivity")
    return np.maximum(product, 0.0)


def _truncated_bivariate_matrix_power_scaled(
    base: np.ndarray, exponent: int, *, max_x_degree: int,
    max_equation_degree: int,
) -> tuple[np.ndarray, float]:
    """Return a scaled truncated power of a nonnegative matrix polynomial."""
    if exponent < 0:
        raise ValueError("matrix-polynomial exponent must be nonnegative")
    size = base.shape[2]
    result = np.zeros((1, 1, size, size), dtype=np.float64)
    result[0, 0] = np.eye(size, dtype=np.float64)
    result_log_scale = 0.0
    power, power_log_scale = _normalized_nonnegative_tensor(base)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _truncated_bivariate_matrix_product(
                result,
                power,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
            result, local_scale = _normalized_nonnegative_tensor(result)
            result_log_scale += power_log_scale + local_scale
        remaining >>= 1
        if remaining:
            power = _truncated_bivariate_matrix_product(
                power,
                power,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
            power, local_scale = _normalized_nonnegative_tensor(power)
            power_log_scale = 2.0 * power_log_scale + local_scale
    return result, result_log_scale


def _normalize_matrix_polynomial_coefficients(
    polynomial: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize every coefficient matrix independently."""
    normalized_polynomial = np.zeros_like(polynomial)
    scales = np.full(polynomial.shape[:2], -math.inf, dtype=np.float64)
    for x_degree in range(polynomial.shape[0]):
        for equation_degree in range(polynomial.shape[1]):
            matrix, log_scale = normalized(
                polynomial[x_degree, equation_degree]
            )
            normalized_polynomial[x_degree, equation_degree] = matrix
            scales[x_degree, equation_degree] = log_scale
    return normalized_polynomial, scales


def _scaled_bivariate_matrix_product(
    left: np.ndarray, left_scales: np.ndarray,
    right: np.ndarray, right_scales: np.ndarray, *, max_x_degree: int,
    max_equation_degree: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Multiply independently scaled matrix-polynomial coefficients."""
    size = left.shape[2]
    output_x = min(max_x_degree + 1, left.shape[0] + right.shape[0] - 1)
    output_a = min(
        max_equation_degree + 1, left.shape[1] + right.shape[1] - 1
    )
    result = np.zeros((output_x, output_a, size, size), dtype=np.float64)
    result_scales = np.full((output_x, output_a), -math.inf, dtype=np.float64)
    for x_degree in range(output_x):
        x_left_lo = max(0, x_degree - right.shape[0] + 1)
        x_left_hi = min(x_degree, left.shape[0] - 1)
        for equation_degree in range(output_a):
            a_left_lo = max(0, equation_degree - right.shape[1] + 1)
            a_left_hi = min(equation_degree, left.shape[1] - 1)
            left_terms: list[np.ndarray] = []
            right_terms: list[np.ndarray] = []
            term_scales: list[float] = []
            for x_left in range(x_left_lo, x_left_hi + 1):
                x_right = x_degree - x_left
                for a_left in range(a_left_lo, a_left_hi + 1):
                    a_right = equation_degree - a_left
                    log_scale = (
                        float(left_scales[x_left, a_left])
                        + float(right_scales[x_right, a_right])
                    )
                    if math.isfinite(log_scale):
                        left_terms.append(left[x_left, a_left])
                        right_terms.append(right[x_right, a_right])
                        term_scales.append(log_scale)
            if not term_scales:
                continue
            maximum = max(term_scales)
            products = np.matmul(
                np.stack(left_terms), np.stack(right_terms)
            )
            matrix = np.tensordot(
                np.exp(np.asarray(term_scales) - maximum),
                products,
                axes=(0, 0),
            )
            matrix, local_scale = normalized(matrix)
            result[x_degree, equation_degree] = matrix
            result_scales[x_degree, equation_degree] = maximum + local_scale
    return result, result_scales


def _scaled_bivariate_matrix_power(
    base: np.ndarray, exponent: int, *, max_x_degree: int,
    max_equation_degree: int, base_scales: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Raise a matrix polynomial with a log scale for each coefficient."""
    if exponent < 0:
        raise ValueError("matrix-polynomial exponent must be nonnegative")
    size = base.shape[2]
    result = np.zeros((1, 1, size, size), dtype=np.float64)
    result[0, 0] = np.eye(size, dtype=np.float64)
    result_scales = np.asarray([[0.0]], dtype=np.float64)
    if base_scales is None:
        power, power_scales = _normalize_matrix_polynomial_coefficients(base)
    else:
        if base_scales.shape != base.shape[:2]:
            raise ValueError("base scales must match polynomial degrees")
        power = base.copy()
        power_scales = base_scales.copy()
    remaining = exponent
    while remaining:
        if remaining & 1:
            result, result_scales = _scaled_bivariate_matrix_product(
                result,
                result_scales,
                power,
                power_scales,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
        remaining >>= 1
        if remaining:
            power, power_scales = _scaled_bivariate_matrix_product(
                power,
                power_scales,
                power,
                power_scales,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
    return result, result_scales


def singleton_constraint_low_equation_exact_point_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int, output_marker: float,
) -> float:
    """Directly sum the structural branch with fewer than ``r`` equations.

    This removes the equation-marker Chernoff slack from the sparse-support
    branch.  The output-weight event is still bounded by its positive marker,
    so the result is a valid floating-point diagnostic rather than an interval
    certificate.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    if not 0.0 < output_marker <= 1.0:
        raise ValueError("output_marker must lie in (0,1]")
    max_a = r - 1
    size = memory + 1

    def trace_polynomial(kind: str) -> np.ndarray:
        matrix = np.zeros((2, size, size), dtype=np.float64)
        for state in range(memory):
            matrix[0, state, 0] = output_marker
            matrix[1, state, state + 1] = 1.0
        if kind == "empty":
            matrix[0, memory, memory] = 1.0
        else:
            matrix[0, memory, 0] = output_marker
            if kind == "collision":
                matrix[1, memory, memory] = 1.0
        return matrix

    empty = trace_polynomial("empty")
    singleton = trace_polynomial("singleton")
    collision = trace_polynomial("collision")
    group = np.zeros(
        (min(group_size, r) + 1, 2, size, size), dtype=np.float64
    )
    group[0] = empty
    if r >= 1:
        group[1] = group_size * singleton
    for weight in range(2, min(group_size, r) + 1):
        group[weight] = math.comb(group_size, weight) * collision

    region, region_scales = _scaled_bivariate_matrix_power(
        group,
        region_length,
        max_x_degree=r,
        max_equation_degree=max_a,
    )
    regional = region[r:r + 1]
    regional_scales = region_scales[r:r + 1] - log_binom(k, r)
    whole, whole_scales = _scaled_bivariate_matrix_power(
        regional,
        region_count,
        max_x_degree=0,
        max_equation_degree=max_a,
        base_scales=regional_scales,
    )
    endpoint_logs = []
    for equation_degree in range(max_a + 1):
        endpoint_mass = float(np.sum(whole[0, equation_degree, memory, :]))
        if endpoint_mass:
            endpoint_logs.append(
                float(whole_scales[0, equation_degree])
                + math.log(endpoint_mass)
            )
    if not endpoint_logs:
        return -math.inf
    trace_log = float(np.logaddexp.reduce(np.asarray(endpoint_logs)))
    return (
        log_binom(k, r) + trace_log - cutoff * math.log(output_marker)
    ) / math.log(2.0)


def _log_bivariate_matrix_product(
    left: np.ndarray, right: np.ndarray, *, max_x_degree: int,
    max_equation_degree: int,
) -> np.ndarray:
    """Matrix-polynomial product in the nonnegative log semiring."""
    size = left.shape[2]
    output_x = min(max_x_degree + 1, left.shape[0] + right.shape[0] - 1)
    output_a = min(
        max_equation_degree + 1, left.shape[1] + right.shape[1] - 1
    )
    result = np.full((output_x, output_a, size, size), -math.inf)
    for x_degree in range(output_x):
        x_left_lo = max(0, x_degree - right.shape[0] + 1)
        x_left_hi = min(x_degree, left.shape[0] - 1)
        for equation_degree in range(output_a):
            a_left_lo = max(0, equation_degree - right.shape[1] + 1)
            a_left_hi = min(equation_degree, left.shape[1] - 1)
            left_terms = []
            right_terms = []
            for x_left in range(x_left_lo, x_left_hi + 1):
                x_right = x_degree - x_left
                for a_left in range(a_left_lo, a_left_hi + 1):
                    a_right = equation_degree - a_left
                    left_terms.append(left[x_left, a_left])
                    right_terms.append(right[x_right, a_right])
            if left_terms:
                left_stack = np.stack(left_terms)
                right_stack = np.stack(right_terms)
                # Axes are (coefficient split, row, inner, column).
                products = (
                    left_stack[:, :, :, np.newaxis]
                    + right_stack[:, np.newaxis, :, :]
                )
                result[x_degree, equation_degree] = logsumexp(
                    products, axis=(0, 2)
                )
    return result


def _log_bivariate_matrix_power(
    base: np.ndarray, exponent: int, *, max_x_degree: int,
    max_equation_degree: int,
) -> np.ndarray:
    """Truncated matrix-polynomial power in the log semiring."""
    if exponent < 0:
        raise ValueError("matrix-polynomial exponent must be nonnegative")
    size = base.shape[2]
    result = np.full((1, 1, size, size), -math.inf)
    diagonal = np.arange(size)
    result[0, 0, diagonal, diagonal] = 0.0
    power = base
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _log_bivariate_matrix_product(
                result,
                power,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
        remaining >>= 1
        if remaining:
            power = _log_bivariate_matrix_product(
                power,
                power,
                max_x_degree=max_x_degree,
                max_equation_degree=max_equation_degree,
            )
    return result


def singleton_constraint_low_equation_log_point_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int, output_marker: float,
) -> float:
    """Underflow-safe direct structural bound with fewer than ``r`` equations."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    if not 0.0 < output_marker <= 1.0:
        raise ValueError("output_marker must lie in (0,1]")
    max_a = r - 1
    size = memory + 1
    log_z = math.log(output_marker)

    def trace_log_polynomial(kind: str) -> np.ndarray:
        matrix = np.full((2, size, size), -math.inf)
        for state in range(memory):
            matrix[0, state, 0] = log_z
            matrix[1, state, state + 1] = 0.0
        if kind == "empty":
            matrix[0, memory, memory] = 0.0
        else:
            matrix[0, memory, 0] = log_z
            if kind == "collision":
                matrix[1, memory, memory] = 0.0
        return matrix

    empty = trace_log_polynomial("empty")
    singleton = trace_log_polynomial("singleton")
    collision = trace_log_polynomial("collision")
    group = np.full(
        (min(group_size, r) + 1, 2, size, size), -math.inf
    )
    group[0] = empty
    if r >= 1:
        group[1] = singleton + math.log(group_size)
    for weight in range(2, min(group_size, r) + 1):
        group[weight] = collision + math.log(math.comb(group_size, weight))

    region = _log_bivariate_matrix_power(
        group,
        region_length,
        max_x_degree=r,
        max_equation_degree=max_a,
    )
    regional = region[r:r + 1] - log_binom(k, r)
    whole = _log_bivariate_matrix_power(
        regional,
        region_count,
        max_x_degree=0,
        max_equation_degree=max_a,
    )
    trace_log = float(logsumexp(whole[0, :, memory, :]))
    return (
        log_binom(k, r) + trace_log - cutoff * log_z
    ) / math.log(2.0)


def singleton_constraint_exact_point_at_markers_scaled(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int, output_marker: float,
    equation_marker: float,
) -> float:
    """Underflow-safe structural-branch log2 bound for one support."""
    region_length = n // region_count
    group_size = k // region_length
    r = message_weight
    slices, scales = singleton_constraint_uniform_support_trace_transfers_scaled(
        region_length=region_length,
        group_size=group_size,
        max_weight=r,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    value = (
        log_binom(k, r)
        + log_weight_mgf(slices[r], region_count, memory)
        + region_count * float(scales[r])
        - cutoff * math.log(output_marker)
        - (r - 1) * math.log(equation_marker)
    )
    return value / math.log(2.0)


def constraint_uniform_occupancy_trace_transfers(
    *, length: int, max_occupied: int, memory: int,
    output_marker: float, equation_marker: float,
) -> list[np.ndarray]:
    """Return uniform-subset transfers for the fresh-constraint trace."""
    empty, occupied = nonwrapping_constraint_trace_matrices(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    return _float_uniform_subset_transfers(
        empty=empty,
        occupied=occupied,
        length=length,
        max_occupied=max_occupied,
    )


def constraint_exact_balanced_ec_band_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> ConstraintECBandBound:
    """Bound an exact support band by counting all fresh zero equations."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")

    laws = balanced_occupancy_prefix(
        region_length=region_length,
        group_size=group_size,
        limit=hi,
    )
    support_logs = np.array(
        [log_binom(k, r) for r in range(lo, hi + 1)], dtype=np.float64
    )
    log_s = math.log(prime - 1)

    def branch_base(log_z: float, log_v: float, *, field: bool) -> float:
        # The useful saddles are substochastic.  Keeping z+v<=1 prevents
        # exponentially large intermediate slice matrices during diagnostics;
        # it is a restriction of the marker search and therefore remains a
        # valid (possibly weaker) upper bound.
        if math.exp(log_z) + math.exp(log_v) > 1.0:
            return 1e300
        try:
            slices = constraint_uniform_occupancy_trace_transfers(
                length=region_length,
                max_occupied=hi,
                memory=memory,
                output_marker=math.exp(log_z),
                equation_marker=math.exp(log_v),
            )
        except (ArithmeticError, ValueError):
            return 1e300
        terms: list[float] = []
        for r, support_log in zip(range(lo, hi + 1), support_logs):
            region = np.zeros_like(slices[0])
            for occupied, probability in enumerate(laws[r]):
                if probability:
                    region += probability * slices[occupied]
            try:
                trace_log = log_weight_mgf(region, region_count, memory)
            except (ArithmeticError, ValueError):
                return 1e300
            cap_factor = -log_s - r * log_v if field else -(r - 1) * log_v
            terms.append(support_log + trace_log + cap_factor)
        return float(np.logaddexp.reduce(np.asarray(terms))) - cutoff * log_z

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(
        math.exp(-log_s),
        min(0.8, math.sqrt(lo * hi) / n),
    )
    a, z_a, v_a = optimize_substochastic_two_markers(
        objective=lambda point: branch_base(*map(float, point), field=False),
        z_starts=z_starts,
        v_start=v_start,
        lower_log=lower,
        v_floor=0.0,
    )
    b, z_b, v_b = optimize_substochastic_two_markers(
        objective=lambda point: branch_base(*map(float, point), field=True),
        z_starts=z_starts,
        v_start=v_start,
        lower_log=lower,
        v_floor=math.exp(-log_s),
    )
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=(z_a, v_a),
        field_markers=(z_b, v_b),
    )


def constraint_exact_balanced_ec_band_at_markers(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
    structural_markers: tuple[float, float],
    field_markers: tuple[float, float],
) -> ConstraintECBandBound:
    """Evaluate an exact support band at supplied valid cap markers."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")
    log_s = math.log(prime - 1)
    laws = balanced_occupancy_prefix(
        region_length=region_length,
        group_size=group_size,
        limit=hi,
    )

    def evaluate(markers: tuple[float, float], *, field: bool) -> float:
        z, v = markers
        if z <= 0.0 or v <= 0.0 or z > 1.0 or v > 1.0:
            raise ValueError("trace markers must lie in (0,1]")
        if field and v < 1.0 / (prime - 1):
            raise ValueError("field marker must be at least 1/(p-1)")
        slices = constraint_uniform_occupancy_trace_transfers(
            length=region_length,
            max_occupied=hi,
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        terms = []
        for r in range(lo, hi + 1):
            region = np.zeros_like(slices[0])
            for occupied, probability in enumerate(laws[r]):
                if probability:
                    region += probability * slices[occupied]
            cap_factor = -log_s - r * math.log(v) if field else -(r - 1) * math.log(v)
            terms.append(
                log_binom(k, r)
                + log_weight_mgf(region, region_count, memory)
                + cap_factor
            )
        return float(np.logaddexp.reduce(np.asarray(terms))) - cutoff * math.log(z)

    a = evaluate(structural_markers, field=False)
    b = evaluate(field_markers, field=True)
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=structural_markers,
        field_markers=field_markers,
    )


def constraint_balanced_ec_block_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> ConstraintECBandBound:
    """Bound a support block using the positive regional coefficient saddle."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")

    log_s = math.log(prime - 1)
    common = math.log(hi - lo + 1)

    def base(log_z: float, log_v: float, log_x: float) -> float:
        matrix = balanced_constraint_slot_trace_matrix(
            memory=memory,
            group_size=group_size,
            slot_marker=math.exp(log_x),
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        log_b = region_count * log_x + log_v
        endpoint = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        return (
            common
            + endpoint
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    def structural(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) + log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) - log_s

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    midpoint = 0.5 * (lo + hi)
    x_start = midpoint / max(1.0, k - midpoint)
    v_start = max(math.exp(-log_s), min(0.8, midpoint / n))
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    a, z_a, v_a, x_a = optimize_three_marker_modes(
        objective=structural,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=lower,
        upper_x=upper_x,
    )
    b, z_b, v_b, x_b = optimize_three_marker_modes(
        objective=field,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=-log_s,
        upper_x=upper_x,
    )
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=(z_a, v_a, x_a),
        field_markers=(z_b, v_b, x_b),
    )


def constraint_full_support_ec_logbound(
    *, prime: int, n: int, cutoff: int, memory: int, message_weight: int,
) -> ConstraintECBandBound:
    """Bound the full-support family using the occupied constraint trace."""
    r = message_weight
    log_s = math.log(prime - 1)

    def base(log_z: float, log_v: float) -> float:
        _, occupied = nonwrapping_constraint_trace_matrices(
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        return log_weight_mgf(occupied, n, memory) - cutoff * log_z

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(math.exp(-log_s), min(0.8, r / n))
    a, z_a, v_a = optimize_two_marker_modes(
        objective=lambda point: base(*map(float, point))
        - (r - 1) * float(point[1]),
        z_starts=z_starts,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    b, z_b, v_b = optimize_two_marker_modes(
        objective=lambda point: -log_s + base(*map(float, point))
        - r * float(point[1]),
        z_starts=z_starts,
        v_start=v_start,
        lower_z=lower,
        lower_v=-log_s,
    )
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=r,
        support_limit=r,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=(z_a, v_a),
        field_markers=(z_b, v_b),
    )


def singleton_constraint_exact_balanced_ec_band_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> ConstraintECBandBound:
    """Exact sparse-support EC bound retaining singleton coordinates."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")
    support_logs = np.array(
        [log_binom(k, r) for r in range(lo, hi + 1)], dtype=np.float64
    )
    log_s = math.log(prime - 1)

    def branch_base(log_z: float, log_v: float, *, field: bool) -> float:
        if math.exp(log_z) + math.exp(log_v) > 1.0:
            return 1e300
        try:
            slices = singleton_constraint_uniform_support_trace_transfers(
                region_length=region_length,
                group_size=group_size,
                max_weight=hi,
                memory=memory,
                output_marker=math.exp(log_z),
                equation_marker=math.exp(log_v),
            )
        except (ArithmeticError, ValueError):
            return 1e300
        terms: list[float] = []
        for r, support_log in zip(range(lo, hi + 1), support_logs):
            try:
                trace_log = log_weight_mgf(slices[r], region_count, memory)
            except (ArithmeticError, ValueError):
                return 1e300
            cap_factor = -log_s - r * log_v if field else -(r - 1) * log_v
            terms.append(support_log + trace_log + cap_factor)
        return float(np.logaddexp.reduce(np.asarray(terms))) - cutoff * log_z

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(math.exp(-log_s), min(0.8, math.sqrt(lo * hi) / n))
    a, z_a, v_a = optimize_substochastic_two_markers(
        objective=lambda point: branch_base(*map(float, point), field=False),
        z_starts=z_starts,
        v_start=v_start,
        lower_log=lower,
        v_floor=0.0,
    )
    b, z_b, v_b = optimize_substochastic_two_markers(
        objective=lambda point: branch_base(*map(float, point), field=True),
        z_starts=z_starts,
        v_start=v_start,
        lower_log=lower,
        v_floor=math.exp(-log_s),
    )
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=(z_a, v_a),
        field_markers=(z_b, v_b),
    )


def singleton_constraint_exact_balanced_ec_band_at_markers(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
    structural_markers: tuple[float, float],
    field_markers: tuple[float, float],
) -> ConstraintECBandBound:
    """Evaluate the singleton-refined exact band at fixed valid markers."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")
    log_s = math.log(prime - 1)

    def evaluate(markers: tuple[float, float], *, field: bool) -> float:
        z, v = markers
        if z <= 0.0 or v <= 0.0 or z > 1.0 or v > 1.0:
            raise ValueError("trace markers must lie in (0,1]")
        if field and v < 1.0 / (prime - 1):
            raise ValueError("field marker must be at least 1/(p-1)")
        slices = singleton_constraint_uniform_support_trace_transfers(
            region_length=region_length,
            group_size=group_size,
            max_weight=hi,
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        terms = []
        for r in range(lo, hi + 1):
            cap_factor = -log_s - r * math.log(v) if field else -(r - 1) * math.log(v)
            terms.append(
                log_binom(k, r)
                + log_weight_mgf(slices[r], region_count, memory)
                + cap_factor
            )
        return float(np.logaddexp.reduce(np.asarray(terms))) - cutoff * math.log(z)

    a = evaluate(structural_markers, field=False)
    b = evaluate(field_markers, field=True)
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=structural_markers,
        field_markers=field_markers,
    )


def singleton_constraint_exact_balanced_ec_band_at_markers_scaled(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
    structural_markers: tuple[float, float],
    field_markers: tuple[float, float],
) -> ConstraintECBandBound:
    """Underflow-safe exact singleton-refined band at fixed markers."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")
    log_s = math.log(prime - 1)

    def evaluate(markers: tuple[float, float], *, field: bool) -> float:
        z, v = markers
        if z <= 0.0 or v <= 0.0 or z > 1.0 or v > 1.0:
            raise ValueError("trace markers must lie in (0,1]")
        if field and v < 1.0 / (prime - 1):
            raise ValueError("field marker must be at least 1/(p-1)")
        slices, scales = singleton_constraint_uniform_support_trace_transfers_scaled(
            region_length=region_length,
            group_size=group_size,
            max_weight=hi,
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        terms = []
        for r in range(lo, hi + 1):
            cap_factor = -log_s - r * math.log(v) if field else -(r - 1) * math.log(v)
            terms.append(
                log_binom(k, r)
                + log_weight_mgf(slices[r], region_count, memory)
                + region_count * float(scales[r])
                + cap_factor
            )
        return float(np.logaddexp.reduce(np.asarray(terms))) - cutoff * math.log(z)

    a = evaluate(structural_markers, field=False)
    b = evaluate(field_markers, field=True)
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=structural_markers,
        field_markers=field_markers,
    )


def singleton_constraint_balanced_ec_block_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> ConstraintECBandBound:
    """Positive coefficient-saddle EC block bound retaining singletons."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")
    log_s = math.log(prime - 1)
    common = math.log(hi - lo + 1)

    def base(log_z: float, log_v: float, log_x: float) -> float:
        matrix = balanced_singleton_constraint_slot_trace_matrix(
            memory=memory,
            group_size=group_size,
            slot_marker=math.exp(log_x),
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        log_b = region_count * log_x + log_v
        endpoint = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        return (
            common + endpoint
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    def structural(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) + log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) - log_s

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    midpoint = 0.5 * (lo + hi)
    x_start = midpoint / max(1.0, k - midpoint)
    v_start = max(math.exp(-log_s), min(0.8, midpoint / n))
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    a, z_a, v_a, x_a = optimize_three_marker_modes(
        objective=structural,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=lower,
        upper_x=upper_x,
    )
    b, z_b, v_b, x_b = optimize_three_marker_modes(
        objective=field,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=-log_s,
        upper_x=upper_x,
    )
    inverse_log_two = 1.0 / math.log(2.0)
    return ConstraintECBandBound(
        support_start=lo,
        support_limit=hi,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        structural_markers=(z_a, v_a, x_a),
        field_markers=(z_b, v_b, x_b),
    )


def balanced_occupancy_prefix(
    *, region_length: int, group_size: int, limit: int
) -> list[np.ndarray]:
    """Return occupancy-size laws after every draw count through ``limit``."""
    total_slots = region_length * group_size
    if not 0 <= limit <= total_slots:
        raise ValueError("invalid occupancy prefix limit")
    laws = [np.array([1.0], dtype=np.float64)]
    for used in range(limit):
        current = laws[-1]
        following = np.zeros(len(current) + 1, dtype=np.float64)
        remaining = total_slots - used
        for occupied, probability in enumerate(current):
            following[occupied] += (
                probability * (group_size * occupied - used) / remaining
            )
            following[occupied + 1] += (
                probability
                * group_size
                * (region_length - occupied)
                / remaining
            )
        laws.append(following)
    return laws


def coarse_shell_stack_ec_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int,
) -> ShellStackECTerm:
    """Apply one coarse cap after the complete fixed-line EC tail.

    This is the direct regional-shell analogue of the uniform-interleaver
    composition bound.  It is valid but can lose rare-placement information
    because the cap is applied after averaging the regional permutation.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    if r > region_length:
        raise ValueError("exact shell-stack diagnostic is for sparse supports")

    occupancy = balanced_occupancy_size_distribution(
        region_length, group_size, r
    )
    log_s = math.log(prime - 1)
    equation_probability_bound = 1.0 / (prime - 1)

    def objective(log_z: float) -> float:
        try:
            empty, occupied = nonwrapping_trace_matrices(
                prime=prime,
                memory=memory,
                output_marker=math.exp(log_z),
                equation_marker=equation_probability_bound,
            )
            slices = _float_uniform_subset_transfers(
                empty=empty,
                occupied=occupied,
                length=region_length,
                max_occupied=r,
            )
            region = np.zeros_like(slices[0])
            for occupied_count, probability in enumerate(occupancy):
                if probability:
                    region += probability * slices[occupied_count]
            return (
                log_weight_mgf(region, region_count, memory)
                - cutoff * log_z
            )
        except (ArithmeticError, ValueError):
            # Extreme Chernoff markers can underflow a complete endpoint mass.
            # They cannot be minimizers once the endpoint has disappeared.
            return 1.0e300

    lower = -max(120.0, log_s + 12.0)
    value, z = optimize_output_marker_modes(
        objective=objective,
        z_starts=tuple(
            z for z in output_marker_starts(
                prime=prime,
                n=n,
                cutoff=cutoff,
                memory=memory,
                lower_log=lower,
            )
            if z >= 1.0e-6
        ),
        lower_z=lower,
    )
    projective = (r - 1) * log_s + value
    capped = min(0.0, projective)
    inverse_log_two = 1.0 / math.log(2.0)
    return ShellStackECTerm(
        message_weight=r,
        line_tail_log2=value * inverse_log_two,
        projective_union_log2=capped * inverse_log_two,
        total_log2=(log_binom(k, r) + capped) * inverse_log_two,
        output_marker=z,
        cap_saturated=projective >= 0.0,
    )


def tail_subtracted_region_matrix(
    *, prime: int, k: int, region_length: int, group_size: int,
    memory: int, slot_marker: float, output_marker: float,
    equation_marker: float, prefix_limit: int,
    occupancy_laws: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Evaluate the regional polynomial after removing degrees below a limit."""
    if prefix_limit < 0:
        raise ValueError("prefix_limit must be nonnegative")
    if occupancy_laws is None:
        occupancy_laws = balanced_occupancy_prefix(
            region_length=region_length,
            group_size=group_size,
            limit=max(0, prefix_limit - 1),
        )
    empty, occupied = nonwrapping_trace_matrices(
        prime=prime,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    x = slot_marker
    log_normalizer = k * math.log1p(x)
    group_normalizer = math.exp(group_size * math.log1p(x))
    group = (
        empty + math.expm1(group_size * math.log1p(x)) * occupied
    ) / group_normalizer
    full = np.linalg.matrix_power(group, region_length)
    if prefix_limit == 0:
        return full

    slices = uniform_occupancy_trace_transfers(
        length=region_length,
        max_occupied=min(region_length, prefix_limit - 1),
        prime=prime,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    prefix = np.zeros_like(full)
    for draws in range(prefix_limit):
        region = np.zeros_like(full)
        for count, probability in enumerate(occupancy_laws[draws]):
            if probability:
                region += probability * slices[count]
        log_weight = (
            log_binom(k, draws)
            + draws * math.log(x)
            - log_normalizer
        )
        prefix += math.exp(log_weight) * region
    tail = full - prefix
    # The exact difference is coefficientwise nonnegative.  Small negative
    # entries are floating-point cancellation in this diagnostic.
    tolerance = 2e-10 * max(1.0, float(np.max(np.abs(full))))
    if float(np.min(tail)) < -tolerance:
        raise ArithmeticError("unstable tail subtraction")
    return np.maximum(tail, 0.0)


def exact_balanced_ec_projective_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int,
) -> BiregularECTerm:
    """Legacy capped-probability diagnostic; do not use as a theorem bound."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    occupancy = balanced_occupancy_size_distribution(
        region_length, group_size, r
    )

    def region_matrix(log_z: float, log_v: float) -> np.ndarray:
        slices = uniform_occupancy_trace_transfers(
            length=region_length,
            max_occupied=r,
            prime=prime,
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        matrix = np.zeros_like(slices[0])
        for occupied, probability in enumerate(occupancy):
            if probability:
                matrix += probability * slices[occupied]
        return matrix

    def base(log_z: float, log_v: float) -> float:
        return (
            log_weight_mgf(region_matrix(log_z, log_v), region_count, memory)
            - cutoff * log_z
        )

    log_s = math.log(prime - 1)
    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(1e-12, min(0.8, r / n))

    def structural(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return base(log_z, log_v) - (r - 1) * log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return -log_s + base(log_z, log_v) - r * log_v

    a, _, _ = optimize_two_marker_modes(
        objective=structural,
        z_starts=z_starts,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    b, _, _ = optimize_two_marker_modes(
        objective=field,
        z_starts=z_starts,
        v_start=max(math.exp(-log_s), v_start),
        lower_z=lower,
        lower_v=-log_s,
    )
    support_count = log_binom(k, r)
    a += support_count
    b += support_count
    inverse_log_two = 1.0 / math.log(2.0)
    return BiregularECTerm(
        message_weight=r,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
    )


def exact_balanced_ec_prefix_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_limit: int, support_start: int = 1,
) -> BiregularECPrefixBound:
    """Legacy capped-probability diagnostic; do not use as a theorem bound.

    A single pair of output/equation markers is shared by the prefix.  Each
    objective evaluation constructs the uniform-subset transfers only once,
    then evaluates every exact balanced regional coefficient in
    ``[support_start,support_limit]``.  This is both tighter and substantially
    cheaper than applying a relaxed coefficient saddle separately to those
    supports.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    if not 1 <= support_start <= support_limit <= k:
        raise ValueError("support band must lie in [1,k]")

    laws = balanced_occupancy_prefix(
        region_length=region_length,
        group_size=group_size,
        limit=support_limit,
    )
    support_logs = np.array(
        [log_binom(k, r) for r in range(support_start, support_limit + 1)],
        dtype=np.float64,
    )
    log_s = math.log(prime - 1)

    def objective(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        slices = uniform_occupancy_trace_transfers(
            length=region_length,
            max_occupied=support_limit,
            prime=prime,
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        terms: list[float] = []
        for r, support_log in enumerate(support_logs, support_start):
            region = np.zeros_like(slices[0])
            for occupied, probability in enumerate(laws[r]):
                if probability:
                    region += probability * slices[occupied]
            trace_log = log_weight_mgf(region, region_count, memory)
            terms.append(support_log + trace_log - (r - 1) * log_v)
            terms.append(
                support_log + trace_log - log_s - r * log_v
            )
        return (
            float(np.logaddexp.reduce(np.asarray(terms)))
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    value, z, v = optimize_two_marker_modes(
        objective=objective,
        z_starts=z_starts,
        v_start=max(
            math.exp(-log_s),
            min(0.8, math.sqrt(support_start * support_limit) / n),
        ),
        lower_z=lower,
        lower_v=-log_s,
    )
    return BiregularECPrefixBound(
        support_start=support_start,
        support_limit=support_limit,
        output_marker=z,
        equation_marker=v,
        total_log2=value / math.log(2.0),
    )


def uncapped_exact_balanced_ec_band_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> BiregularECPrefixBound:
    """Valid first-moment bound for an exact support band.

    The factor ``(p-1)^(r-1-a)`` union-bounds the projective messages before
    averaging over the random convolution.  Unlike the EA-only cap, this
    factor remains valid when different messages impose different feedback
    equations.  Setting the equation marker to ``1/(p-1)`` evaluates it.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi <= k:
        raise ValueError("support band must lie in [1,k]")

    laws = balanced_occupancy_prefix(
        region_length=region_length,
        group_size=group_size,
        limit=hi,
    )
    support_logs = np.array(
        [log_binom(k, r) for r in range(lo, hi + 1)], dtype=np.float64
    )
    log_s = math.log(prime - 1)
    equation_marker = 1.0 / (prime - 1)

    def objective(log_z: float) -> float:
        slices = uniform_occupancy_trace_transfers(
            length=region_length,
            max_occupied=hi,
            prime=prime,
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=equation_marker,
        )
        terms = []
        for r, support_log in zip(range(lo, hi + 1), support_logs):
            region = np.zeros_like(slices[0])
            for occupied, probability in enumerate(laws[r]):
                if probability:
                    region += probability * slices[occupied]
            terms.append(
                support_log
                + (r - 1) * log_s
                + log_weight_mgf(region, region_count, memory)
            )
        return (
            float(np.logaddexp.reduce(np.asarray(terms)))
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    value, z = optimize_output_marker_modes(
        objective=objective, z_starts=z_starts, lower_z=lower
    )
    return BiregularECPrefixBound(
        support_start=lo,
        support_limit=hi,
        output_marker=z,
        equation_marker=equation_marker,
        total_log2=value / math.log(2.0),
    )


def uncapped_balanced_ec_block_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> tuple[BiregularECPrefixBound, tuple[float, float]]:
    """Valid convex block bound using the uncapped projective first moment."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")

    log_s = math.log(prime - 1)
    equation_marker = 1.0 / (prime - 1)
    common = math.log(hi - lo + 1) - log_s

    def objective(point: np.ndarray) -> float:
        log_z, log_x = map(float, point)
        matrix = balanced_slot_trace_matrix(
            prime=prime,
            memory=memory,
            group_size=group_size,
            slot_marker=math.exp(log_x),
            output_marker=math.exp(log_z),
            equation_marker=equation_marker,
        )
        log_b = region_count * log_x - log_s
        endpoint = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        return (
            common
            + endpoint
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    midpoint = 0.5 * (lo + hi)
    x_start = midpoint / max(1.0, k - midpoint)
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    value, z, x = optimize_output_slot_marker_modes(
        objective=objective,
        z_starts=z_starts,
        x_start=x_start,
        lower_z=lower,
        upper_x=upper_x,
    )
    return (
        BiregularECPrefixBound(
            support_start=lo,
            support_limit=hi,
            output_marker=z,
            equation_marker=equation_marker,
            total_log2=value / math.log(2.0),
        ),
        (x, z),
    )


def uncapped_full_support_ec_logbound(
    *, prime: int, n: int, cutoff: int, memory: int, message_weight: int
) -> tuple[float, float]:
    """Return the uncapped full-support first moment and output marker."""
    log_s = math.log(prime - 1)
    equation_marker = 1.0 / (prime - 1)

    def objective(log_z: float) -> float:
        _, occupied = nonwrapping_trace_matrices(
            prime=prime,
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=equation_marker,
        )
        return (
            (message_weight - 1) * log_s
            + log_weight_mgf(occupied, n, memory)
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    value, z = optimize_output_marker_modes(
        objective=objective, z_starts=z_starts, lower_z=lower
    )
    return value / math.log(2.0), z


def balanced_empty_group_probability(
    *, k: int, group_size: int, message_weight: int
) -> float:
    """Probability that a fixed balanced output group misses the support."""
    r = message_weight
    if not 0 <= r <= k or not 1 <= group_size <= k:
        raise ValueError("invalid balanced occupancy parameters")
    if k - r < group_size:
        return 0.0
    return math.exp(sum(
        math.log(k - r - offset) - math.log(k - offset)
        for offset in range(group_size)
    ))


def nonzero_support_region_saddle_matrix(
    *, prime: int, memory: int, region_length: int, group_size: int,
    slot_marker: float, output_marker: float,
) -> tuple[np.ndarray, float]:
    """Evaluate the positive-degree regional polynomial without subtraction."""
    equation_marker = 1.0 / (prime - 1)
    empty, occupied = nonwrapping_trace_matrices(
        prime=prime,
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    log_group_normalizer = group_size * math.log1p(slot_marker)
    occupied_factor = math.expm1(log_group_normalizer)
    normalizer = math.exp(log_group_normalizer)
    full = (empty + occupied_factor * occupied) / normalizer
    zero_degree = empty / normalizer
    difference = occupied_factor * occupied / normalizer
    size = memory + 1

    def product(
        left: tuple[np.ndarray, float], right: tuple[np.ndarray, float]
    ) -> tuple[np.ndarray, float]:
        if not math.isfinite(left[1]) or not math.isfinite(right[1]):
            return np.zeros((size, size), dtype=np.float64), -math.inf
        matrix, local = normalized(left[0] @ right[0])
        return matrix, left[1] + right[1] + local

    def add(
        left: tuple[np.ndarray, float], right: tuple[np.ndarray, float]
    ) -> tuple[np.ndarray, float]:
        if not math.isfinite(left[1]):
            return right
        if not math.isfinite(right[1]):
            return left
        scale = max(left[1], right[1])
        matrix = (
            math.exp(left[1] - scale) * left[0]
            + math.exp(right[1] - scale) * right[0]
        )
        matrix, local = normalized(matrix)
        return matrix, scale + local

    def singleton(matrix: np.ndarray) -> tuple[np.ndarray, float]:
        return normalized(matrix)

    # A segment stores (P^length, A^length, P^length-A^length), with an
    # independent logarithmic scale for each matrix.  Independent scales are
    # essential because the positive-degree block can be exponentially
    # smaller than the zero-degree block.
    identity = np.eye(size, dtype=np.float64)
    zero = np.zeros((size, size), dtype=np.float64)
    result = (singleton(identity), singleton(identity), (zero, -math.inf))
    power = (singleton(full), singleton(zero_degree), singleton(difference))

    def compose(left, right):
        left_full, left_zero, left_tail = left
        right_full, right_zero, right_tail = right
        return (
            product(left_full, right_full),
            product(left_zero, right_zero),
            add(
                product(left_full, right_tail),
                product(left_tail, right_zero),
            ),
        )

    remaining = region_length
    while remaining:
        if remaining & 1:
            result = compose(result, power)
        remaining >>= 1
        if remaining:
            power = compose(power, power)
    return result[2]


def uncapped_nonzero_region_ec_block_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> tuple[BiregularECPrefixBound, tuple[float, float]]:
    """Convex block bound after removing each region's zero-support term."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")
    log_s = math.log(prime - 1)
    common = math.log(hi - lo + 1) - log_s

    def objective(point: np.ndarray) -> float:
        log_z, log_x = map(float, point)
        x = math.exp(log_x)
        try:
            region, region_log_scale = nonzero_support_region_saddle_matrix(
                prime=prime,
                memory=memory,
                region_length=region_length,
                group_size=group_size,
                slot_marker=x,
                output_marker=math.exp(log_z),
            )
            trace_log = (
                log_weight_mgf(region, region_count, memory)
                + region_count * region_log_scale
            )
        except (OverflowError, ValueError):
            return 1e300
        log_b = region_count * log_x - log_s
        endpoint = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        return (
            common
            + endpoint
            + region_count * k * math.log1p(x)
            + trace_log
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    midpoint = 0.5 * (lo + hi)
    x_start = midpoint / max(1.0, k - midpoint)
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    value, z, x = optimize_output_slot_marker_modes(
        objective=objective,
        z_starts=z_starts,
        x_start=x_start,
        lower_z=lower,
        upper_x=upper_x,
    )
    return (
        BiregularECPrefixBound(
            support_start=lo,
            support_limit=hi,
            output_marker=z,
            equation_marker=1.0 / (prime - 1),
            total_log2=value / math.log(2.0),
        ),
        (x, z),
    )


def uncapped_balanced_ec_dense_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int,
) -> tuple[float, float]:
    """Uncapped first moment from balanced empty-set domination."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    if not 1 <= r <= k:
        raise ValueError("message_weight must lie in [1,k]")
    beta = balanced_empty_group_probability(
        k=k, group_size=group_size, message_weight=r
    )
    log_s = math.log(prime - 1)
    equation_marker = 1.0 / (prime - 1)

    def objective(log_z: float) -> float:
        empty, occupied = nonwrapping_trace_matrices(
            prime=prime,
            memory=memory,
            output_marker=math.exp(log_z),
            equation_marker=equation_marker,
        )
        matrix = beta * empty + occupied
        return (
            log_binom(k, r)
            + (r - 1) * log_s
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    value, z = optimize_output_marker_modes(
        objective=objective, z_starts=z_starts, lower_z=lower
    )
    return value / math.log(2.0), z


def balanced_ec_projective_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int,
) -> BiregularECTerm:
    """Legacy capped-probability diagnostic; do not use as a theorem bound."""
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    if not 1 <= message_weight <= k:
        raise ValueError("message_weight must lie in [1,k]")

    r = message_weight
    log_s = math.log(prime - 1)
    lower = -max(120.0, log_s + 12.0)
    # Rare runs of ``memory`` zero outputs can turn the convolution off and
    # create a very low-weight component.  That component places the useful
    # global saddle close to one, even though the active-state saddle alone
    # would be near 1/(p-1).
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(1e-12, min(0.8, r / n))

    def finish(base_objective, *, has_slot_marker: bool) -> BiregularECTerm:
        def structural(point: np.ndarray) -> float:
            return base_objective(*map(float, point)) - (r - 1) * float(point[1])

        def field(point: np.ndarray) -> float:
            return -log_s + base_objective(*map(float, point)) - r * float(point[1])

        if has_slot_marker:
            x_start = r / max(1.0, k - r)
            upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
            a, _, _, _ = optimize_three_marker_modes(
                objective=structural,
                z_starts=z_starts,
                v_start=v_start,
                x_start=x_start,
                lower_z=lower,
                lower_v=lower,
                upper_x=upper_x,
            )
            b, _, _, _ = optimize_three_marker_modes(
                objective=field,
                z_starts=z_starts,
                v_start=max(math.exp(-log_s), v_start),
                x_start=x_start,
                lower_z=lower,
                lower_v=-log_s,
                upper_x=upper_x,
            )
        else:
            a, _, _ = optimize_two_marker_modes(
                objective=structural,
                z_starts=z_starts,
                v_start=v_start,
                lower_z=lower,
                lower_v=lower,
            )
            b, _, _ = optimize_two_marker_modes(
                objective=field,
                z_starts=z_starts,
                v_start=max(math.exp(-log_s), v_start),
                lower_z=lower,
                lower_v=-log_s,
            )
        support_count = log_binom(k, r)
        a += support_count
        b += support_count
        inverse_log_two = 1.0 / math.log(2.0)
        return BiregularECTerm(
            message_weight=r,
            structural_log2=a * inverse_log_two,
            field_log2=b * inverse_log_two,
            total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
        )

    if r == k:
        def full_base(log_z: float, log_v: float) -> float:
            _, occupied = nonwrapping_trace_matrices(
                prime=prime,
                memory=memory,
                output_marker=math.exp(log_z),
                equation_marker=math.exp(log_v),
            )
            return (
                log_weight_mgf(occupied, n, memory)
                - cutoff * log_z
            )

        return finish(full_base, has_slot_marker=False)

    normalization = region_count * log_binom(k, r)

    def coefficient_base(log_z: float, log_v: float, log_x: float) -> float:
        matrix = balanced_slot_trace_matrix(
            prime=prime,
            memory=memory,
            group_size=group_size,
            slot_marker=math.exp(log_x),
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        return (
            -normalization
            - region_count * r * log_x
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    return finish(coefficient_base, has_slot_marker=True)


def balanced_ec_block_logbound(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, support_start: int, support_limit: int,
) -> tuple[BiregularECPrefixBound, dict[str, tuple[float, float, float]]]:
    """Legacy capped-probability diagnostic; do not use as a theorem bound.

    This is the block form of the positive regional coefficient saddle.  The
    dependence on the support size is log-convex, so its maximum on the block
    occurs at an endpoint; multiplying by the block length then bounds the
    whole interval.
    """
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")

    log_s = math.log(prime - 1)
    common = math.log(hi - lo + 1)

    def base(log_z: float, log_v: float, log_x: float) -> float:
        matrix = balanced_slot_trace_matrix(
            prime=prime,
            memory=memory,
            group_size=group_size,
            slot_marker=math.exp(log_x),
            output_marker=math.exp(log_z),
            equation_marker=math.exp(log_v),
        )
        log_b = region_count * log_x + log_v
        endpoint = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        return (
            common
            + endpoint
            + log_weight_mgf(matrix, n, memory)
            - cutoff * log_z
        )

    def structural(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) + log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) - log_s

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    midpoint = 0.5 * (lo + hi)
    x_start = midpoint / max(1.0, k - midpoint)
    v_start = max(math.exp(-log_s), min(0.8, midpoint / n))
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    a, z_a, v_a, x_a = optimize_three_marker_modes(
        objective=structural,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=lower,
        upper_x=upper_x,
    )
    b, z_b, v_b, x_b = optimize_three_marker_modes(
        objective=field,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=-log_s,
        upper_x=upper_x,
    )
    value = float(np.logaddexp(a, b))
    return (
        BiregularECPrefixBound(
            support_start=lo,
            support_limit=hi,
            output_marker=z_a,
            equation_marker=v_a,
            total_log2=value / math.log(2.0),
        ),
        {
            "structural": (x_a, z_a, v_a),
            "field": (x_b, z_b, v_b),
        },
    )


def tail_subtracted_ec_projective_logterm(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, message_weight: int, prefix_limit: int,
) -> BiregularECTerm:
    """Optimize a saddle after removing impossible lower support degrees."""
    if not 1 <= prefix_limit <= message_weight:
        raise ValueError("prefix_limit must lie in [1,message_weight]")
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    group_size = k // region_length
    r = message_weight
    log_s = math.log(prime - 1)
    log_choose = log_binom(k, r)
    occupancy_laws = balanced_occupancy_prefix(
        region_length=region_length,
        group_size=group_size,
        limit=prefix_limit - 1,
    )

    def base(log_z: float, log_v: float, log_x: float) -> float:
        x = math.exp(log_x)
        try:
            region = tail_subtracted_region_matrix(
                prime=prime,
                k=k,
                region_length=region_length,
                group_size=group_size,
                memory=memory,
                slot_marker=x,
                output_marker=math.exp(log_z),
                equation_marker=math.exp(log_v),
                prefix_limit=prefix_limit,
                occupancy_laws=occupancy_laws,
            )
            mgf = log_weight_mgf(region, region_count, memory)
        except (ArithmeticError, ValueError):
            return 1e300
        return (
            region_count
            * (k * math.log1p(x) - r * log_x - log_choose)
            + mgf
            - cutoff * log_z
        )

    def structural(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return base(log_z, log_v, log_x) - (r - 1) * log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        return -log_s + base(log_z, log_v, log_x) - r * log_v

    lower = -max(120.0, log_s + 12.0)
    z_starts = output_marker_starts(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        lower_log=lower,
    )
    v_start = max(1e-12, min(0.8, r / n))
    x_start = r / max(1.0, k - r)
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    a, _, _, _ = optimize_three_marker_modes(
        objective=structural,
        z_starts=z_starts,
        v_start=v_start,
        x_start=x_start,
        lower_z=lower,
        lower_v=lower,
        upper_x=upper_x,
    )
    b, _, _, _ = optimize_three_marker_modes(
        objective=field,
        z_starts=z_starts,
        v_start=max(math.exp(-log_s), v_start),
        x_start=x_start,
        lower_z=lower,
        lower_v=-log_s,
        upper_x=upper_x,
    )
    a += log_choose
    b += log_choose
    inverse_log_two = 1.0 / math.log(2.0)
    return BiregularECTerm(
        message_weight=r,
        structural_log2=a * inverse_log_two,
        field_log2=b * inverse_log_two,
        total_log2=float(np.logaddexp(a, b)) * inverse_log_two,
    )


def candidate_parameters(target_length: int, degree: int) -> tuple[int, int, int]:
    """Return (k,n,ell) at rate one half with ``n`` divisible by degree."""
    if degree < 2 or degree % 2:
        raise ValueError("two-sided rate-half degree must be positive and even")
    region_length = target_length // degree
    n = degree * region_length
    k = n // 2
    if k % region_length:
        raise AssertionError("balanced rate-half construction is inconsistent")
    return k, n, region_length


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime", type=int,
        default=170141183460469231731687303715884105727,
    )
    parser.add_argument("--target-length", type=int, default=2**21)
    parser.add_argument("--degree", type=int, default=28)
    parser.add_argument("--memory", type=int, default=3)
    parser.add_argument("--exact-limit", type=int, default=8)
    parser.add_argument("--gv-scale", type=float, default=1.0)
    args = parser.parse_args()

    k, n, ell = candidate_parameters(args.target_length, args.degree)
    delta = gv_distance(args.prime, 0.5) * args.gv_scale
    cutoff = math.floor(n * delta)
    weights = diagnostic_weights(k)
    terms: list[tuple[int, ConstraintECBandBound]] = []
    for r in weights:
        if r <= args.exact_limit:
            bound = constraint_exact_balanced_ec_band_logbound(
                prime=args.prime,
                k=k,
                n=n,
                cutoff=cutoff,
                region_count=args.degree,
                memory=args.memory,
                support_start=r,
                support_limit=r,
            )
        elif r == k:
            bound = constraint_full_support_ec_logbound(
                prime=args.prime,
                n=n,
                cutoff=cutoff,
                memory=args.memory,
                message_weight=r,
            )
        else:
            bound = constraint_balanced_ec_block_logbound(
                prime=args.prime,
                k=k,
                n=n,
                cutoff=cutoff,
                region_count=args.degree,
                memory=args.memory,
                support_start=r,
                support_limit=r,
            )
        terms.append((r, bound))
    worst_r, worst = max(terms, key=lambda item: item[1].total_log2)
    print(
        f"p={args.prime} d={args.degree} q={args.degree//2} "
        f"memory={args.memory} ell={ell} n={n} k={k} cutoff={cutoff} "
        f"delta={cutoff/n:.12f}"
    )
    print("r\ttrace\tpoint-mass\tcombined")
    for r, term in terms:
        dense = balanced_point_mass_logterm(
            prime=args.prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=args.degree,
            message_weight=r,
        )
        combined = min(term.total_log2, dense)
        if r <= 64 or combined > -100.0:
            print(
                f"{r}\t{term.total_log2:.6f}\t"
                f"{dense:.6f}\t{combined:.6f}"
            )
    print(
        f"sampled worst: r={worst_r}, "
        f"log2_bound={worst.total_log2:.6f}"
    )


if __name__ == "__main__":
    main()
