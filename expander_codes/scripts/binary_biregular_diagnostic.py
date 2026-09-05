#!/usr/bin/env python3
"""Floating-point diagnostics for binary two-sided regular EA and EC.

Each region is a uniform partition of the ``k`` left vertices into groups of
size ``right_degree``.  For a fixed support of size ``r``, the occupied slots
therefore form a uniform ``r``-subset.  A right coordinate is one precisely
when its group contains an odd number of selected slots.

Small supports use the exact regional parity-shell distribution.  Larger
supports use a positive coefficient saddle.  This program selects candidate
parameters; it is not an interval certificate.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

from expander_bounds import ea_matrix, ec_matrix, log_binom, log_weight_mgf
from regular_ec_diagnostic import (
    diagnostic_weights,
    log_hamming_ball_geometric,
    uniform_slice_transfer_matrices,
)
from regular_ea_diagnostic import accumulator_uniform_slice_transfers


@dataclass(frozen=True)
class BiregularBinaryTerm:
    message_weight: int
    input_probability: float
    input_marker: float
    output_marker: float
    log2_bound: float


@dataclass(frozen=True)
class DegreeFiveDenseResult:
    support_start: int
    support_limit: int
    worst_support: int
    worst_log2_term: float
    summed_log2_bound: float


@dataclass(frozen=True)
class BiregularBinaryBlock:
    support_start: int
    support_limit: int
    input_marker: float
    output_marker: float
    log2_bound: float


@dataclass(frozen=True)
class ExactBiregularBinaryBlock:
    support_start: int
    support_limit: int
    output_marker: float
    worst_support: int
    worst_log2_term: float
    summed_log2_bound: float
    log2_terms: tuple[float, ...]


def _multiply_polynomials(
    left: list[int], right: list[int], limit: int
) -> list[int]:
    """Multiply integer polynomials and discard degrees above ``limit``."""
    result = [0] * min(limit + 1, len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if not a:
            continue
        stop = min(len(right), limit - i + 1)
        for j in range(stop):
            if right[j]:
                result[i + j] += a * right[j]
    return result


def _power_polynomial(base: list[int], exponent: int, limit: int) -> list[int]:
    """Return ``base**exponent`` truncated above ``limit``."""
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    result = [1]
    power = base[: limit + 1]
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _multiply_polynomials(result, power, limit)
        remaining >>= 1
        if remaining:
            power = _multiply_polynomials(power, power, limit)
    return result


def biregular_parity_shell_counts(
    *, left_vertices: int, right_degree: int, support_size: int
) -> tuple[int, dict[int, int]]:
    """Return the denominator and exact parity-shell counts for one region."""
    if right_degree < 1 or left_vertices % right_degree:
        raise ValueError("right_degree must divide left_vertices")
    if not 0 <= support_size <= left_vertices:
        raise ValueError("support_size must lie in [0,left_vertices]")
    groups = left_vertices // right_degree
    r = support_size
    even = [
        math.comb(right_degree, degree) if degree % 2 == 0 else 0
        for degree in range(min(right_degree, r) + 1)
    ]
    odd = [
        math.comb(right_degree, degree) if degree % 2 == 1 else 0
        for degree in range(min(right_degree, r) + 1)
    ]
    denominator = math.comb(left_vertices, r)
    counts: dict[int, int] = {}
    for parity_weight in range(min(groups, r) + 1):
        if parity_weight % 2 != r % 2:
            continue
        odd_power = _power_polynomial(odd, parity_weight, r)
        even_power = _power_polynomial(even, groups - parity_weight, r)
        coefficient = _multiply_polynomials(odd_power, even_power, r)
        numerator = (
            math.comb(groups, parity_weight)
            * (coefficient[r] if r < len(coefficient) else 0)
        )
        if numerator:
            counts[parity_weight] = numerator
    return denominator, counts


def degree_three_parity_shell_counts(
    *, left_vertices: int, support_size: int,
) -> tuple[int, dict[int, int]]:
    """Return the right-degree-three shell counts by a direct finite sum."""
    if left_vertices % 3:
        raise ValueError("three must divide left_vertices")
    if not 0 <= support_size <= left_vertices:
        raise ValueError("support_size must lie in [0,left_vertices]")
    groups = left_vertices // 3
    r = support_size
    denominator = math.comb(left_vertices, r)
    counts: dict[int, int] = {}
    for parity_weight in range(r % 2, min(groups, r) + 1, 2):
        half_excess = (r - parity_weight) // 2
        coefficient = 0
        lo = max(0, half_excess - (groups - parity_weight))
        hi = min(parity_weight, half_excess)
        for triple_groups in range(lo, hi + 1):
            double_groups = half_excess - triple_groups
            coefficient += (
                math.comb(parity_weight, triple_groups)
                * math.comb(groups - parity_weight, double_groups)
                * 3 ** (
                    parity_weight
                    + double_groups
                    - triple_groups
                )
            )
        if coefficient:
            counts[parity_weight] = math.comb(groups, parity_weight) * coefficient
    return denominator, counts


def biregular_parity_shell_distribution(
    *, left_vertices: int, right_degree: int, support_size: int
) -> dict[int, float]:
    """Return the exact parity-weight law in one two-sided regular region."""
    denominator, counts = biregular_parity_shell_counts(
        left_vertices=left_vertices,
        right_degree=right_degree,
        support_size=support_size,
    )
    return {weight: count / denominator for weight, count in counts.items()}


def parity_probability(right_degree: int, input_marker: float) -> float:
    """Odd-group probability under independently marked slots."""
    if right_degree < 1 or input_marker <= 0.0:
        raise ValueError("invalid degree or marker")
    ratio = (1.0 - input_marker) / (1.0 + input_marker)
    return 0.5 * (1.0 - ratio**right_degree)


def central_fourier_l1_log2(*, left_vertices: int, right_degree: int) -> float:
    """Fourier L1 norm for support ``(left_vertices-1)/2``.

    The formula uses the central Krawtchouk identity for odd
    ``left_vertices``.  It bounds a regional point mass by
    ``2**(-region_length)`` times the returned norm.
    """
    k = left_vertices
    if k % 2 != 1 or right_degree % 2 != 1 or k % right_degree:
        raise ValueError("central identity requires odd compatible parameters")
    region_length = k // right_degree
    half = (k - 1) // 2
    groups = np.arange(region_length + 1, dtype=np.float64)
    slots = right_degree * groups

    def array_log_binom(n: int, values: np.ndarray) -> np.ndarray:
        return gammaln(n + 1.0) - gammaln(values + 1.0) - gammaln(n - values + 1.0)

    log_terms = (
        array_log_binom(region_length, groups)
        + array_log_binom(half, np.floor(slots / 2.0))
        - array_log_binom(k, slots)
    )
    return float(logsumexp(log_terms) / math.log(2.0))


def central_dense_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int
) -> float:
    """Invertibility bound for either central support of odd ``k``."""
    region_length = k // right_degree
    n = left_degree * region_length
    r = (k - 1) // 2
    return (
        log_binom(k, r)
        + log_hamming_ball_geometric(n, cutoff)
        - n * math.log(2.0)
        + left_degree
        * central_fourier_l1_log2(
            left_vertices=k, right_degree=right_degree
        )
        * math.log(2.0)
    ) / math.log(2.0)


def degree_five_conditional_variances(input_probability: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Variances of half-counts conditioned on even or odd group parity."""
    p = input_probability
    one_minus_p = 1.0 - p
    odd_probability = 0.5 * (1.0 - (1.0 - 2.0 * p) ** 5)
    even_probability = 1.0 - odd_probability

    even_one = 10.0 * p**2 * one_minus_p**3
    even_two = 5.0 * p**4 * one_minus_p
    even_mean = (even_one + 2.0 * even_two) / even_probability
    even_variance = (
        (even_one + 4.0 * even_two) / even_probability - even_mean**2
    )

    odd_one = 10.0 * p**3 * one_minus_p**2
    odd_two = p**5
    odd_mean = (odd_one + 2.0 * odd_two) / odd_probability
    odd_variance = (
        (odd_one + 4.0 * odd_two) / odd_probability - odd_mean**2
    )
    return even_variance, odd_variance


def degree_three_conditional_variances(
    input_probability: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Variances of half-counts conditioned on even or odd group parity."""
    p = input_probability
    even_success = 3.0 * p**2 / (4.0 * p**2 - 2.0 * p + 1.0)
    odd_success = p**2 / (4.0 * p**2 - 6.0 * p + 3.0)
    return (
        even_success * (1.0 - even_success),
        odd_success * (1.0 - odd_success),
    )


def degree_three_dense_scan(
    *, k: int, left_degree: int, cutoff: int,
    support_start: int, support_limit: int,
) -> DegreeFiveDenseResult:
    """Scan the degree-three conditional local-limit point-mass bound."""
    if k % 3:
        raise ValueError("degree three must divide k")
    if not 1 <= support_start <= support_limit < k:
        raise ValueError("invalid dense support interval")
    region_length = k // 3
    n = left_degree * region_length
    supports = np.arange(support_start, support_limit + 1, dtype=np.float64)
    p = supports / k
    one_minus_p = 1.0 - p
    odd_probability = 0.5 * (1.0 - (1.0 - 2.0 * p) ** 3)
    even_probability = 1.0 - odd_probability
    even_variance, odd_variance = degree_three_conditional_variances(p)
    minimum_variance = np.minimum(even_variance, odd_variance)

    log_choose = (
        gammaln(k + 1.0)
        - gammaln(supports + 1.0)
        - gammaln(k - supports + 1.0)
    )
    log_binomial_mass = (
        log_choose
        + supports * np.log(p)
        + (k - supports) * np.log(one_minus_p)
    )
    log_conditional_mass = np.minimum(
        0.0,
        0.5
        * (
            math.log(math.pi / 8.0)
            - np.log(region_length * minimum_variance)
        ),
    )
    log_regional_point_mass = (
        region_length * np.log(np.maximum(odd_probability, even_probability))
        + log_conditional_mass
        - log_binomial_mass
    )
    log_ball = log_hamming_ball_geometric(n, cutoff)
    log_terms = log_choose + log_ball + left_degree * log_regional_point_mass
    worst_index = int(np.argmax(log_terms))
    worst_log2 = float(log_terms[worst_index] / math.log(2.0))
    return DegreeFiveDenseResult(
        support_start=support_start,
        support_limit=support_limit,
        worst_support=int(supports[worst_index]),
        worst_log2_term=worst_log2,
        summed_log2_bound=worst_log2
        + math.log2(support_limit - support_start + 1),
    )


def degree_five_dense_scan(
    *, k: int, left_degree: int, cutoff: int,
    support_start: int, support_limit: int,
) -> DegreeFiveDenseResult:
    """Scan the proved degree-five local-limit point-mass bound."""
    if k % 5:
        raise ValueError("degree five must divide k")
    if not 1 <= support_start <= support_limit < k:
        raise ValueError("invalid dense support interval")
    region_length = k // 5
    n = left_degree * region_length
    supports = np.arange(support_start, support_limit + 1, dtype=np.float64)
    p = supports / k
    one_minus_p = 1.0 - p
    odd_probability = 0.5 * (1.0 - (1.0 - 2.0 * p) ** 5)
    even_probability = 1.0 - odd_probability
    even_variance, odd_variance = degree_five_conditional_variances(p)
    minimum_variance = np.minimum(even_variance, odd_variance)

    log_choose = (
        gammaln(k + 1.0)
        - gammaln(supports + 1.0)
        - gammaln(k - supports + 1.0)
    )
    log_binomial_mass = (
        log_choose
        + supports * np.log(p)
        + (k - supports) * np.log(one_minus_p)
    )
    log_conditional_mass = np.minimum(
        0.0,
        0.5
        * (
            math.log(math.pi / 8.0)
            - np.log(region_length * minimum_variance)
        ),
    )
    log_regional_point_mass = (
        region_length * np.log(np.maximum(odd_probability, even_probability))
        + log_conditional_mass
        - log_binomial_mass
    )
    log_ball = log_hamming_ball_geometric(n, cutoff)
    log_terms = log_choose + log_ball + left_degree * log_regional_point_mass
    worst_index = int(np.argmax(log_terms))
    worst_log2 = float(log_terms[worst_index] / math.log(2.0))
    return DegreeFiveDenseResult(
        support_start=support_start,
        support_limit=support_limit,
        worst_support=int(supports[worst_index]),
        worst_log2_term=worst_log2,
        summed_log2_bound=worst_log2
        + math.log2(support_limit - support_start + 1),
    )


def _exact_region_matrix(
    *, code: str, region_length: int, support_size: int,
    shell: dict[int, float], output_marker: float, memory: int,
) -> np.ndarray:
    if code == "ea":
        slices = accumulator_uniform_slice_transfers(
            length=region_length,
            max_weight=support_size,
            output_marker=output_marker,
        )
    else:
        slices = uniform_slice_transfer_matrices(
            length=region_length,
            max_weight=support_size,
            output_marker=output_marker,
            memory=memory,
        )
    region = np.zeros_like(slices[0])
    for weight, probability in shell.items():
        region += probability * slices[weight]
    return region


def exact_biregular_logterm(
    *, code: str, k: int, left_degree: int, right_degree: int,
    cutoff: int, memory: int, message_weight: int,
) -> BiregularBinaryTerm:
    """Optimize the exact regional transfer for one small support."""
    if code not in {"ea", "ec"}:
        raise ValueError("code must be ea or ec")
    if k % right_degree or left_degree * k % right_degree:
        raise ValueError("incompatible two-sided degrees")
    region_length = k // right_degree
    n = left_degree * region_length
    shell = biregular_parity_shell_distribution(
        left_vertices=k,
        right_degree=right_degree,
        support_size=message_weight,
    )

    def objective(log_z: float) -> float:
        matrix = _exact_region_matrix(
            code=code,
            region_length=region_length,
            support_size=message_weight,
            shell=shell,
            output_marker=math.exp(log_z),
            memory=memory,
        )
        start_state = 0 if code == "ea" else memory
        return (
            log_weight_mgf(matrix, left_degree, start_state)
            - cutoff * log_z
        )

    optimum = minimize_scalar(
        objective,
        bounds=(-40.0, 0.0),
        method="bounded",
        options={"xatol": 2e-9, "maxiter": 160},
    )
    mean_weight = sum(weight * probability for weight, probability in shell.items())
    return BiregularBinaryTerm(
        message_weight=message_weight,
        input_probability=mean_weight / region_length,
        input_marker=math.nan,
        output_marker=math.exp(float(optimum.x)),
        log2_bound=(log_binom(k, message_weight) + float(optimum.fun))
        / math.log(2.0),
    )


def exact_biregular_fixed_marker_block(
    *, code: str, k: int, left_degree: int, right_degree: int,
    cutoff: int, memory: int, support_start: int, support_limit: int,
    output_marker: float,
) -> ExactBiregularBinaryBlock:
    """Evaluate an exact support block while reusing one slice family.

    The result is a floating-point diagnostic.  Every support in the block
    uses the same positive output marker, but its regional parity-shell law is
    evaluated separately.  Reusing the uniform-slice transfers is much faster
    than rebuilding them for every support.
    """
    if code not in {"ea", "ec"}:
        raise ValueError("code must be ea or ec")
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    if not 1 <= support_start <= support_limit < k:
        raise ValueError("support block must lie in [1,k-1]")
    if not 0.0 < output_marker <= 1.0:
        raise ValueError("output_marker must lie in (0,1]")
    region_length = k // right_degree
    n = left_degree * region_length
    if code == "ea":
        slices = accumulator_uniform_slice_transfers(
            length=region_length,
            max_weight=support_limit,
            output_marker=output_marker,
        )
        start_state = 0
    else:
        slices = uniform_slice_transfer_matrices(
            length=region_length,
            max_weight=support_limit,
            output_marker=output_marker,
            memory=memory,
        )
        start_state = memory

    log2_terms: list[float] = []
    for r in range(support_start, support_limit + 1):
        shell = biregular_parity_shell_distribution(
            left_vertices=k,
            right_degree=right_degree,
            support_size=r,
        )
        region = np.zeros_like(slices[0])
        for weight, probability in shell.items():
            region += probability * slices[weight]
        log_term = (
            log_binom(k, r)
            + log_weight_mgf(region, left_degree, start_state)
            - cutoff * math.log(output_marker)
        )
        log2_terms.append(log_term / math.log(2.0))

    worst_index = int(np.argmax(log2_terms))
    worst = log2_terms[worst_index]
    summed = worst + math.log2(sum(2.0 ** (value - worst) for value in log2_terms))
    return ExactBiregularBinaryBlock(
        support_start=support_start,
        support_limit=support_limit,
        output_marker=output_marker,
        worst_support=support_start + worst_index,
        worst_log2_term=worst,
        summed_log2_bound=summed,
        log2_terms=tuple(log2_terms),
    )


def saddle_biregular_logterm(
    *, code: str, k: int, left_degree: int, right_degree: int,
    cutoff: int, memory: int, message_weight: int,
) -> BiregularBinaryTerm:
    """Optimize the positive regional coefficient bound for one support."""
    if code not in {"ea", "ec"}:
        raise ValueError("code must be ea or ec")
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    region_length = k // right_degree
    n = left_degree * region_length
    if left_degree * k != n * right_degree:
        raise ValueError("left and right degrees do not balance")
    r = message_weight
    if not 1 <= r < k:
        raise ValueError("coefficient saddle requires 1 <= r < k")
    log_choose = log_binom(k, r)

    def objective(point: np.ndarray) -> float:
        log_x, log_z = map(float, point)
        x = math.exp(log_x)
        z = math.exp(log_z)
        q = parity_probability(right_degree, x)
        matrix = ea_matrix(q, z) if code == "ea" else ec_matrix(q, z, memory, True)
        start_state = 0 if code == "ea" else memory
        return (
            (1.0 - left_degree) * log_choose
            + left_degree * (k * math.log1p(x) - r * log_x)
            + log_weight_mgf(matrix, n, start_state)
            - cutoff * log_z
        )

    x_center = math.log(r / (k - r))
    z_center = math.log(max(1e-8, min(0.95, cutoff / (n - cutoff))))
    candidates = []
    for dx in (-1.0, 0.0, 1.0):
        for log_z in (z_center, -0.5, -0.05):
            candidates.append(
                minimize(
                    objective,
                    np.array([x_center + dx, log_z]),
                    method="L-BFGS-B",
                    bounds=((-40.0, 40.0), (-40.0, 0.0)),
                    options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
                )
            )
    optimum = min(candidates, key=lambda result: float(result.fun))
    x = math.exp(float(optimum.x[0]))
    return BiregularBinaryTerm(
        message_weight=r,
        input_probability=parity_probability(right_degree, x),
        input_marker=x,
        output_marker=math.exp(float(optimum.x[1])),
        log2_bound=float(optimum.fun) / math.log(2.0),
    )


def saddle_biregular_block_logterm(
    *, code: str, k: int, left_degree: int, right_degree: int,
    cutoff: int, memory: int, support_start: int, support_limit: int,
    multistart: bool = True,
) -> BiregularBinaryBlock:
    """Bound a support block by convexity at fixed positive markers."""
    if code not in {"ea", "ec"}:
        raise ValueError("code must be ea or ec")
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    lo, hi = support_start, support_limit
    if not 1 <= lo <= hi < k:
        raise ValueError("support block must lie in [1,k-1]")
    region_length = k // right_degree
    n = left_degree * region_length

    def objective(point: np.ndarray) -> float:
        log_x, log_z = map(float, point)
        x = math.exp(log_x)
        z = math.exp(log_z)
        q = parity_probability(right_degree, x)
        matrix = ea_matrix(q, z) if code == "ea" else ec_matrix(q, z, memory, True)
        start_state = 0 if code == "ea" else memory
        endpoint = max(
            (1.0 - left_degree) * log_binom(k, lo) - left_degree * lo * log_x,
            (1.0 - left_degree) * log_binom(k, hi) - left_degree * hi * log_x,
        )
        return (
            math.log(hi - lo + 1)
            + endpoint
            + left_degree * k * math.log1p(x)
            + log_weight_mgf(matrix, n, start_state)
            - cutoff * log_z
        )

    midpoint = 0.5 * (lo + hi)
    x_center = math.log(midpoint / (k - midpoint))
    z_center = math.log(max(1e-8, min(0.95, cutoff / (n - cutoff))))
    candidates = []
    starts = (
        [(dx, log_z) for dx in (-0.5, 0.0, 0.5)
         for log_z in (z_center, -0.5, -0.05)]
        if multistart
        else [(0.0, z_center)]
    )
    for dx, log_z in starts:
        candidates.append(
            minimize(
                objective,
                np.array([x_center + dx, log_z]),
                method="L-BFGS-B",
                bounds=((-40.0, 40.0), (-40.0, 0.0)),
                options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500},
            )
        )
    optimum = min(candidates, key=lambda result: float(result.fun))
    return BiregularBinaryBlock(
        support_start=lo,
        support_limit=hi,
        input_marker=math.exp(float(optimum.x[0])),
        output_marker=math.exp(float(optimum.x[1])),
        log2_bound=float(optimum.fun) / math.log(2.0),
    )


def geometric_blocks(start: int, limit: int, relative_width: float) -> list[tuple[int, int]]:
    """Partition an integer interval into short multiplicative blocks."""
    if not 1 <= start <= limit or not 0.0 < relative_width < 1.0:
        raise ValueError("invalid block interval or relative width")
    blocks = []
    lo = start
    while lo <= limit:
        hi = min(limit, max(lo, int(math.floor(lo * (1.0 + relative_width)))))
        blocks.append((lo, hi))
        lo = hi + 1
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code", choices=("ea", "ec"), default="ec")
    parser.add_argument("--k", type=int, default=1_048_567)
    parser.add_argument("--left-degree", type=int, default=26)
    parser.add_argument("--right-degree", type=int, default=13)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=9)
    parser.add_argument("--exact-limit", type=int, default=8)
    parser.add_argument("--scan-blocks", action="store_true")
    parser.add_argument("--block-relative-width", type=float, default=0.02)
    args = parser.parse_args()
    if args.k % args.right_degree:
        raise ValueError("right-degree must divide k")
    region_length = args.k // args.right_degree
    n = args.left_degree * region_length
    if args.left_degree * args.k != n * args.right_degree:
        raise ValueError("degree equation is not satisfied")
    if args.right_degree % 2 == 0:
        raise ValueError("even binary right degree puts the all-one message in the kernel")

    print(
        f"code={args.code} k={args.k} n={n} rate={args.k/n:.12f} "
        f"left_degree={args.left_degree} right_degree={args.right_degree} "
        f"region_length={region_length} memory={args.memory} "
        f"cutoff={args.cutoff} relative_cutoff={args.cutoff/n:.12f}"
    )
    print("exact regional terms")
    exact_terms = []
    for r in range(1, args.exact_limit + 1):
        term = exact_biregular_logterm(
            code=args.code,
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
            memory=args.memory,
            message_weight=r,
        )
        exact_terms.append(term)
        print(
            f"{r}\tq={term.input_probability:.9g}\t"
            f"z={term.output_marker:.9g}\tlog2={term.log2_bound:.6f}"
        )

    print("sampled coefficient-saddle terms")
    saddle_terms = []
    for r in diagnostic_weights(args.k):
        if r <= args.exact_limit or r == args.k:
            continue
        term = saddle_biregular_logterm(
            code=args.code,
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
            memory=args.memory,
            message_weight=r,
        )
        saddle_terms.append(term)
        if r <= 64 or term.log2_bound > -100.0:
            print(
                f"{r}\tq={term.input_probability:.9g}\t"
                f"x={term.input_marker:.9g}\tz={term.output_marker:.9g}\t"
                f"log2={term.log2_bound:.6f}"
            )
    worst = max(exact_terms + saddle_terms, key=lambda term: term.log2_bound)
    print(
        f"sampled worst: r={worst.message_weight}, "
        f"log2_bound={worst.log2_bound:.6f}"
    )
    if args.k % 2 and args.right_degree % 2:
        fourier_l1 = central_fourier_l1_log2(
            left_vertices=args.k, right_degree=args.right_degree
        )
        central = central_dense_logterm(
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
        )
        print(
            f"central Fourier: log2(L1)={fourier_l1:.9f}, "
            f"log2(one-support term)<={central:.6f}"
        )
    if args.right_degree == 5:
        dense = degree_five_dense_scan(
            k=args.k,
            left_degree=args.left_degree,
            cutoff=args.cutoff,
            support_start=max(1, args.k // 5),
            support_limit=args.k - args.k // 5,
        )
        print(
            f"degree-five dense scan: r={dense.support_start}..{dense.support_limit}, "
            f"worst_r={dense.worst_support}, "
            f"worst_log2={dense.worst_log2_term:.6f}, "
            f"summed_log2<={dense.summed_log2_bound:.6f}"
        )
        if args.scan_blocks:
            low_blocks = geometric_blocks(
                args.exact_limit + 1,
                args.k // 5 - 1,
                args.block_relative_width,
            )
            complement_blocks = geometric_blocks(
                1,
                args.k // 5 - 1,
                args.block_relative_width,
            )
            high_blocks = [
                (args.k - complement_hi, args.k - complement_lo)
                for complement_lo, complement_hi in reversed(complement_blocks)
            ]
            block_terms = [
                saddle_biregular_block_logterm(
                    code=args.code,
                    k=args.k,
                    left_degree=args.left_degree,
                    right_degree=args.right_degree,
                    cutoff=args.cutoff,
                    memory=args.memory,
                    support_start=lo,
                    support_limit=hi,
                )
                for lo, hi in low_blocks + high_blocks
            ]
            maximum = max(block_terms, key=lambda term: term.log2_bound)
            logs = np.array(
                [term.log2_bound * math.log(2.0) for term in block_terms]
            )
            total_log2 = float(logsumexp(logs) / math.log(2.0))
            print(
                f"coefficient blocks: count={len(block_terms)}, "
                f"worst={maximum.support_start}..{maximum.support_limit}, "
                f"worst_log2={maximum.log2_bound:.6f}, "
                f"summed_log2={total_log2:.6f}"
            )


if __name__ == "__main__":
    main()
