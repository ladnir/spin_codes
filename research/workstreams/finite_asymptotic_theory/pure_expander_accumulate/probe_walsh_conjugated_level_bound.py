#!/usr/bin/env python3
"""Probe a Walsh-conjugated level bound for one-stage sparse EA.

The original covariance operator acts on Fourier subsets.  Conjugating it
by the normalized Walsh transform moves the calculation to the primal
pre-accumulator shell indicator.  That indicator has a sparse and exactly
known level profile.  This program tests whether a level-norm relaxation in
that basis retains enough cancellation to be useful.

The identities are exact.  The implementation uses ordinary binary64 and
is therefore a diagnostic, not an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np

from analyze_dual_walk_and_accumulator_energy import accumulator_joint_counts
from probe_radial_block_bound import marginal_mean, radial_walk_table
from verify_pair_kernel_small import accumulator
from verify_pair_kernel_small import krawtchouk
from verify_subset_covariance_blocks_small import (
    covariance_kernel as exact_covariance_kernel,
    region_counts,
    shell_coefficient,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "walsh_conjugated_level_bound_probe.json"
Kernel = Callable[[int, int, int], float]


def invariant_blocks(
    output_bits: int, kernel: Kernel, *, exact_integer_terms: bool = False
) -> list[np.ndarray]:
    """Return normalized two-row Specht blocks of an invariant kernel."""
    blocks: list[np.ndarray] = []
    for sector in range(output_bits // 2 + 1):
        outside = output_bits - 2 * sector
        levels = range(sector, output_bits - sector + 1)
        size = output_bits - 2 * sector + 1
        matrix = np.zeros((size, size), dtype=float)
        for row, first_weight in enumerate(levels):
            first_outside = first_weight - sector
            first_count = math.comb(outside, first_outside)
            for column, second_weight in enumerate(levels):
                second_outside = second_weight - sector
                second_count = math.comb(outside, second_outside)
                terms = []
                for same_pairs in range(sector + 1):
                    pair_factor = (
                        (-1) ** (sector - same_pairs)
                        * math.comb(sector, same_pairs)
                    )
                    low = max(0, second_outside - (outside - first_outside))
                    high = min(first_outside, second_outside)
                    for outside_intersection in range(low, high + 1):
                        ways = math.comb(
                            first_outside, outside_intersection
                        ) * math.comb(
                            outside - first_outside,
                            second_outside - outside_intersection,
                        )
                        intersection = same_pairs + outside_intersection
                        terms.append(
                            pair_factor
                            * ways
                            * kernel(first_weight, second_weight, intersection)
                        )
                total = sum(terms) if exact_integer_terms else math.fsum(terms)
                matrix[row, column] = math.sqrt(
                    first_count / second_count
                ) * total
        residual = float(np.max(np.abs(matrix - matrix.T)))
        scale = max(1.0, float(np.max(np.abs(matrix))))
        if residual > 2e-11 * scale:
            raise ArithmeticError(
                f"sector {sector} symmetry residual {residual} is too large"
            )
        blocks.append((matrix + matrix.T) / 2)
    return blocks


def covariance_blocks(
    message_bits: int, output_bits: int, right_degree: int
) -> tuple[list[np.ndarray], int]:
    walk = radial_walk_table(message_bits, right_degree, output_bits)
    multiplicities = np.asarray(
        [math.comb(message_bits, weight) for weight in range(message_bits + 1)],
        dtype=float,
    )

    @functools.lru_cache(maxsize=None)
    def covariance(first_weight: int, second_weight: int, intersection: int) -> float:
        first_only = first_weight - intersection
        second_only = second_weight - intersection
        joint = float(
            np.dot(
                multiplicities * walk[first_only] * walk[second_only],
                walk[intersection],
            )
        )
        return joint - walk[first_weight, 0] * walk[second_weight, 0]

    blocks = invariant_blocks(output_bits, covariance)
    return blocks, covariance.cache_info().currsize


def walsh_blocks(output_bits: int) -> tuple[list[np.ndarray], float]:
    blocks = invariant_blocks(
        output_bits,
        lambda _first, _second, intersection: (-1) ** intersection,
        exact_integer_terms=True,
    )
    target = math.ldexp(1.0, output_bits)
    maximum_residual = 0.0
    for block in blocks:
        residual = np.max(np.abs(block @ block.T - target * np.eye(len(block))))
        maximum_residual = max(maximum_residual, float(residual / target))
    return blocks, maximum_residual


def primal_level_counts(output_bits: int, shell_weight: int) -> list[int]:
    """Count z by wt(z) subject to wt(A z)=shell_weight."""
    joint = accumulator_joint_counts(output_bits)
    counts = [0] * (output_bits + 1)
    for derivative_weight, count in joint[shell_weight].items():
        counts[derivative_weight] = count
    if sum(counts) != math.comb(output_bits, shell_weight):
        raise AssertionError("primal shell level counts have the wrong mass")
    return counts


def conditional_determinant_maxima(
    message_bits: int, right_degree: int
) -> tuple[list[float], list[float]]:
    """Return beta(h) and max_y |det(P_{x,y})| for every wt(x)=h."""
    denominator = math.comb(message_bits, right_degree)
    biases = [
        krawtchouk(message_bits, right_degree, weight) / denominator
        for weight in range(message_bits + 1)
    ]
    maxima = [0.0] * (message_bits + 1)
    for first_weight in range(message_bits + 1):
        first_bias = biases[first_weight]
        maximum = 0.0
        for overlap in range(first_weight + 1):
            for second_only in range(message_bits - first_weight + 1):
                second_weight = overlap + second_only
                difference_weight = first_weight - overlap + second_only
                determinant = abs(
                    (
                        biases[difference_weight]
                        - first_bias * biases[second_weight]
                    )
                    / 4
                )
                maximum = max(maximum, determinant)
        maxima[first_weight] = maximum
    return biases, maxima


def higher_sector_diagonal_envelope(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    maximum_level: int,
    minimum_sector: int = 2,
) -> tuple[np.ndarray, dict[str, object]]:
    """Bound max_j D^(j)_{p,p} for j >= minimum_sector.

    The bound uses the exact conditional determinant second moment and the
    largest determinant available at each first-message weight.
    """
    biases, determinant_maxima = conditional_determinant_maxima(
        message_bits, right_degree
    )
    row_count = math.comb(message_bits, right_degree)
    determinant_second_scale = math.ldexp(1.0, message_bits) / (16 * row_count)
    message_multiplicities = [
        math.comb(message_bits, weight) for weight in range(message_bits + 1)
    ]
    envelope = np.zeros(maximum_level + 1, dtype=float)
    maximizing_sector = np.zeros(maximum_level + 1, dtype=int)
    for level in range(minimum_sector, maximum_level + 1):
        maximum = 0.0
        argmax = minimum_sector
        for sector in range(minimum_sector, level + 1):
            degree = output_bits - 2 * sector
            symmetric_weight = level - sector
            terms = []
            for weight in range(message_bits + 1):
                bias = biases[weight]
                maximum_determinant = determinant_maxima[weight]
                if maximum_determinant == 0 and sector > 2:
                    continue
                one_probability = (1 - bias) / 2
                marginal = (
                    one_probability**symmetric_weight
                    * (1 - one_probability) ** (degree - symmetric_weight)
                )
                terms.append(
                    message_multiplicities[weight]
                    * (1 - bias * bias)
                    * marginal
                    * maximum_determinant ** (sector - 2)
                )
            value = determinant_second_scale * math.fsum(terms)
            if value > maximum:
                maximum = value
                argmax = sector
        envelope[level] = maximum
        maximizing_sector[level] = argmax
    return envelope, {
        "minimum_sector": minimum_sector,
        "maximum_conditional_determinant": max(determinant_maxima),
        "maximizing_sectors": {
            str(level): int(maximizing_sector[level])
            for level in range(minimum_sector, maximum_level + 1)
        },
    }


def known_exact_ratio(
    message_bits: int, output_bits: int, right_degree: int, shell_weight: int
) -> float | None:
    pattern = f"one_stage_variance_K{message_bits}_B{output_bits}_r{right_degree}_probe.json"
    path = HERE / pattern
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = payload["result"]
    row = result["shells"][shell_weight]
    return row["variance_to_mean"]


def verify_small_conjugation(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> dict[str, float] | None:
    """Check the primal/dual Walsh identity by full enumeration of subsets."""
    if output_bits > 10:
        return None
    _returns, exact_kernel = exact_covariance_kernel(
        message_bits, output_bits, right_degree
    )
    size = 1 << output_bits
    covariance = np.asarray(
        [
            [
                float(exact_kernel[region_counts(first, second, output_bits)])
                for second in range(size)
            ]
            for first in range(size)
        ],
        dtype=float,
    )
    walsh = np.asarray(
        [
            [float((-1) ** ((first & second).bit_count())) for second in range(size)]
            for first in range(size)
        ],
        dtype=float,
    )
    primal = np.asarray(
        [
            float(accumulator(word, output_bits).bit_count() == shell_weight)
            for word in range(size)
        ]
    )
    dual = np.asarray(
        [
            float(shell_coefficient(output_bits, shell_weight, subset))
            for subset in range(size)
        ]
    )
    coefficient_residual = float(np.max(np.abs(walsh @ primal - dual)))
    dual_quadratic = float(dual @ covariance @ dual)
    transformed = math.ldexp(1.0, -output_bits) * walsh @ covariance @ walsh
    primal_quadratic = float(primal @ transformed @ primal)
    scaled_residual = abs(dual_quadratic - math.ldexp(primal_quadratic, output_bits))
    relative_residual = scaled_residual / max(1.0, abs(dual_quadratic))
    if coefficient_residual != 0 or relative_residual > 3e-12:
        raise AssertionError("small Walsh conjugation check failed")
    return {
        "coefficient_maximum_residual": coefficient_residual,
        "quadratic_relative_residual": relative_residual,
    }


def evaluate(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    shell_weight: int,
) -> dict[str, object]:
    covariance, cached_entries = covariance_blocks(
        message_bits, output_bits, right_degree
    )
    walsh, walsh_residual = walsh_blocks(output_bits)
    normalization = math.ldexp(1.0, -output_bits)

    transformed = [
        normalization * transform @ block @ transform.T
        for transform, block in zip(walsh, covariance, strict=True)
    ]
    transform_symmetry_residual = max(
        float(np.max(np.abs(block - block.T))) for block in transformed
    )

    level_counts = primal_level_counts(output_bits, shell_weight)
    level_norms = np.sqrt(np.asarray(level_counts, dtype=float))
    entry_norm = np.zeros((output_bits + 1, output_bits + 1), dtype=float)
    diagonal_norm = np.zeros(output_bits + 1, dtype=float)
    diagonal_sector = np.zeros(output_bits + 1, dtype=int)
    low_sector_diagonal_norm = np.zeros(output_bits + 1, dtype=float)
    level_trace = np.zeros(output_bits + 1, dtype=float)
    for sector, block in enumerate(transformed):
        levels = slice(sector, output_bits - sector + 1)
        entry_norm[levels, levels] = np.maximum(
            entry_norm[levels, levels], np.abs(block)
        )
        block_diagonal = np.diag(block)
        if np.min(block_diagonal) < -1e-12 * max(1.0, np.max(np.abs(block))):
            raise ArithmeticError("transformed covariance block has negative diagonal")
        candidate = np.maximum(block_diagonal, 0)
        current = diagonal_norm[levels]
        improved = candidate > current
        current[improved] = candidate[improved]
        diagonal_sector[levels][improved] = sector
        if sector <= 2:
            low_sector_diagonal_norm[levels] = np.maximum(
                low_sector_diagonal_norm[levels], candidate
            )
        dimension = math.comb(output_bits, sector) - (
            math.comb(output_bits, sector - 1) if sector else 0
        )
        level_trace[levels] += dimension * candidate

    quadratic_bound = float(level_norms @ entry_norm @ level_norms)
    diagonal_quadratic_bound = float(
        np.dot(level_norms, np.sqrt(diagonal_norm)) ** 2
    )
    trace_envelope = low_sector_diagonal_norm.copy()
    dimension_three = math.comb(output_bits, 3) - math.comb(output_bits, 2)
    trace_envelope[3 : output_bits - 2] = np.maximum(
        trace_envelope[3 : output_bits - 2],
        level_trace[3 : output_bits - 2] / dimension_three,
    )
    trace_quadratic_bound = float(
        np.dot(level_norms, np.sqrt(trace_envelope)) ** 2
    )
    mean = marginal_mean(message_bits, output_bits, right_degree, shell_weight)
    variance_scale = math.ldexp(1.0, 2 * message_bits - output_bits)
    ratio_bound = variance_scale * quadratic_bound / mean
    diagonal_ratio_bound = variance_scale * diagonal_quadratic_bound / mean
    trace_ratio_bound = variance_scale * trace_quadratic_bound / mean
    exact_ratio = known_exact_ratio(
        message_bits, output_bits, right_degree, shell_weight
    )
    small_check = verify_small_conjugation(
        message_bits, output_bits, right_degree, shell_weight
    )
    active_levels = [index for index, count in enumerate(level_counts) if count]
    higher_envelope, higher_metadata = higher_sector_diagonal_envelope(
        message_bits,
        output_bits,
        right_degree,
        max(active_levels),
    )
    exact_low_sector_envelope = np.zeros(output_bits + 1, dtype=float)
    for sector in range(min(2, len(transformed))):
        block = transformed[sector]
        levels = slice(sector, output_bits - sector + 1)
        exact_low_sector_envelope[levels] = np.maximum(
            exact_low_sector_envelope[levels], np.maximum(np.diag(block), 0)
        )
    determinant_envelope = exact_low_sector_envelope.copy()
    determinant_envelope[: len(higher_envelope)] = np.maximum(
        determinant_envelope[: len(higher_envelope)], higher_envelope
    )
    determinant_quadratic_bound = float(
        np.dot(level_norms, np.sqrt(determinant_envelope)) ** 2
    )
    determinant_ratio_bound = (
        variance_scale * determinant_quadratic_bound / mean
    )

    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "shell_weight": shell_weight,
        "shell_mean": mean,
        "active_primal_level_minimum": min(active_levels),
        "active_primal_level_maximum": max(active_levels),
        "active_primal_level_count": len(active_levels),
        "psd_diagonal_maximizing_sectors": {
            str(level): int(diagonal_sector[level]) for level in active_levels
        },
        "walsh_conjugated_level_variance_to_mean_upper": ratio_bound,
        "psd_diagonal_variance_to_mean_upper": diagonal_ratio_bound,
        "sector_0_2_plus_trace_variance_to_mean_upper": trace_ratio_bound,
        "sector_0_1_plus_determinant_variance_to_mean_upper": (
            determinant_ratio_bound
        ),
        "higher_sector_determinant_envelope": higher_metadata,
        "known_exact_variance_to_mean": exact_ratio,
        "loss_over_exact": None if exact_ratio is None else ratio_bound / exact_ratio,
        "cached_covariance_entries": cached_entries,
        "normalized_walsh_involution_residual": walsh_residual,
        "transformed_block_symmetry_residual": transform_symmetry_residual,
        "small_full_walsh_check": small_check,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=24)
    parser.add_argument("--output-bits", type=int, default=48)
    parser.add_argument("--right-degree", type=int, default=5)
    parser.add_argument("--shell-weight", type=int, default=31)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits > 128:
        parser.error("this binary64 gate is intentionally limited to output_bits <= 128")
    if args.output_bits != 2 * args.message_bits:
        parser.error("the present probe expects output_bits = 2 * message_bits")
    if not 0 <= args.shell_weight <= args.output_bits:
        parser.error("invalid shell weight")

    result = evaluate(
        args.message_bits,
        args.output_bits,
        args.right_degree,
        args.shell_weight,
    )
    payload = {
        "schema": "pure-ea-walsh-conjugated-level-bound-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "result": result,
        "scope": [
            "The covariance, Specht-block, Walsh-conjugation, and primal-level formulas are exact identities.",
            "The displayed upper bound is valid in exact arithmetic.",
            "This implementation uses nondirected binary64 and is not a rigorous certificate.",
            "The receipt does not extrapolate to length 512.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
