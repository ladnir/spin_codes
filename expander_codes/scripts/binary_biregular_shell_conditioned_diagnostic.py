#!/usr/bin/env python3
"""Shell-conditioned coefficient diagnostic for binary biregular EC.

For one message support, the right-degree-three expander gives an exact law
for the parity weight in each region.  Conditional on that weight, the
convolution input is a uniform fixed-weight string.  This script bounds its
transfer entrywise with one positive input marker and then averages under the
exact parity-weight law.  It is a floating-point diagnostic, not an interval
certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.special import logsumexp

from binary_biregular_activation_diagnostic import (
    activation_bands,
    exact_logterm,
    log_matrix_power_entries,
)
from binary_biregular_diagnostic import degree_three_parity_shell_counts
from expander_bounds import ec_enumerator_matrix, log_binom, log_weight_mgf
from expander_bounds import normalized


def convolution_activation_matrix(
    *, input_marker: float, output_marker: float,
    activation_marker: float, memory: int,
) -> np.ndarray:
    """Return the convolution transfer marked by input and activation count."""
    matrix = ec_enumerator_matrix(
        input_marker, output_marker, memory, True
    )
    matrix[memory, 0] *= activation_marker
    return matrix


def convolution_duration_matrix(
    *, input_marker: float, output_marker: float,
    duration_marker: float, memory: int,
) -> np.ndarray:
    """Return the convolution transfer marked by inactive zero-input steps."""
    matrix = ec_enumerator_matrix(
        input_marker, output_marker, memory, True
    )
    matrix[memory, memory] *= duration_marker
    return matrix


def convolution_duration_class_matrix(
    *, input_marker: float, output_marker: float, memory: int,
) -> np.ndarray:
    """Mark whether a path uses an inactive zero-input self-loop.

    The first layer contains paths that have not used the self-loop.  That
    transition moves a path to the second layer, which then follows the
    original convolution transfer.  Consequently, the upper-left and
    upper-right blocks of a power enumerate the zero-duration and
    positive-duration path classes without subtraction.
    """
    base = ec_enumerator_matrix(
        input_marker, output_marker, memory, True
    )
    size = memory + 1
    result = np.zeros((2 * size, 2 * size), dtype=np.float64)
    result[:size, :size] = base
    result[memory, memory] = 0.0
    result[memory, size + memory] = 1.0
    result[size:, size:] = base
    return result


def scaled_matrix_power(
    matrix: np.ndarray, exponent: int
) -> tuple[np.ndarray, float]:
    """Return a normalized nonnegative matrix power and its logarithmic scale."""
    result = np.eye(matrix.shape[0], dtype=np.float64)
    result_log_scale = 0.0
    power, power_log_scale = normalized(matrix)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result, local_scale = normalized(result @ power)
            result_log_scale += power_log_scale + local_scale
        remaining >>= 1
        if remaining:
            power, local_scale = normalized(power @ power)
            power_log_scale = 2.0 * power_log_scale + local_scale
    return result, result_log_scale


def first_wait_class_power_logs(
    *, input_marker: float, output_marker: float, memory: int,
    length: int, cut: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return logs for paths with and without an inactive wait before ``cut``."""
    if not 0 <= cut <= length:
        raise ValueError("cut must lie in [0,length]")
    size = memory + 1
    augmented = convolution_duration_class_matrix(
        input_marker=input_marker,
        output_marker=output_marker,
        memory=memory,
    )
    prefix, prefix_scale = scaled_matrix_power(augmented, cut)
    suffix_base = ec_enumerator_matrix(
        input_marker, output_marker, memory, True
    )
    suffix, suffix_scale = scaled_matrix_power(suffix_base, length - cut)
    late, late_local_scale = normalized(prefix[:size, :size] @ suffix)
    early, early_local_scale = normalized(prefix[:size, size:] @ suffix)
    with np.errstate(divide="ignore"):
        late_logs = (
            np.log(late) + prefix_scale + suffix_scale + late_local_scale
        )
        early_logs = (
            np.log(early) + prefix_scale + suffix_scale + early_local_scale
        )
    return early_logs, late_logs


def geometric_duration_bands(
    length: int, relative_width: float
) -> list[tuple[int, int]]:
    """Partition ``0..length`` into a singleton and geometric intervals."""
    if length < 0 or relative_width <= 0.0:
        raise ValueError("invalid duration-band parameters")
    result = [(0, 0)]
    lo = 1
    while lo <= length:
        width = max(1, int(math.ceil(relative_width * lo)))
        hi = min(length, lo + width - 1)
        result.append((lo, hi))
        lo = hi + 1
    return result


def shell_conditioned_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_marker_min: float, log_marker_max: float, marker_steps: int,
) -> tuple[float, np.ndarray, dict[int, np.ndarray]]:
    """Bound one first-moment term after conditioning on parity weight."""
    if right_degree != 3 or k % right_degree:
        raise ValueError("this diagnostic requires compatible right degree three")
    if marker_steps < 1:
        raise ValueError("marker grid must be nonempty")
    region_length = k // right_degree
    denominator, shell_counts = degree_three_parity_shell_counts(
        left_vertices=k, support_size=message_weight
    )
    size = memory + 1
    slice_logs = {
        weight: np.full((size, size), math.inf, dtype=np.float64)
        for weight in shell_counts
    }
    for log_y in np.linspace(log_marker_min, log_marker_max, marker_steps):
        y = math.exp(float(log_y))
        matrix = ec_enumerator_matrix(y, output_marker, memory, True)
        powered = log_matrix_power_entries(matrix, region_length)
        for weight, current in slice_logs.items():
            candidate = (
                powered
                - weight * float(log_y)
                - log_binom(region_length, weight)
            )
            slice_logs[weight] = np.minimum(current, candidate)

    regional_terms = []
    for weight, count in shell_counts.items():
        log_probability = math.log(count) - math.log(denominator)
        regional_terms.append(slice_logs[weight] + log_probability)
    regional_logs = logsumexp(np.stack(regional_terms, axis=0), axis=0)
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    return log_term / math.log(2.0), regional_logs, slice_logs


def activation_refined_shell_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_marker_min: float, log_marker_max: float, marker_steps: int,
    refined_shells: int, activation_band_width: int,
    log_activation_limit: float, activation_steps: int,
    coupled_activation_input: bool = True,
) -> tuple[float, list[int], np.ndarray, dict[str, float]]:
    """Refine dominant parity shells by activation-count bands."""
    simple_bound, _, slice_logs = shell_conditioned_logterm(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        memory=memory,
        message_weight=message_weight,
        output_marker=output_marker,
        log_marker_min=log_marker_min,
        log_marker_max=log_marker_max,
        marker_steps=marker_steps,
    )
    del simple_bound
    denominator, shell_counts = degree_three_parity_shell_counts(
        left_vertices=k, support_size=message_weight
    )
    selected = sorted(shell_counts, reverse=True)[:refined_shells]
    size = memory + 1
    band_logs = {
        weight: np.full(
            (len(activation_bands(weight, activation_band_width)), size, size),
            math.inf,
            dtype=np.float64,
        )
        for weight in selected
    }
    input_grid = np.linspace(log_marker_min, log_marker_max, marker_steps)
    activation_grid = np.linspace(
        -log_activation_limit, log_activation_limit, activation_steps
    )
    region_length = k // right_degree
    for input_coordinate in input_grid:
        for log_v in activation_grid:
            log_y = (
                float(input_coordinate) - float(log_v)
                if coupled_activation_input
                else float(input_coordinate)
            )
            matrix = convolution_activation_matrix(
                input_marker=math.exp(log_y),
                output_marker=output_marker,
                activation_marker=math.exp(float(log_v)),
                memory=memory,
            )
            powered = log_matrix_power_entries(matrix, region_length)
            for weight in selected:
                denominator_log = log_binom(region_length, weight)
                for index, (lo, hi) in enumerate(
                    activation_bands(weight, activation_band_width)
                ):
                    endpoint = max(-lo * log_v, -hi * log_v)
                    candidate = (
                        powered
                        - weight * log_y
                        - denominator_log
                        + math.log(hi - lo + 1)
                        + endpoint
                    )
                    band_logs[weight][index] = np.minimum(
                        band_logs[weight][index], candidate
                    )

    improvement_values = []
    for weight in selected:
        # A positive-weight input cannot leave an inactive start unchanged
        # without first taking an activation transition.
        if weight > 0:
            band_logs[weight][0, memory, :] = -math.inf
        refined = logsumexp(band_logs[weight], axis=0)
        finite = np.isfinite(slice_logs[weight]) & np.isfinite(refined)
        improvements = np.maximum(
            0.0, slice_logs[weight][finite] - refined[finite]
        ) / math.log(2.0)
        improvement_values.extend(improvements.tolist())
        slice_logs[weight] = np.minimum(slice_logs[weight], refined)

    regional_terms = [
        slice_logs[weight] + math.log(count) - math.log(denominator)
        for weight, count in shell_counts.items()
    ]
    regional_logs = logsumexp(np.stack(regional_terms, axis=0), axis=0)
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    positive = np.array(
        [value for value in improvement_values if value > 0.0],
        dtype=np.float64,
    )
    stats = {
        "improved_entries": float(len(positive)),
        "median_improvement_bits": (
            float(np.median(positive)) if len(positive) else 0.0
        ),
        "max_improvement_bits": (
            float(np.max(positive)) if len(positive) else 0.0
        ),
    }
    return log_term / math.log(2.0), selected, regional_logs, stats


def duration_refined_shell_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_marker_min: float, log_marker_max: float, marker_steps: int,
    refined_shells: int, duration_relative_width: float,
    log_duration_limit: float, duration_steps: int,
) -> tuple[float, list[int], list[tuple[int, int]], np.ndarray, dict[str, float]]:
    """Refine dominant parity shells by inactive-duration bands."""
    _, _, slice_logs = shell_conditioned_logterm(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        memory=memory,
        message_weight=message_weight,
        output_marker=output_marker,
        log_marker_min=log_marker_min,
        log_marker_max=log_marker_max,
        marker_steps=marker_steps,
    )
    denominator, shell_counts = degree_three_parity_shell_counts(
        left_vertices=k, support_size=message_weight
    )
    selected = sorted(shell_counts, reverse=True)[:refined_shells]
    region_length = k // right_degree
    bands = geometric_duration_bands(region_length, duration_relative_width)
    size = memory + 1
    band_logs = {
        weight: np.full(
            (len(bands), size, size), math.inf, dtype=np.float64
        )
        for weight in selected
    }
    input_grid = np.linspace(log_marker_min, log_marker_max, marker_steps)
    duration_grid = np.linspace(
        -log_duration_limit, log_duration_limit, duration_steps
    )
    for log_y in input_grid:
        y = math.exp(float(log_y))
        for log_u in duration_grid:
            matrix = convolution_duration_matrix(
                input_marker=y,
                output_marker=output_marker,
                duration_marker=math.exp(float(log_u)),
                memory=memory,
            )
            powered = log_matrix_power_entries(matrix, region_length)
            for weight in selected:
                base = (
                    powered
                    - weight * float(log_y)
                    - log_binom(region_length, weight)
                )
                for index, (lo, hi) in enumerate(bands):
                    endpoint = max(-lo * log_u, -hi * log_u)
                    candidate = base + math.log(hi - lo + 1) + endpoint
                    band_logs[weight][index] = np.minimum(
                        band_logs[weight][index], candidate
                    )

    improvement_values = []
    for weight in selected:
        refined = logsumexp(band_logs[weight], axis=0)
        finite = np.isfinite(slice_logs[weight]) & np.isfinite(refined)
        improvements = np.maximum(
            0.0, slice_logs[weight][finite] - refined[finite]
        ) / math.log(2.0)
        improvement_values.extend(improvements.tolist())
        slice_logs[weight] = np.minimum(slice_logs[weight], refined)

    regional_terms = [
        slice_logs[weight] + math.log(count) - math.log(denominator)
        for weight, count in shell_counts.items()
    ]
    regional_logs = logsumexp(np.stack(regional_terms, axis=0), axis=0)
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    positive = np.array(
        [value for value in improvement_values if value > 0.0],
        dtype=np.float64,
    )
    stats = {
        "improved_entries": float(len(positive)),
        "median_improvement_bits": (
            float(np.median(positive)) if len(positive) else 0.0
        ),
        "max_improvement_bits": (
            float(np.max(positive)) if len(positive) else 0.0
        ),
    }
    return log_term / math.log(2.0), selected, bands, regional_logs, stats


def duration_class_refined_shell_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_marker_min: float, log_marker_max: float, marker_steps: int,
    refined_shells: int,
) -> tuple[float, list[int], np.ndarray, dict[str, float]]:
    """Refine dominant shells by exact zero/positive inactive duration."""
    _, _, slice_logs = shell_conditioned_logterm(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        memory=memory,
        message_weight=message_weight,
        output_marker=output_marker,
        log_marker_min=log_marker_min,
        log_marker_max=log_marker_max,
        marker_steps=marker_steps,
    )
    denominator, shell_counts = degree_three_parity_shell_counts(
        left_vertices=k, support_size=message_weight
    )
    selected = sorted(shell_counts, reverse=True)[:refined_shells]
    region_length = k // right_degree
    size = memory + 1
    class_logs = {
        weight: np.full((2, size, size), math.inf, dtype=np.float64)
        for weight in selected
    }
    for log_y in np.linspace(log_marker_min, log_marker_max, marker_steps):
        matrix = convolution_duration_class_matrix(
            input_marker=math.exp(float(log_y)),
            output_marker=output_marker,
            memory=memory,
        )
        powered = log_matrix_power_entries(matrix, region_length)
        zero_duration = powered[:size, :size]
        positive_duration = powered[:size, size:]
        for weight in selected:
            adjustment = (
                -weight * float(log_y)
                - log_binom(region_length, weight)
            )
            class_logs[weight][0] = np.minimum(
                class_logs[weight][0], zero_duration + adjustment
            )
            class_logs[weight][1] = np.minimum(
                class_logs[weight][1], positive_duration + adjustment
            )

    improvement_values = []
    class_improvement_values = [[], []]
    for weight in selected:
        refined = logsumexp(class_logs[weight], axis=0)
        finite = np.isfinite(slice_logs[weight]) & np.isfinite(refined)
        improvements = np.maximum(
            0.0, slice_logs[weight][finite] - refined[finite]
        ) / math.log(2.0)
        improvement_values.extend(improvements.tolist())
        for class_index in range(2):
            class_finite = (
                np.isfinite(slice_logs[weight])
                & np.isfinite(class_logs[weight][class_index])
            )
            values = np.maximum(
                0.0,
                slice_logs[weight][class_finite]
                - class_logs[weight][class_index][class_finite],
            ) / math.log(2.0)
            class_improvement_values[class_index].extend(values.tolist())
        slice_logs[weight] = np.minimum(slice_logs[weight], refined)

    regional_terms = [
        slice_logs[weight] + math.log(count) - math.log(denominator)
        for weight, count in shell_counts.items()
    ]
    regional_logs = logsumexp(np.stack(regional_terms, axis=0), axis=0)
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    positive = np.array(
        [value for value in improvement_values if value > 0.0],
        dtype=np.float64,
    )
    stats = {
        "improved_entries": float(len(positive)),
        "median_improvement_bits": (
            float(np.median(positive)) if len(positive) else 0.0
        ),
        "max_improvement_bits": (
            float(np.max(positive)) if len(positive) else 0.0
        ),
    }
    for index, name in enumerate(("zero_duration", "positive_duration")):
        values = np.array(
            [value for value in class_improvement_values[index] if value > 0.0],
            dtype=np.float64,
        )
        stats[f"{name}_max_improvement_bits"] = (
            float(np.max(values)) if len(values) else 0.0
        )
    return log_term / math.log(2.0), selected, regional_logs, stats


def first_wait_refined_shell_logterm(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: float,
    log_marker_min: float, log_marker_max: float, marker_steps: int,
    refined_shells: int, cut_fraction: float,
) -> tuple[float, list[int], int, np.ndarray, dict[str, float]]:
    """Refine dominant shells by whether the first inactive wait is early."""
    if not 0.0 <= cut_fraction <= 1.0:
        raise ValueError("cut_fraction must lie in [0,1]")
    _, _, slice_logs = shell_conditioned_logterm(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        memory=memory,
        message_weight=message_weight,
        output_marker=output_marker,
        log_marker_min=log_marker_min,
        log_marker_max=log_marker_max,
        marker_steps=marker_steps,
    )
    denominator, shell_counts = degree_three_parity_shell_counts(
        left_vertices=k, support_size=message_weight
    )
    selected = sorted(shell_counts, reverse=True)[:refined_shells]
    region_length = k // right_degree
    cut = int(round(cut_fraction * region_length))
    size = memory + 1
    class_logs = {
        weight: np.full((2, size, size), math.inf, dtype=np.float64)
        for weight in selected
    }
    for log_y in np.linspace(log_marker_min, log_marker_max, marker_steps):
        early, late = first_wait_class_power_logs(
            input_marker=math.exp(float(log_y)),
            output_marker=output_marker,
            memory=memory,
            length=region_length,
            cut=cut,
        )
        for weight in selected:
            adjustment = (
                -weight * float(log_y)
                - log_binom(region_length, weight)
            )
            class_logs[weight][0] = np.minimum(
                class_logs[weight][0], early + adjustment
            )
            class_logs[weight][1] = np.minimum(
                class_logs[weight][1], late + adjustment
            )

    improvement_values = []
    for weight in selected:
        refined = logsumexp(class_logs[weight], axis=0)
        finite = np.isfinite(slice_logs[weight]) & np.isfinite(refined)
        improvements = np.maximum(
            0.0, slice_logs[weight][finite] - refined[finite]
        ) / math.log(2.0)
        improvement_values.extend(improvements.tolist())
        slice_logs[weight] = np.minimum(slice_logs[weight], refined)

    regional_terms = [
        slice_logs[weight] + math.log(count) - math.log(denominator)
        for weight, count in shell_counts.items()
    ]
    regional_logs = logsumexp(np.stack(regional_terms, axis=0), axis=0)
    scale = float(np.max(regional_logs))
    regional = np.exp(regional_logs - scale)
    log_term = (
        log_binom(k, message_weight)
        + left_degree * scale
        + log_weight_mgf(regional, left_degree, memory)
        - cutoff * math.log(output_marker)
    )
    positive = np.array(
        [value for value in improvement_values if value > 0.0],
        dtype=np.float64,
    )
    stats = {
        "improved_entries": float(len(positive)),
        "median_improvement_bits": (
            float(np.median(positive)) if len(positive) else 0.0
        ),
        "max_improvement_bits": (
            float(np.max(positive)) if len(positive) else 0.0
        ),
    }
    return log_term / math.log(2.0), selected, cut, regional_logs, stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--left-degree", type=int, default=6)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=80)
    parser.add_argument("--message-weight", type=int, default=256)
    parser.add_argument("--output-marker", type=float, default=0.992)
    parser.add_argument("--log-marker-min", type=float, default=-12.0)
    parser.add_argument("--log-marker-max", type=float, default=-2.0)
    parser.add_argument("--marker-steps", type=int, default=101)
    parser.add_argument("--refined-shells", type=int, default=0)
    parser.add_argument("--activation-band-width", type=int, default=8)
    parser.add_argument("--log-activation-limit", type=float, default=8.0)
    parser.add_argument("--activation-steps", type=int, default=33)
    parser.add_argument(
        "--independent-activation-input",
        action="store_true",
        help="vary the activation marker while keeping the input marker fixed",
    )
    parser.add_argument("--compare-exact-slice", action="store_true")
    parser.add_argument("--duration-refined-shells", type=int, default=0)
    parser.add_argument("--duration-relative-width", type=float, default=0.25)
    parser.add_argument("--log-duration-limit", type=float, default=8.0)
    parser.add_argument("--duration-steps", type=int, default=33)
    parser.add_argument("--duration-class-refined-shells", type=int, default=0)
    parser.add_argument("--first-wait-refined-shells", type=int, default=0)
    parser.add_argument("--first-wait-cut", type=float, default=0.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    common = dict(
        k=args.k,
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        cutoff=args.cutoff,
        memory=args.memory,
        message_weight=args.message_weight,
        output_marker=args.output_marker,
    )
    exact = exact_logterm(**common)
    if args.first_wait_refined_shells:
        bound, selected, cut, regional_logs, refinement_stats = (
            first_wait_refined_shell_logterm(
                **common,
                log_marker_min=args.log_marker_min,
                log_marker_max=args.log_marker_max,
                marker_steps=args.marker_steps,
                refined_shells=args.first_wait_refined_shells,
                cut_fraction=args.first_wait_cut,
            )
        )
        shell_count = len(
            degree_three_parity_shell_counts(
                left_vertices=args.k, support_size=args.message_weight
            )[1]
        )
        print(f"first_wait_refined_shells={selected}")
        print(f"first_wait_cut={cut}")
        print(f"refinement_stats={refinement_stats}")
    elif args.duration_class_refined_shells:
        bound, selected, regional_logs, refinement_stats = (
            duration_class_refined_shell_logterm(
                **common,
                log_marker_min=args.log_marker_min,
                log_marker_max=args.log_marker_max,
                marker_steps=args.marker_steps,
                refined_shells=args.duration_class_refined_shells,
            )
        )
        shell_count = len(
            degree_three_parity_shell_counts(
                left_vertices=args.k, support_size=args.message_weight
            )[1]
        )
        print(f"duration_class_refined_shells={selected}")
        print(f"refinement_stats={refinement_stats}")
    elif args.duration_refined_shells:
        bound, selected, duration_bands, regional_logs, refinement_stats = (
            duration_refined_shell_logterm(
                **common,
                log_marker_min=args.log_marker_min,
                log_marker_max=args.log_marker_max,
                marker_steps=args.marker_steps,
                refined_shells=args.duration_refined_shells,
                duration_relative_width=args.duration_relative_width,
                log_duration_limit=args.log_duration_limit,
                duration_steps=args.duration_steps,
            )
        )
        shell_count = len(
            degree_three_parity_shell_counts(
                left_vertices=args.k, support_size=args.message_weight
            )[1]
        )
        print(f"duration_refined_shells={selected}")
        print(f"duration_bands={len(duration_bands)}")
        print(f"refinement_stats={refinement_stats}")
    elif args.refined_shells:
        bound, selected, regional_logs, refinement_stats = (
            activation_refined_shell_logterm(
            **common,
            log_marker_min=args.log_marker_min,
            log_marker_max=args.log_marker_max,
            marker_steps=args.marker_steps,
            refined_shells=args.refined_shells,
            activation_band_width=args.activation_band_width,
            log_activation_limit=args.log_activation_limit,
            activation_steps=args.activation_steps,
            coupled_activation_input=not args.independent_activation_input,
            )
        )
        shell_count = len(
            degree_three_parity_shell_counts(
                left_vertices=args.k, support_size=args.message_weight
            )[1]
        )
        print(f"activation_refined_shells={selected}")
        print(f"refinement_stats={refinement_stats}")
    else:
        bound, regional_logs, slice_logs = shell_conditioned_logterm(
            **common,
            log_marker_min=args.log_marker_min,
            log_marker_max=args.log_marker_max,
            marker_steps=args.marker_steps,
        )
        shell_count = len(slice_logs)
    print(f"exact_log2_term={exact:.12f}")
    print(f"shell_conditioned_log2_bound={bound:.12f}")
    print(f"gap_bits={bound - exact:.12f}")
    print(f"parity_shells={shell_count}")
    print(f"regional_log2_max={float(np.max(regional_logs)) / math.log(2):.12f}")
    if args.compare_exact_slice:
        from binary_biregular_hybrid_diagnostic import (
            scaled_uniform_slice_transfers,
        )

        largest_shell = max(
            degree_three_parity_shell_counts(
                left_vertices=args.k, support_size=args.message_weight
            )[1]
        )
        slices = scaled_uniform_slice_transfers(
            length=args.k // args.right_degree,
            max_weight=largest_shell,
            output_marker=args.output_marker,
            memory=args.memory,
        )
        exact_matrix, exact_scale = slices[largest_shell]
        with np.errstate(divide="ignore"):
            exact_slice_logs = np.log(exact_matrix) + exact_scale
        if any((
            args.refined_shells,
            args.duration_refined_shells,
            args.duration_class_refined_shells,
            args.first_wait_refined_shells,
        )):
            simple_slice = shell_conditioned_logterm(
                **common,
                log_marker_min=args.log_marker_min,
                log_marker_max=args.log_marker_max,
                marker_steps=args.marker_steps,
            )[2][largest_shell]
        else:
            simple_slice = slice_logs[largest_shell]
        finite = np.isfinite(exact_slice_logs)
        gaps = (simple_slice[finite] - exact_slice_logs[finite]) / math.log(2)
        print(f"exact_slice={largest_shell}")
        print(f"exact_slice_nonzero_entries={int(np.count_nonzero(finite))}")
        print(f"slice_gap_min_bits={float(np.min(gaps)):.12f}")
        print(f"slice_gap_median_bits={float(np.median(gaps)):.12f}")
        print(f"slice_gap_max_bits={float(np.max(gaps)):.12f}")


if __name__ == "__main__":
    main()
