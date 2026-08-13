#!/usr/bin/env python3
"""Diagnostic g=4 inner bound sharing systematic state columns globally.

The committed split-cap relaxation applies the right-half marginal
``C(64,r)`` independently in every emitted-weight row.  In the true
systematic map those columns are shared: across all 64-bit systematic inputs,
exactly ``C(64,r)`` have parity/state weight ``r``.  This probe builds a
second valid upper operator using that global marginal (while dropping the
row-specific split caps), then takes the componentwise minimum with the
existing row-specific shared-drive operator.  The minimum still dominates
the exact image and can be strictly sharper than either relaxation alone.

Binary64 evaluation is diagnostic; outward arithmetic is required for a
certificate.
"""

from __future__ import annotations

import argparse
import bisect
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from packet_group_drive_stratified import block_histograms, point_caps
from packet_group_outer_profile import D, normalization_log2
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS
from probe_packet_group_pattern_exact_transfer import (
    pattern_emission_counts,
    polynomial_for_pattern,
)
from packet_group_drive_stratified import compatible_pair_polynomials


BITS = 64
GROUP_BITS = 4
CLASSES = GROUP_BITS + 1


class GlobalColumnKernel:
    """Point-cap majorization with one shared C(64,r) state marginal."""

    def __init__(self, histograms, caps, pole: float):
        self.prefix_units: list[tuple[float, ...] | None] = [None] * (BITS + 1)
        self.prefix_mass: list[np.ndarray | None] = [None] * (BITS + 1)
        self.source_scores: list[np.ndarray | None] = [None] * (BITS + 1)
        for incoming in range(BITS + 1):
            denominator = float(math.comb(BITS, incoming))
            sources = []
            for drive in range(BITS + 1):
                cap = float(caps[incoming, drive])
                if cap <= 0.0:
                    if np.any(histograms[incoming, drive]):
                        raise RuntimeError("global-column transfer: positive mass has zero cap")
                    continue
                for emitted in range(BITS + 1):
                    mass = float(histograms[incoming, drive, emitted]) / denominator
                    if mass:
                        sources.append((cap * pole**emitted, mass / cap))
            sources.sort(reverse=True)
            if not sources:
                continue
            units = [float(row[1]) for row in sources]
            scores = np.asarray([row[0] for row in sources], dtype=np.float64)
            self.prefix_units[incoming] = tuple(itertools.accumulate(units))
            self.prefix_mass[incoming] = np.cumsum(np.asarray(units) * scores)
            self.source_scores[incoming] = scores
        self.column_capacities = tuple(math.comb(BITS, state) for state in range(BITS + 1))

    def apply(self, values: np.ndarray) -> np.ndarray:
        result = np.zeros(BITS + 1, dtype=np.float64)
        order = np.argsort(values)[::-1]
        boundaries = tuple(
            itertools.accumulate(self.column_capacities[state] for state in order)
        )
        ordered_values = values[order]
        for incoming in range(BITS + 1):
            units = self.prefix_units[incoming]
            mass = self.prefix_mass[incoming]
            scores = self.source_scores[incoming]
            if units is None or mass is None or scores is None:
                continue
            total_units = units[-1]
            total_mass = float(mass[-1])
            filled = []
            for boundary in boundaries:
                clipped = min(float(boundary), total_units)
                index = min(bisect.bisect_left(units, clipped), len(units) - 1)
                prior_units = 0.0 if index == 0 else units[index - 1]
                prior_mass = 0.0 if index == 0 else float(mass[index - 1])
                filled.append(
                    min(total_mass, prior_mass + (clipped - prior_units) * float(scores[index]))
                )
            allocations = np.diff(np.asarray([0.0, *filled], dtype=np.float64))
            result[incoming] = float(allocations @ ordered_values)
        return result


class MinimumUpperKernel:
    def __init__(self, *kernels):
        self.kernels = kernels

    def apply(self, values: np.ndarray) -> np.ndarray:
        images = [kernel.apply(values) for kernel in self.kernels]
        return np.minimum.reduce(images)


def parse_vector(text: str, expected: int, cast=float):
    values = [cast(value) for value in text.split(",")]
    if len(values) != expected:
        raise ValueError(f"expected {expected} comma-separated values")
    return values


def row_specific_column_usage(kernel, values: np.ndarray, incoming: int) -> dict:
    """Recover the diagnostic sink-unit usage of the row-wise greedy."""

    order = np.argsort(values)[::-1]
    usage = np.zeros(BITS + 1, dtype=np.float64)
    by_emitted = np.zeros((BITS + 1, BITS + 1), dtype=np.float64)
    for emitted in range(BITS + 1):
        sources = kernel.sources[incoming][emitted]
        if not sources:
            continue
        capacities = [float(kernel.split_caps[emitted][state]) for state in order]
        source_index = 0
        sink_index = 0
        source_remaining = sources[0][1]
        sink_remaining = capacities[0]
        while source_index < len(sources) and sink_index < len(capacities):
            while sink_index < len(capacities) and sink_remaining <= 0.0:
                sink_index += 1
                if sink_index < len(capacities):
                    sink_remaining = capacities[sink_index]
            if sink_index == len(capacities):
                break
            take = min(source_remaining, sink_remaining)
            state = int(order[sink_index])
            usage[state] += take
            by_emitted[emitted, state] += take
            if source_remaining <= sink_remaining:
                sink_remaining -= source_remaining
                source_index += 1
                if source_index < len(sources):
                    source_remaining = sources[source_index][1]
            else:
                source_remaining -= sink_remaining
                sink_index += 1
                if sink_index < len(capacities):
                    sink_remaining = capacities[sink_index]
    capacities = np.asarray([math.comb(BITS, state) for state in range(BITS + 1)])
    ratios = np.divide(usage, capacities, out=np.zeros_like(usage), where=capacities > 0)
    worst = np.argsort(ratios)[::-1][:10]
    return {
        "incoming_state": incoming,
        "worst_global_column_reuse": [
            {
                "next_state": int(state),
                "used_units": float(usage[state]),
                "global_capacity": int(capacities[state]),
                "reuse_ratio": float(ratios[state]),
                "emitted_rows_used": int(np.count_nonzero(by_emitted[:, state])),
            }
            for state in worst
            if usage[state] > 0.0
        ],
    }


def row_specific_emitted_contributions(kernel, values: np.ndarray, incoming: int) -> list[dict]:
    order = np.argsort(values)[::-1]
    ordered_values = values[order]
    rows = []
    for emitted in range(BITS + 1):
        sources = kernel.sources[incoming][emitted]
        if not sources:
            continue
        capacities = [float(kernel.split_caps[emitted][state]) for state in order]
        source_index = 0
        sink_index = 0
        source_remaining = sources[0][1]
        sink_remaining = capacities[0]
        value = 0.0
        while source_index < len(sources):
            while sink_index < len(capacities) and sink_remaining <= 0.0:
                sink_index += 1
                if sink_index < len(capacities):
                    sink_remaining = capacities[sink_index]
            if sink_index == len(capacities):
                break
            take = min(source_remaining, sink_remaining)
            value += take * sources[source_index][0] * float(ordered_values[sink_index])
            if source_remaining <= sink_remaining:
                sink_remaining -= source_remaining
                source_index += 1
                if source_index < len(sources):
                    source_remaining = sources[source_index][1]
            else:
                source_remaining -= sink_remaining
                sink_index += 1
                if sink_index < len(capacities):
                    sink_remaining = capacities[sink_index]
        contribution = float(kernel.pole_powers[emitted]) * value
        rows.append(
            {
                "emitted_weight": emitted,
                "contribution": contribution,
                "exact_split_row": emitted <= 10 or emitted >= 54,
            }
        )
    rows.sort(key=lambda row: row["contribution"], reverse=True)
    total = math.fsum(float(row["contribution"]) for row in rows)
    cumulative = 0.0
    for row in rows:
        cumulative += float(row["contribution"])
        row["fraction_of_image"] = float(row["contribution"]) / total
        row["cumulative_fraction"] = cumulative / total
    return rows


def combined_constraint_lp(
    sources,
    split_caps,
    pole: float,
    values: np.ndarray,
    incoming: int,
) -> dict:
    """Maximize the point-cap relaxation with both capacity families."""

    variables = []
    by_source = [[] for _ in sources]
    by_row_column: dict[tuple[int, int], list[int]] = {}
    by_column = [[] for _ in range(BITS + 1)]
    objective = []
    for source_index, (emitted, cap, _units) in enumerate(sources):
        for state in range(BITS + 1):
            if not split_caps[emitted][state]:
                continue
            variable = len(variables)
            variables.append((source_index, emitted, state))
            by_source[source_index].append(variable)
            by_row_column.setdefault((emitted, state), []).append(variable)
            by_column[state].append(variable)
            objective.append(-cap * pole**emitted * float(values[state]))

    equality_rows = []
    equality_columns = []
    equality_coefficients = []
    equality_right = []
    row_indices = []
    column_indices = []
    coefficients = []
    right = []
    scale = float(1 << BITS)
    equality = 0
    for source_index, (_emitted, _cap, units) in enumerate(sources):
        for variable in by_source[source_index]:
            equality_rows.append(equality)
            equality_columns.append(variable)
            equality_coefficients.append(1.0)
        equality_right.append(float(units) / scale)
        equality += 1
    constraint = 0
    for (emitted, state), members in sorted(by_row_column.items()):
        for variable in members:
            row_indices.append(constraint)
            column_indices.append(variable)
            coefficients.append(1.0)
        right.append(float(split_caps[emitted][state]) / scale)
        constraint += 1
    for state, members in enumerate(by_column):
        for variable in members:
            row_indices.append(constraint)
            column_indices.append(variable)
            coefficients.append(1.0)
        right.append(float(math.comb(BITS, state)) / scale)
        constraint += 1
    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(constraint, len(variables)),
    ).tocsr()
    equality_matrix = coo_matrix(
        (equality_coefficients, (equality_rows, equality_columns)),
        shape=(equality, len(variables)),
    ).tocsr()
    result = linprog(
        np.asarray(objective, dtype=np.float64),
        A_ub=matrix,
        b_ub=np.asarray(right, dtype=np.float64),
        A_eq=equality_matrix,
        b_eq=np.asarray(equality_right, dtype=np.float64),
        bounds=(0.0, None),
        method="highs-ds",
        options={"primal_feasibility_tolerance": 1e-8, "dual_feasibility_tolerance": 1e-8},
    )
    if not result.success:
        raise RuntimeError("combined capacity LP failed: " + result.message)
    return {
        "incoming_state": incoming,
        "sources": len(sources),
        "variables": len(variables),
        "equality_constraints": equality,
        "capacity_constraints": constraint,
        "maximum_image": -float(result.fun) * scale,
        "total_source_units": sum(int(source[2]) for source in sources),
    }


def pattern_quantile_sources(
    fugacities: np.ndarray, incoming: int, bins: int
) -> list[tuple[int, float, int]]:
    """Upper-round exact drive-pattern weights into top-resolving quantiles."""

    local = compatible_pair_polynomials(GROUP_BITS, fugacities)
    denominator = float(math.comb(BITS, incoming))
    by_emitted: list[list[tuple[float, int]]] = [[] for _ in range(BITS + 1)]
    for counts, emitted_counts in pattern_emission_counts().items():
        polynomial = polynomial_for_pattern(local, counts)
        coefficient = float(polynomial[incoming]) / denominator
        if coefficient <= 0.0:
            continue
        for emitted, units in enumerate(emitted_counts):
            if units:
                by_emitted[emitted].append((coefficient, int(units)))

    quantiles = []
    for emitted, rows in enumerate(by_emitted):
        rows.sort(reverse=True)
        total = sum(units for _cap, units in rows)
        if total != math.comb(BITS, emitted):
            raise RuntimeError("pattern quantiles: emitted-weight mass mismatch")
        current_cap = None
        current_units = 0
        target = 1
        emitted_bins = 0
        for cap, units in rows:
            remaining = units
            while remaining:
                if current_cap is None:
                    current_cap = cap
                take = min(remaining, target - current_units)
                current_units += take
                remaining -= take
                if current_units == target:
                    quantiles.append((emitted, float(current_cap), current_units))
                    emitted_bins += 1
                    current_cap = None
                    current_units = 0
                    # Resolve the dangerous high-cap tail word by word, then
                    # grow rank buckets geometrically.  ``bins`` caps the
                    # number of doublings, not the final bucket count.
                    if emitted_bins < bins:
                        target = min(total, max(1, target * 2))
                    else:
                        target = total
        if current_units:
            quantiles.append((emitted, float(current_cap), current_units))
    return quantiles


def pattern_ratio_sources(
    fugacities: np.ndarray, incoming: int, relative_width: float
) -> list[tuple[int, float, int]]:
    """Group exact pattern masses by cap ratio, rounding each bucket upward."""

    if relative_width <= 0.0:
        raise ValueError("pattern ratio buckets require positive width")
    local = compatible_pair_polynomials(GROUP_BITS, fugacities)
    denominator = float(math.comb(BITS, incoming))
    log_step = math.log1p(relative_width)
    buckets: dict[tuple[int, int], list[float | int]] = {}
    for counts, emitted_counts in pattern_emission_counts().items():
        polynomial = polynomial_for_pattern(local, counts)
        coefficient = float(polynomial[incoming]) / denominator
        if coefficient <= 0.0:
            continue
        cap_bucket = math.floor(math.log(coefficient) / log_step)
        for emitted, units in enumerate(emitted_counts):
            if not units:
                continue
            row = buckets.setdefault((emitted, cap_bucket), [coefficient, 0])
            row[0] = max(float(row[0]), coefficient)
            row[1] = int(row[1]) + int(units)
    sources = [
        (emitted, float(row[0]), int(row[1]))
        for (emitted, _bucket), row in buckets.items()
    ]
    for emitted in range(BITS + 1):
        if sum(units for row_emitted, _cap, units in sources if row_emitted == emitted) != math.comb(BITS, emitted):
            raise RuntimeError("pattern ratio buckets: emitted-weight mass mismatch")
    return sources


def top_mass(sources: list[tuple[float, float]], capacity: float) -> float:
    remaining = max(0.0, capacity)
    result = 0.0
    for score, units in sources:
        take = min(remaining, units)
        result += take * score
        remaining -= take
        if remaining <= 0.0:
            break
    return result


def weighted_mass_capacity_lp(
    row_kernel,
    split_caps,
    pole: float,
    values: np.ndarray,
    incoming: int,
) -> dict:
    """Small LP intersecting row and global weighted top-mass caps."""

    row_sources = []
    all_sources = []
    for emitted in range(BITS + 1):
        rows = [
            (float(cap) * pole**emitted, float(units))
            for cap, units in row_kernel.sources[incoming][emitted]
        ]
        rows.sort(reverse=True)
        row_sources.append(rows)
        all_sources.extend(rows)
    all_sources.sort(reverse=True)

    variables = []
    bounds = []
    by_emitted = [[] for _ in range(BITS + 1)]
    by_state = [[] for _ in range(BITS + 1)]
    objective = []
    for emitted in range(BITS + 1):
        if not row_sources[emitted]:
            continue
        for state in range(BITS + 1):
            capacity = int(split_caps[emitted][state])
            if not capacity:
                continue
            variable = len(variables)
            variables.append((emitted, state))
            by_emitted[emitted].append(variable)
            by_state[state].append(variable)
            bounds.append((0.0, top_mass(row_sources[emitted], float(capacity))))
            objective.append(-float(values[state]))

    eq_rows = []
    eq_columns = []
    eq_values = []
    eq_right = []
    eq_index = 0
    for emitted, members in enumerate(by_emitted):
        if not members:
            continue
        for variable in members:
            eq_rows.append(eq_index)
            eq_columns.append(variable)
            eq_values.append(1.0)
        eq_right.append(top_mass(row_sources[emitted], math.inf))
        eq_index += 1
    ub_rows = []
    ub_columns = []
    ub_values = []
    ub_right = []
    ub_index = 0
    for state, members in enumerate(by_state):
        if not members:
            continue
        for variable in members:
            ub_rows.append(ub_index)
            ub_columns.append(variable)
            ub_values.append(1.0)
        ub_right.append(top_mass(all_sources, float(math.comb(BITS, state))))
        ub_index += 1
    equality = coo_matrix(
        (eq_values, (eq_rows, eq_columns)), shape=(eq_index, len(variables))
    ).tocsr()
    upper = coo_matrix(
        (ub_values, (ub_rows, ub_columns)), shape=(ub_index, len(variables))
    ).tocsr()
    result = linprog(
        np.asarray(objective),
        A_ub=upper,
        b_ub=np.asarray(ub_right),
        A_eq=equality,
        b_eq=np.asarray(eq_right),
        bounds=bounds,
        method="highs-ds",
    )
    if not result.success:
        raise RuntimeError("weighted mass capacity LP failed: " + result.message)
    return {
        "incoming_state": incoming,
        "variables": len(variables),
        "row_equalities": eq_index,
        "global_column_caps": ub_index,
        "maximum_image": -float(result.fun),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--fugacities", required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--cap-ratio-epsilon", type=float, default=0.1)
    parser.add_argument("--skip-pattern-lp", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = parse_vector(args.profile, CLASSES, int)
    fugacities = np.asarray(parse_vector(args.fugacities, CLASSES), dtype=np.float64)

    histograms = block_histograms(GROUP_BITS, fugacities)
    caps = point_caps(GROUP_BITS, fugacities)
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    row_kernel = SharedDriveStratifiedKernel(histograms, caps, split_caps, args.output_pole)
    global_kernel = GlobalColumnKernel(histograms, caps, args.output_pole)

    reports = {}
    row_values = None
    row_worst_state = None
    for name, kernel in (
        ("row_specific", row_kernel),
        ("global_column", global_kernel),
        ("minimum", MinimumUpperKernel(row_kernel, global_kernel)),
    ):
        eigenvalue, domination, values, worst_state = witness(kernel, args.iterations)
        mgf = (
            math.log2(domination)
            + INNER_BLOCKS * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        charge = sum(
            count * math.log2(float(fugacity))
            for count, fugacity in zip(profile, fugacities)
        )
        normalization = float(
            normalization_log2(GROUP_BITS, np.asarray(profile, dtype=np.float64))[0]
        )
        probability = mgf - charge - normalization - D * math.log2(args.output_pole)
        reports[name] = {
            "inner_probability_log2": probability,
            "inner_mgf_log2": mgf,
            "lambda_log2": math.log2(eigenvalue),
            "domination_log2": math.log2(domination),
            "worst_state": worst_state,
        }
        if name == "row_specific":
            row_values = values.copy()
            row_worst_state = worst_state
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_GLOBAL_SYSTEMATIC_COLUMN_BOUND",
        "profile": profile,
        "fugacities": fugacities.tolist(),
        "pole": args.output_pole,
        "branches": reports,
        "row_specific_usage_at_worst_state": row_specific_column_usage(
            row_kernel, row_values, int(row_worst_state)
        ),
        "row_specific_emitted_contributions_at_worst_state": (
            row_specific_emitted_contributions(
                row_kernel, row_values, int(row_worst_state)
            )
        ),
        "weighted_mass_capacity_lp_at_worst_state": weighted_mass_capacity_lp(
            row_kernel,
            split_caps,
            args.output_pole,
            row_values,
            int(row_worst_state),
        ),
    }
    if not args.skip_pattern_lp:
        report["combined_capacity_lp_at_worst_state"] = combined_constraint_lp(
            pattern_ratio_sources(
                fugacities, int(row_worst_state), args.cap_ratio_epsilon
            ),
            split_caps,
            args.output_pole,
            row_values,
            int(row_worst_state),
        )
        report["combined_capacity_lp_at_worst_state"]["row_specific_image"] = float(
            row_kernel.apply(row_values)[int(row_worst_state)]
        )
        report["combined_capacity_lp_at_worst_state"]["one_step_improvement_log2"] = math.log2(
            report["combined_capacity_lp_at_worst_state"]["row_specific_image"]
            / report["combined_capacity_lp_at_worst_state"]["maximum_image"]
        )
        report["combined_capacity_lp_at_worst_state"]["cap_ratio_epsilon"] = (
            args.cap_ratio_epsilon
        )
    report["weighted_mass_capacity_lp_at_worst_state"]["row_specific_image"] = float(
        row_kernel.apply(row_values)[int(row_worst_state)]
    )
    report["weighted_mass_capacity_lp_at_worst_state"]["one_step_improvement_log2"] = math.log2(
        report["weighted_mass_capacity_lp_at_worst_state"]["row_specific_image"]
        / report["weighted_mass_capacity_lp_at_worst_state"]["maximum_image"]
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
