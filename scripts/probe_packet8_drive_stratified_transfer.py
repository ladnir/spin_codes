#!/usr/bin/env python3
"""Drive-stratified 65-state transfer for one packet-weight profile.

This replaces the robust kernel's single ``max(f)^8`` point-mass cap by the
exact-local cap c[q,d] for incoming state weight q and accumulator-drive
weight d.  Each drive stratum is pessimistically allowed to reuse the entire
systematic EBCH split cap.  Consequently the operator is a valid relaxation
of the exact local map even though it does not yet share column capacities
between drive weights.

Local tables and cap coverage are exact finite combinatorics evaluated in
binary64.  Perron iteration, fugacities, and logarithms remain diagnostic
until converted to outward or rational arithmetic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_drive_stratified_caps import BITS, block_histograms, exact_point_caps
from packet_group_native import load_shared_drive_apply
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_small_alphabet import log2_multinomial
from probe_packet8_weight_profile_scalar import CLASSES, parse_profile
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import (
    D,
    INNER_BLOCKS,
    RobustProfileKernel,
)


class DriveStratifiedKernel:
    def __init__(self, histograms, point_caps, split_caps, pole: float):
        self.histograms = histograms
        self.point_caps = point_caps
        self.split_caps = split_caps
        self.pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
        self.denominators = np.array(
            [float(math.comb(BITS, state)) for state in range(BITS + 1)]
        )

    def apply(self, values: np.ndarray) -> np.ndarray:
        result = np.zeros(BITS + 1, dtype=np.float64)
        order = np.argsort(values)[::-1]
        scores = values[order]
        for incoming_state in range(BITS + 1):
            denominator = self.denominators[incoming_state]
            total = 0.0
            for drive_weight in range(BITS + 1):
                point_cap = float(self.point_caps[incoming_state, drive_weight])
                if point_cap == 0.0:
                    if np.any(self.histograms[incoming_state, drive_weight]):
                        raise SystemExit(
                            "drive-stratified transfer: positive mass has zero point cap"
                        )
                    continue
                for emitted in range(BITS + 1):
                    mass = float(
                        self.histograms[incoming_state, drive_weight, emitted]
                    ) / denominator
                    if mass == 0.0:
                        continue
                    capacities = np.array(
                        [
                            point_cap * float(self.split_caps[emitted][next_state])
                            for next_state in order
                        ]
                    )
                    positive = capacities > 0.0
                    capacities = capacities[positive]
                    active_scores = scores[positive]
                    cumulative = np.cumsum(capacities)
                    weighted = np.cumsum(capacities * active_scores)
                    if len(cumulative) == 0:
                        raise SystemExit(
                            "drive-stratified transfer: positive mass has no split capacity"
                        )
                    tolerance = 3e-11 * max(1.0, float(cumulative[-1]))
                    if mass > float(cumulative[-1]) + tolerance:
                        raise SystemExit(
                            "drive-stratified transfer: split caps do not cover stratum "
                            f"q={incoming_state} d={drive_weight} y={emitted} "
                            f"mass={mass} capacity={float(cumulative[-1])}"
                        )
                    index = min(
                        int(np.searchsorted(cumulative, mass, side="left")),
                        len(cumulative) - 1,
                    )
                    prior_mass = 0.0 if index == 0 else float(cumulative[index - 1])
                    prior_weighted = 0.0 if index == 0 else float(weighted[index - 1])
                    greedy = prior_weighted + (mass - prior_mass) * float(
                        active_scores[index]
                    )
                    total += float(self.pole_powers[emitted]) * greedy
            result[incoming_state] = total
        return result


class SharedDriveStratifiedKernel:
    """Share each BCH split-capacity column across all drive weights.

    Put n[d,r]=x[d,r]/c[q,d].  Row d has fixed supply
    mass[q,d,y]/c[q,d], column r has capacity split[y,r], and the profit is
    c[q,d]*values[r].  This rank-one transportation problem is maximized by
    monotonically pairing the sorted caps and sorted witness values.
    """

    def __init__(self, histograms, point_caps, split_caps, pole: float):
        self.split_caps = split_caps
        self.pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
        self.sources: list[list[list[tuple[float, float]]]] = [
            [[] for _ in range(BITS + 1)] for _ in range(BITS + 1)
        ]
        for incoming_state in range(BITS + 1):
            denominator = float(math.comb(BITS, incoming_state))
            for emitted in range(BITS + 1):
                sources = []
                for drive_weight in range(BITS + 1):
                    mass = float(
                        histograms[incoming_state, drive_weight, emitted]
                    ) / denominator
                    if mass == 0.0:
                        continue
                    cap = float(point_caps[incoming_state, drive_weight])
                    if cap <= 0.0:
                        raise SystemExit(
                            "shared drive transfer: positive mass has zero point cap"
                        )
                    sources.append((cap, mass / cap))
                sources.sort(reverse=True)
                self.sources[incoming_state][emitted] = sources

        self.native_apply = load_shared_drive_apply()
        self.native_source_caps = None
        self.native_source_units = None
        self.native_source_counts = None
        self.native_split_caps = None
        if self.native_apply is not None:
            shape = (BITS + 1, BITS + 1, BITS + 1)
            source_caps = np.zeros(shape, dtype=np.float64)
            source_units = np.zeros(shape, dtype=np.float64)
            source_counts = np.zeros((BITS + 1, BITS + 1), dtype=np.int32)
            for incoming_state in range(BITS + 1):
                for emitted in range(BITS + 1):
                    sources = self.sources[incoming_state][emitted]
                    source_counts[incoming_state, emitted] = len(sources)
                    for source_index, (cap, units) in enumerate(sources):
                        source_caps[incoming_state, emitted, source_index] = cap
                        source_units[incoming_state, emitted, source_index] = units
            self.native_source_caps = np.ascontiguousarray(source_caps.ravel())
            self.native_source_units = np.ascontiguousarray(source_units.ravel())
            self.native_source_counts = np.ascontiguousarray(source_counts.ravel())
            self.native_split_caps = np.ascontiguousarray(
                np.asarray(split_caps, dtype=np.float64).reshape(-1)
            )

    def _apply_native(self, values: np.ndarray) -> np.ndarray:
        order = np.ascontiguousarray(np.argsort(values)[::-1], dtype=np.int32)
        contiguous_values = np.ascontiguousarray(values, dtype=np.float64)
        result = np.empty(BITS + 1, dtype=np.float64)
        status = self.native_apply(
            self.native_source_caps,
            self.native_source_units,
            self.native_source_counts,
            self.native_split_caps,
            np.ascontiguousarray(self.pole_powers),
            contiguous_values,
            order,
            result,
        )
        if status:
            raise RuntimeError(
                "native shared drive transfer: split columns do not cover mass"
            )
        return result

    def apply(self, values: np.ndarray) -> np.ndarray:
        if self.native_apply is not None:
            return self._apply_native(values)
        result = np.zeros(BITS + 1, dtype=np.float64)
        order = np.argsort(values)[::-1]
        ordered_values = values[order]
        ordered_capacities = [
            [float(self.split_caps[emitted][state]) for state in order]
            for emitted in range(BITS + 1)
        ]
        for incoming_state in range(BITS + 1):
            total = 0.0
            for emitted in range(BITS + 1):
                sources = self.sources[incoming_state][emitted]
                if not sources:
                    continue
                source_mass = sum(cap * units for cap, units in sources)
                capacities = ordered_capacities[emitted]
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
                        residual_mass = source_remaining * sources[source_index][0]
                        residual_mass += sum(
                            cap * units
                            for cap, units in sources[source_index + 1 :]
                        )
                        if residual_mass > 1e-10 * max(1.0, source_mass):
                            raise SystemExit(
                                "shared drive transfer: split columns do not cover mass "
                                f"q={incoming_state} y={emitted} "
                                f"residual_mass={residual_mass}"
                            )
                        # Binary64 exhaustion residue: assign it to the largest
                        # witness value, which is the pessimistic repair.
                        value += residual_mass * float(ordered_values[0])
                        source_index = len(sources)
                        break
                    take = min(source_remaining, sink_remaining)
                    value += (
                        take
                        * sources[source_index][0]
                        * float(ordered_values[sink_index])
                    )
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
                total += float(self.pole_powers[emitted]) * value
            result[incoming_state] = total
        return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--fugacities", required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument(
        "--capacity-mode", choices=("independent", "shared"), default="shared"
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"drive-stratified transfer: {error}") from error
    fugacities = np.array(
        [float(value) for value in args.fugacities.split(",")], dtype=np.float64
    )
    if (
        len(fugacities) != 9
        or np.any(fugacities < 0.0)
        or fugacities[0] <= 0.0
        or not 0.0 < args.output_pole < 1.0
    ):
        raise SystemExit("drive-stratified transfer: invalid parameters")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    histograms = block_histograms(fugacities)
    point_caps = exact_point_caps(fugacities)
    kernel_type = (
        DriveStratifiedKernel
        if args.capacity_mode == "independent"
        else SharedDriveStratifiedKernel
    )
    kernel = kernel_type(histograms, point_caps, split_caps, args.output_pole)
    normalization = log2_multinomial(profile) + sum(
        count * math.log2(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )
    charge = sum(
        profile[weight] * math.log2(float(fugacities[weight]))
        for weight in range(9)
        if profile[weight]
    )
    def evaluate(active_kernel):
        eigenvalue, domination, values, worst_state = witness(
            active_kernel, args.iterations
        )
        mgf = (
            math.log2(domination)
            + INNER_BLOCKS * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        probability = (
            mgf - charge - normalization - D * math.log2(args.output_pole)
        )
        return probability, eigenvalue, domination, worst_state

    probability, eigenvalue, domination, worst_state = evaluate(kernel)
    baseline_kernel = RobustProfileKernel(
        np.sum(histograms, axis=1),
        split_caps,
        args.output_pole,
        float(np.max(fugacities)) ** 8,
    )
    baseline_probability, _base_eigenvalue, _base_domination, _base_worst = evaluate(
        baseline_kernel
    )

    print("packet-8 drive-stratified 65-state transfer probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"pole={args.output_pole:.12f}")
    print(f"capacity_mode={args.capacity_mode}")
    print(f"baseline_robust_probability_log2={baseline_probability:.12f}")
    print(f"conditional_probability_log2={probability:.12f}")
    print(f"conditional_probability_per_group={probability / INNER_BLOCKS:.12f}")
    print(
        f"improvement_log2={baseline_probability - probability:.12f} "
        f"improvement_per_group={(baseline_probability - probability) / INNER_BLOCKS:.12f}"
    )
    print(f"lambda_log2={math.log2(eigenvalue):.12f}")
    print(f"domination_log2={math.log2(domination):.12f}")
    print(f"worst_state={worst_state}")
    print("status=DIAGNOSTIC_EXACT_LOCAL_DRIVE_STRATIFIED_BINARY64")


if __name__ == "__main__":
    main()
