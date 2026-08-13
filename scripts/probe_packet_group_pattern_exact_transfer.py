#!/usr/bin/env python3
"""Exact drive-pattern rearrangement probe for the g=4 inner transfer.

For g=4 there are only C(16+5-1,4)=4845 multisets of the sixteen drive-atom
weights.  The compatible-pair mass is constant within each multiset.  This
probe keeps that full mass distribution, computes the exact number of drive
words in each multiset and accumulator-output-weight class, and only then
applies the existing systematic-EBCH split-cap rearrangement.  It therefore
removes the point-cap relaxation while retaining the rigorous split caps.

Binary64 evaluation is diagnostic; an outward implementation is required
before this operator can be used in a certificate.
"""

from __future__ import annotations

import argparse
import bisect
import functools
import itertools
import json
import math
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from packet_group_drive_stratified import (
    BLOCK_BITS,
    accumulator_pattern_counts,
    compatible_pair_polynomials,
)
from packet_group_outer_profile import D, normalization_log2
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


GROUP_BITS = 4
ATOMS = BLOCK_BITS // GROUP_BITS
CLASSES = GROUP_BITS + 1


@functools.lru_cache(maxsize=1)
def pattern_emission_counts() -> dict[tuple[int, ...], tuple[int, ...]]:
    """Map an atom-weight multiset to exact Acc output-weight counts."""

    local = accumulator_pattern_counts(GROUP_BITS)
    transitions = {
        parity: [
            (outgoing, weight, emitted, int(local[parity, outgoing, weight, emitted]))
            for outgoing in range(2)
            for weight in range(CLASSES)
            for emitted in range(CLASSES)
            if local[parity, outgoing, weight, emitted]
        ]
        for parity in range(2)
    }
    zero = (0,) * CLASSES
    current: dict[tuple[int, ...], list[list[int]]] = {
        zero: [[1] + [0] * BLOCK_BITS, [0] * (BLOCK_BITS + 1)]
    }
    for position in range(ATOMS):
        next_rows: dict[tuple[int, ...], list[list[int]]] = {}
        maximum_emitted = position * GROUP_BITS
        for counts, parity_rows in current.items():
            for parity in range(2):
                row = parity_rows[parity]
                for emitted_so_far, multiplicity in enumerate(row[: maximum_emitted + 1]):
                    if not multiplicity:
                        continue
                    for outgoing, weight, emitted, local_count in transitions[parity]:
                        updated = list(counts)
                        updated[weight] += 1
                        key = tuple(updated)
                        target = next_rows.get(key)
                        if target is None:
                            target = [[0] * (BLOCK_BITS + 1), [0] * (BLOCK_BITS + 1)]
                            next_rows[key] = target
                        target[outgoing][emitted_so_far + emitted] += (
                            multiplicity * local_count
                        )
        current = next_rows
    result = {
        counts: tuple(left + right for left, right in zip(rows[0], rows[1]))
        for counts, rows in current.items()
    }
    if len(result) != math.comb(ATOMS + CLASSES - 1, CLASSES - 1):
        raise RuntimeError("pattern-exact transfer: multiset enumeration mismatch")
    if sum(sum(row) for row in result.values()) != 1 << BLOCK_BITS:
        raise RuntimeError("pattern-exact transfer: drive-word mass mismatch")
    return result


def polynomial_for_pattern(local: np.ndarray, counts: tuple[int, ...]) -> np.ndarray:
    polynomial = np.array([1.0], dtype=np.float64)
    powers = []
    for weight, count in enumerate(counts):
        for _ in range(count):
            polynomial = np.convolve(polynomial, local[weight])
    return polynomial


class PatternExactSharedDriveKernel:
    """Shared split-column relaxation with exact drive-pattern supplies."""

    def __init__(self, fugacities: np.ndarray, split_caps, pole: float):
        local = compatible_pair_polynomials(GROUP_BITS, fugacities)
        emissions = pattern_emission_counts()
        raw_sources: list[list[list[tuple[float, int]]]] = [
            [[] for _ in range(BLOCK_BITS + 1)] for _ in range(BLOCK_BITS + 1)
        ]
        for counts, emitted_counts in emissions.items():
            polynomial = polynomial_for_pattern(local, counts)
            for state_weight, coefficient in enumerate(polynomial):
                if coefficient <= 0.0:
                    continue
                cap = float(coefficient) / float(math.comb(BLOCK_BITS, state_weight))
                for emitted, units in enumerate(emitted_counts):
                    if units:
                        raw_sources[state_weight][emitted].append((cap, units))

        self.prefix_units: list[list[tuple[int, ...] | None]] = [
            [None for _ in range(BLOCK_BITS + 1)] for _ in range(BLOCK_BITS + 1)
        ]
        self.prefix_mass: list[list[np.ndarray | None]] = [
            [None for _ in range(BLOCK_BITS + 1)] for _ in range(BLOCK_BITS + 1)
        ]
        self.source_caps: list[list[np.ndarray | None]] = [
            [None for _ in range(BLOCK_BITS + 1)] for _ in range(BLOCK_BITS + 1)
        ]
        for state_weight in range(BLOCK_BITS + 1):
            for emitted in range(BLOCK_BITS + 1):
                sources = sorted(raw_sources[state_weight][emitted], reverse=True)
                if not sources:
                    continue
                source_units = [int(row[1]) for row in sources]
                caps = np.asarray([row[0] for row in sources], dtype=np.float64)
                cumulative_units = tuple(itertools.accumulate(source_units))
                self.prefix_units[state_weight][emitted] = cumulative_units
                self.prefix_mass[state_weight][emitted] = np.cumsum(
                    np.asarray(source_units, dtype=np.float64) * caps
                )
                self.source_caps[state_weight][emitted] = caps
        self.split_caps = split_caps
        self.pole_powers = np.asarray(
            [pole**weight for weight in range(BLOCK_BITS + 1)], dtype=np.float64
        )

    def apply(self, values: np.ndarray) -> np.ndarray:
        result = np.zeros(BLOCK_BITS + 1, dtype=np.float64)
        order = np.argsort(values)[::-1]
        ordered_values = values[order]
        for emitted in range(BLOCK_BITS + 1):
            cumulative_sinks = tuple(
                itertools.accumulate(int(self.split_caps[emitted][state]) for state in order)
            )
            for incoming_state in range(BLOCK_BITS + 1):
                units = self.prefix_units[incoming_state][emitted]
                mass = self.prefix_mass[incoming_state][emitted]
                caps = self.source_caps[incoming_state][emitted]
                if units is None or mass is None or caps is None:
                    continue
                total_mass = float(mass[-1])
                total_units = units[-1]
                filled_rows = []
                for boundary in cumulative_sinks:
                    clipped = min(boundary, total_units)
                    index = min(bisect.bisect_left(units, clipped), len(units) - 1)
                    prior_units = 0 if index == 0 else units[index - 1]
                    prior_mass = 0.0 if index == 0 else float(mass[index - 1])
                    filled_rows.append(
                        min(
                            total_mass,
                            prior_mass + float(clipped - prior_units) * float(caps[index]),
                        )
                    )
                allocations = np.diff(np.asarray([0.0, *filled_rows], dtype=np.float64))
                result[incoming_state] += self.pole_powers[emitted] * float(
                    allocations @ ordered_values
                )
        return result


def parse_vector(text: str, expected: int, cast=float):
    values = [cast(value) for value in text.split(",")]
    if len(values) != expected:
        raise ValueError(f"expected {expected} comma-separated values")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--fugacities", required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = parse_vector(args.profile, CLASSES, int)
    fugacities = np.asarray(parse_vector(args.fugacities, CLASSES), dtype=np.float64)
    if sum(profile) != (1 << 21) // GROUP_BITS or any(value <= 0 for value in profile):
        parser.error("profile must be a full-support g=4 profile of mass 2^19")
    if not 0.0 < args.output_pole < 1.0 or np.any(fugacities <= 0.0):
        parser.error("pole and fugacities must be positive")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    kernel = PatternExactSharedDriveKernel(fugacities, split_caps, args.output_pole)
    eigenvalue, domination, values, worst_state = witness(kernel, args.iterations)
    mgf = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    charge = sum(count * math.log2(float(fugacity)) for count, fugacity in zip(profile, fugacities))
    normalization = float(
        normalization_log2(GROUP_BITS, np.asarray(profile, dtype=np.float64))[0]
    )
    probability = mgf - charge - normalization - D * math.log2(args.output_pole)
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_PATTERN_EXACT_SHARED_DRIVE",
        "profile": profile,
        "fugacities": fugacities.tolist(),
        "pole": args.output_pole,
        "inner_probability_log2": probability,
        "inner_mgf_log2": mgf,
        "lambda_log2": math.log2(eigenvalue),
        "domination_log2": math.log2(domination),
        "worst_state": worst_state,
        "drive_pattern_multisets": len(pattern_emission_counts()),
        "remaining_relaxation": "systematic EBCH split caps and 65-state Perron witness",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
