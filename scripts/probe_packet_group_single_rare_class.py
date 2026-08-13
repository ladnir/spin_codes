#!/usr/bin/env python3
"""Exact degree-one rare-class extraction for a g=4 inner profile.

For a profile with exactly one class-4 atom, set its fugacity to zero in the
base transfer and form a second block operator for exactly one marked
class-4 input atom.  Sublinearity and homogeneity give

    [x4] T(x4)^B 1 <= B * domination * mu * lambda^(B-1) * v,

when K0(v)<=lambda*v and K1(v)<=mu*v.  This removes the irrelevant all-zero
term that floors the ordinary positive-fugacity Cauchy MGF.

Binary64 evaluation is diagnostic.  The same finite DP and operator
inequalities need outward evaluation before theorem use.
"""

from __future__ import annotations

import argparse
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
)
from packet_group_outer_profile import D, normalization_log2
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


GROUP_BITS = 4
CLASSES = 5
ATOMS = BLOCK_BITS // GROUP_BITS
MARKED_WEIGHT = 4


def local_pair_polynomials(fugacities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    base = np.zeros((CLASSES, CLASSES), dtype=np.float64)
    marked = np.zeros_like(base)
    for drive_weight in range(CLASSES):
        for overlap in range(drive_weight + 1):
            for outside in range(GROUP_BITS - drive_weight + 1):
                input_weight = overlap + outside
                state_weight = drive_weight - overlap + outside
                multiplicity = math.comb(drive_weight, overlap) * math.comb(
                    GROUP_BITS - drive_weight, outside
                )
                if input_weight == MARKED_WEIGHT:
                    marked[drive_weight, state_weight] += multiplicity
                else:
                    base[drive_weight, state_weight] += (
                        multiplicity * float(fugacities[input_weight])
                    )
    return base, marked


def marked_block_histograms(
    fugacities: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    patterns = accumulator_pattern_counts(GROUP_BITS)
    pair_base, pair_marked = local_pair_polynomials(fugacities)
    entries = []
    for incoming in range(2):
        for outgoing in range(2):
            for selected in range(CLASSES):
                for drive in range(CLASSES):
                    for emitted in range(CLASSES):
                        count = int(patterns[incoming, outgoing, drive, emitted])
                        if not count:
                            continue
                        base = count * pair_base[drive, selected]
                        marked = count * pair_marked[drive, selected]
                        if base or marked:
                            entries.append(
                                (incoming, outgoing, selected, drive, emitted, base, marked)
                            )
    shape = (2, BLOCK_BITS + 1, BLOCK_BITS + 1, BLOCK_BITS + 1)
    dp0 = np.zeros(shape, dtype=np.float64)
    dp1 = np.zeros(shape, dtype=np.float64)
    dp0[0, 0, 0, 0] = 1.0
    maximum = 0
    for _position in range(ATOMS):
        next0 = np.zeros_like(dp0)
        next1 = np.zeros_like(dp1)
        for incoming, outgoing, selected, drive, emitted, base, marked in entries:
            source0 = dp0[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            source1 = dp1[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            target = (
                outgoing,
                slice(selected, selected + maximum + 1),
                slice(drive, drive + maximum + 1),
                slice(emitted, emitted + maximum + 1),
            )
            if base:
                next0[target] += base * source0
                next1[target] += base * source1
            if marked:
                next1[target] += marked * source0
        dp0, dp1 = next0, next1
        maximum += GROUP_BITS
    return dp0.sum(axis=0), dp1.sum(axis=0)


def marked_point_caps(fugacities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    local0, local1 = local_pair_polynomials(fugacities)
    maxima0 = np.zeros((BLOCK_BITS + 1, BLOCK_BITS + 1), dtype=np.float64)
    maxima1 = np.zeros_like(maxima0)
    for weights in itertools.combinations_with_replacement(range(CLASSES), ATOMS):
        polynomial0 = np.array([1.0])
        polynomial1 = np.array([0.0])
        for weight in weights:
            next1 = np.convolve(polynomial1, local0[weight]) + np.convolve(
                polynomial0, local1[weight]
            )
            polynomial0 = np.convolve(polynomial0, local0[weight])
            polynomial1 = next1
        drive = sum(weights)
        maxima0[:, drive] = np.maximum(maxima0[:, drive], polynomial0)
        maxima1[:, drive] = np.maximum(maxima1[:, drive], polynomial1)
    for state in range(BLOCK_BITS + 1):
        denominator = float(math.comb(BLOCK_BITS, state))
        maxima0[state] /= denominator
        maxima1[state] /= denominator
    return maxima0, maxima1


def parse_vector(text: str, expected: int, cast=float):
    values = [cast(value) for value in text.split(",")]
    if len(values) != expected:
        raise ValueError(f"expected {expected} comma-separated values")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--base-fugacities", required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = parse_vector(args.profile, CLASSES, int)
    fugacities = np.asarray(
        parse_vector(args.base_fugacities, CLASSES), dtype=np.float64
    )
    if profile[MARKED_WEIGHT] != 1:
        parser.error("single-rare-class probe requires class-4 count exactly one")
    fugacities[MARKED_WEIGHT] = 0.0

    hist0, hist1 = marked_block_histograms(fugacities)
    caps0, caps1 = marked_point_caps(fugacities)
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    kernel0 = SharedDriveStratifiedKernel(hist0, caps0, split_caps, args.output_pole)
    kernel1 = SharedDriveStratifiedKernel(hist1, caps1, split_caps, args.output_pole)
    eigenvalue, domination, values, worst_state = witness(kernel0, args.iterations)
    marked_image = kernel1.apply(values)
    marked_ratio = float(np.max(marked_image / values))
    marked_worst_state = int(np.argmax(marked_image / values))
    mgf = (
        math.log2(domination)
        + (INNER_BLOCKS - 1) * math.log2(eigenvalue)
        + math.log2(INNER_BLOCKS)
        + math.log2(marked_ratio)
        + math.log2(float(values[0]))
    )
    charge = sum(
        count * math.log2(float(fugacity))
        for index, (count, fugacity) in enumerate(zip(profile, fugacities))
        if index != MARKED_WEIGHT and count
    )
    normalization = float(
        normalization_log2(GROUP_BITS, np.asarray(profile, dtype=np.float64))[0]
    )
    probability = mgf - charge - normalization - D * math.log2(args.output_pole)
    report = {
        "status": "DIAGNOSTIC_BINARY64_G4_SINGLE_RARE_CLASS_EXTRACTION",
        "profile": profile,
        "base_fugacities": fugacities.tolist(),
        "pole": args.output_pole,
        "inner_probability_log2": probability,
        "inner_mgf_log2": mgf,
        "base_lambda_log2": math.log2(eigenvalue),
        "base_domination_log2": math.log2(domination),
        "base_worst_state": worst_state,
        "marked_ratio_log2": math.log2(marked_ratio),
        "marked_worst_state": marked_worst_state,
        "insertion_positions_log2": math.log2(INNER_BLOCKS),
        "proof_formula": "B*domination*mu*lambda^(B-1)*v0",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
