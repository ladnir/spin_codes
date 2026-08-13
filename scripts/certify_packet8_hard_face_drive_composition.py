#!/usr/bin/env python3
"""End-to-end outward certificate for the first packet-8 hard face.

The outer certificate uses the symmetric p=1/2 three-band linear
Brascamp--Lieb inequality at explicit integer packet variables.  It compares
the physical graph-hole word to the unpunctured 16384-data-word layout by
paying the worst adjacent packet-variable ratio for each of the 128 bit
replacements.  The inner term is the outward shared-drive Collatz certificate.
All logarithmic comparisons use the standard-library outward interval layer.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from certify_packet8_hard_face_drive_inner import exact_float, outward_transport_image
from outward_log2 import LN2, Interval, ln_factorial, log2_fraction, log2_int, self_check
from probe_packet8_drive_stratified_caps import (
    block_histograms,
    block_histograms_outward,
    exact_point_caps,
    exact_point_caps_outward,
)
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_hard_face_drive_witness import DEFAULT_ANCHORS, DEFAULT_NAME, load_anchor
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_three_band_profile_enumerator import D, GROUPS, K, M
from probe_packet8_weight_profile_scalar import CLASSES
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


# Rounded from the symmetric p=1/2 binary64 optimizer, then frozen as exact
# positive integers.  Optimization is not part of the certificate.
OUTER_PACKET_VARIABLES = (
    1,
    86987344005864,
    30731199332287,
    16164074455794,
    1,
    24537302422689,
    25116885168823,
    35605959782668,
    307830577111409,
)
GRAPH_REPLACEMENTS = 128


def exact_convolve(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return result


def outer_s2() -> Fraction:
    packet_polynomial = [
        math.comb(8, weight) * OUTER_PACKET_VARIABLES[weight]
        for weight in range(9)
    ]
    power = [1]
    for _ in range(8):
        power = exact_convolve(power, packet_polynomial)
    if len(power) != 65:
        raise SystemExit("hard-face composition: bad packet polynomial degree")
    return sum(
        Fraction(coefficient * coefficient, math.comb(64, weight) * (1 << 64))
        for weight, coefficient in enumerate(power)
    )


def exact_normalization(profile: tuple[int, ...]) -> Interval:
    natural = ln_factorial(M)
    for count in profile:
        natural = natural - ln_factorial(count)
    result = natural / LN2
    for count, classes in zip(profile, CLASSES):
        if count:
            result = result + log2_int(classes).times_int(count)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, default=DEFAULT_ANCHORS)
    parser.add_argument("--name", default=DEFAULT_NAME)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    self_check()

    anchor = load_anchor(args.anchors, args.name)
    profile = anchor.profile
    if sum(profile) != M:
        raise SystemExit("hard-face composition: profile has wrong mass")

    # Exact symmetric three-band outer witness.
    s2 = outer_s2()
    outer = Interval.exact(K) + log2_fraction(s2).times_int(GROUPS // 2)
    for count, variable in zip(profile, OUTER_PACKET_VARIABLES):
        if count:
            outer = outer - log2_int(variable).times_int(count)
    adjacent_ratio = max(
        Fraction(OUTER_PACKET_VARIABLES[new], OUTER_PACKET_VARIABLES[old])
        for old in range(9)
        for new in range(9)
        if abs(new - old) <= 1
    )
    graph_correction = log2_fraction(adjacent_ratio).times_int(GRAPH_REPLACEMENTS)
    physical_outer = outer + graph_correction

    # Frozen outward inner Collatz witness.
    fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    diagnostic_histograms = block_histograms(fugacities)
    diagnostic_caps = exact_point_caps(fugacities)
    diagnostic_kernel = SharedDriveStratifiedKernel(
        diagnostic_histograms, diagnostic_caps, split_caps, anchor.pole
    )
    _eigenvalue, _domination, values, _worst = witness(
        diagnostic_kernel, args.iterations
    )
    histograms_upper = block_histograms_outward(fugacities)
    caps_upper = exact_point_caps_outward(fugacities)
    image, transport_stats = outward_transport_image(
        values, histograms_upper, caps_upper, split_caps, anchor.pole
    )
    value_fractions = [exact_float(value) for value in values]
    ratios = [entry / value for entry, value in zip(image, value_fractions)]
    eigenvalue_upper = max(ratios)
    domination = max(Fraction(1) / value for value in value_fractions)
    inner_mgf = (
        log2_fraction(domination)
        + log2_fraction(eigenvalue_upper).times_int(INNER_BLOCKS)
        + log2_fraction(value_fractions[0])
    )

    fugacity_charge = Interval.exact(0)
    for count, fugacity in zip(profile, fugacities):
        if count:
            fugacity_charge = fugacity_charge + log2_fraction(
                exact_float(float(fugacity))
            ).times_int(count)
    normalization = exact_normalization(profile)
    distance_penalty = log2_fraction(exact_float(anchor.pole)).times_int(D)
    combined = (
        physical_outer
        + inner_mgf
        - fugacity_charge
        - normalization
        - distance_penalty
    )
    profile_count = log2_int(math.comb(M + 8, 8))
    union_ledger = combined + profile_count
    passed = union_ledger.hi <= -40
    if not passed:
        raise SystemExit(
            "hard-face composition: outward ledger does not close: "
            f"upper={union_ledger.hi}"
        )

    report = {
        "status": "OUTWARD_CERTIFIED_FIRST_HARD_FACE_2^-40",
        "anchor": anchor.name,
        "profile": list(profile),
        "outward_float_invariant": (
            "positive dyadic multiplication must remain positive in "
            "binary64; any underflow aborts certification"
        ),
        "outer_packet_variables": list(OUTER_PACKET_VARIABLES),
        "outer_unpunctured_log2_interval": [str(outer.lo), str(outer.hi)],
        "graph_replacement_adjacent_ratio": str(adjacent_ratio),
        "graph_replacement_log2_interval": [
            str(graph_correction.lo),
            str(graph_correction.hi),
        ],
        "physical_outer_log2_interval": [str(physical_outer.lo), str(physical_outer.hi)],
        "inner_mgf_log2_interval": [str(inner_mgf.lo), str(inner_mgf.hi)],
        "fugacity_charge_log2_interval": [
            str(fugacity_charge.lo),
            str(fugacity_charge.hi),
        ],
        "profile_normalization_log2_interval": [
            str(normalization.lo),
            str(normalization.hi),
        ],
        "distance_pole_term_log2_interval": [
            str(distance_penalty.lo),
            str(distance_penalty.hi),
        ],
        "combined_hard_face_log2_interval": [str(combined.lo), str(combined.hi)],
        "profile_count_log2_interval": [str(profile_count.lo), str(profile_count.hi)],
        "union_ledger_log2_interval": [str(union_ledger.lo), str(union_ledger.hi)],
        "margin_below_2^-40_log2": float(-40 - union_ledger.hi),
        "passed": passed,
        **transport_stats,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
