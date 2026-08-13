#!/usr/bin/env python3
"""Outward certificate for one generalized-g packet profile witness."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from certify_packet8_hard_face_drive_inner import exact_float, outward_transport_image
from outward_log2 import LN2, Interval, ln_factorial, log2_fraction, log2_int, self_check
from packet_group_drive_stratified import (
    block_histograms,
    block_histograms_outward,
    point_caps,
    point_caps_outward,
    profile_classes,
    profile_count,
)
from packet_group_outer_profile import (
    D,
    GROUPS,
    K,
    N,
    atom_count,
    graph_replacement_ratio,
    optimize_outer,
    symmetric_s2,
)
from packet_group_profile_bound import GRAPH_REPLACEMENTS, split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


def load_row(paths: list[Path], group_bits: int) -> dict:
    candidates = []
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        candidates.extend(data if isinstance(data, list) else [data])
    matches = [row for row in candidates if int(row["group_bits"]) == group_bits]
    if not matches:
        raise SystemExit(f"generalized profile certificate: no g={group_bits} row")
    return matches[-1]


def outward_normalization(group_bits: int, profile: list[int]) -> Interval:
    natural = ln_factorial(atom_count(group_bits))
    for count in profile:
        natural = natural - ln_factorial(count)
    result = natural / LN2
    for count, classes in zip(profile, profile_classes(group_bits)):
        if count:
            result = result + log2_int(classes).times_int(count)
    return result


def integer_outer_variables(group_bits: int, profile: list[int]) -> tuple[int, ...]:
    _value, result = optimize_outer(
        group_bits, profile, optimize_band_coefficients=False
    )
    logs = np.concatenate(([0.0], result.x[:group_bits]))
    logs -= float(np.min(logs))
    # Exact positive integers are the certificate parameters; exp/round only
    # proposes them and has no role in soundness.
    return tuple(max(1, int(round(math.exp(float(value))))) for value in logs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    self_check()
    row = load_row(args.input, args.group)
    group_bits = args.group
    profile = [int(value) for value in row["profile"]]
    fugacities = np.asarray(row["fugacities"], dtype=np.float64)
    pole = float(row["pole"])
    split_caps = split_cap_table()

    diagnostic_histograms = block_histograms(group_bits, fugacities)
    diagnostic_caps = point_caps(group_bits, fugacities)
    kernel = SharedDriveStratifiedKernel(
        diagnostic_histograms, diagnostic_caps, split_caps, pole
    )
    _eigenvalue, _domination, values, _worst = witness(kernel, args.iterations)
    histograms_upper = block_histograms_outward(group_bits, fugacities)
    caps_upper = point_caps_outward(group_bits, fugacities)
    if not np.all(histograms_upper >= diagnostic_histograms):
        raise SystemExit("generalized profile certificate: histogram enclosure failed")
    if not np.all(caps_upper >= diagnostic_caps):
        raise SystemExit("generalized profile certificate: cap enclosure failed")
    image, transport_stats = outward_transport_image(
        values, histograms_upper, caps_upper, split_caps, pole
    )
    value_fractions = [exact_float(value) for value in values]
    eigenvalue_upper = max(
        entry / value for entry, value in zip(image, value_fractions)
    )
    domination = max(Fraction(1) / value for value in value_fractions)
    inner_mgf = (
        log2_fraction(domination)
        + log2_fraction(eigenvalue_upper).times_int(INNER_BLOCKS)
        + log2_fraction(value_fractions[0])
    )

    variables = integer_outer_variables(group_bits, profile)
    outer = Interval.exact(K) + log2_fraction(
        symmetric_s2(group_bits, variables)
    ).times_int(GROUPS // 2)
    for count, variable in zip(profile, variables):
        if count:
            outer = outer - log2_int(variable).times_int(count)
    graph = log2_fraction(graph_replacement_ratio(variables)).times_int(
        GRAPH_REPLACEMENTS
    )

    charge = Interval.exact(0)
    for count, fugacity in zip(profile, fugacities):
        if count:
            charge = charge + log2_fraction(
                exact_float(float(fugacity))
            ).times_int(count)
    normalization = outward_normalization(group_bits, profile)
    distance = log2_fraction(exact_float(pole)).times_int(D)
    combined = outer + graph + inner_mgf - charge - normalization - distance
    count_log = log2_int(profile_count(group_bits, N))
    union = combined + count_log
    margin = -40 - union.hi
    report = {
        "status": "OUTWARD_CERTIFIED_GENERALIZED_G_MAPPED_PROFILE",
        "group_bits": group_bits,
        "profile": profile,
        "outer_variables": list(variables),
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "outward_float_invariant": (
            "positive dyadic multiplication must remain positive in "
            "binary64; any underflow aborts certification"
        ),
        "outer_log2_interval": [str(outer.lo), str(outer.hi)],
        "graph_log2_interval": [str(graph.lo), str(graph.hi)],
        "inner_mgf_log2_interval": [str(inner_mgf.lo), str(inner_mgf.hi)],
        "charge_log2_interval": [str(charge.lo), str(charge.hi)],
        "normalization_log2_interval": [str(normalization.lo), str(normalization.hi)],
        "distance_log2_interval": [str(distance.lo), str(distance.hi)],
        "combined_log2_interval": [str(combined.lo), str(combined.hi)],
        "profile_count_log2_interval": [str(count_log.lo), str(count_log.hi)],
        "union_log2_interval": [str(union.lo), str(union.hi)],
        "certified_margin_log2": float(margin),
        "passed": bool(margin >= 0),
        **transport_stats,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
