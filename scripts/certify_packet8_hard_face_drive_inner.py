#!/usr/bin/env python3
"""Outward certificate for the hard-face shared-drive inner operator.

The binary64 fugacities and frozen Collatz vector are treated as exact dyadic
rationals.  Local histogram and point-cap tables are generated using one-ulp
upward rounding after every positive IEEE-754 addition, multiplication, and
division.  The final rank-one transportation problems are then solved with
``Fraction`` arithmetic.  Thus the reported Collatz inequality is an outward
upper bound, not a Monte Carlo estimate or an optimized floating-point value.
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
from outward_log2 import log2_fraction, self_check
from probe_packet8_drive_stratified_caps import (
    BITS,
    block_histograms,
    block_histograms_outward,
    exact_point_caps,
    exact_point_caps_outward,
)
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_hard_face_drive_witness import DEFAULT_ANCHORS, DEFAULT_NAME, load_anchor
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


def exact_float(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def outward_transport_image(
    values_float: np.ndarray,
    histograms_upper: np.ndarray,
    caps_upper: np.ndarray,
    split_caps,
    pole_float: float,
) -> tuple[list[Fraction], dict[str, int]]:
    """Evaluate the enlarged transport operator exactly over dyadics.

    Row and column constraints are upper bounds.  If the outward row supplies
    exceed the split-column capacity, the maximizing relaxation simply keeps
    the largest source caps; actual local mass remains feasible because the
    local tables and caps enclose their exact counterparts.
    """

    values = [exact_float(value) for value in values_float]
    state_order = sorted(range(BITS + 1), key=lambda state: values[state], reverse=True)
    pole = exact_float(pole_float)
    pole_powers = [Fraction(1)]
    for _ in range(BITS):
        pole_powers.append(pole_powers[-1] * pole)

    image: list[Fraction] = []
    total_transports = 0
    source_limited = 0
    sink_limited = 0
    for incoming_state in range(BITS + 1):
        denominator = math.comb(BITS, incoming_state)
        total = Fraction(0)
        for emitted in range(BITS + 1):
            sources: list[tuple[Fraction, Fraction]] = []
            for drive_weight in range(BITS + 1):
                histogram = float(histograms_upper[incoming_state, drive_weight, emitted])
                if histogram == 0.0:
                    continue
                cap_float = float(caps_upper[incoming_state, drive_weight])
                if cap_float <= 0.0:
                    raise SystemExit(
                        "hard-face drive inner: outward positive mass has zero cap"
                    )
                cap = exact_float(cap_float)
                mass_upper = exact_float(histogram) / denominator
                sources.append((cap, mass_upper / cap))
            if not sources:
                continue
            sources.sort(key=lambda item: item[0], reverse=True)
            sinks = [
                (values[state], Fraction(int(split_caps[emitted][state])))
                for state in state_order
                if int(split_caps[emitted][state]) > 0
            ]
            if not sinks:
                raise SystemExit(
                    "hard-face drive inner: positive outward mass has no split column"
                )

            source_index = 0
            sink_index = 0
            source_remaining = sources[0][1]
            sink_remaining = sinks[0][1]
            value = Fraction(0)
            while source_index < len(sources) and sink_index < len(sinks):
                take = min(source_remaining, sink_remaining)
                value += take * sources[source_index][0] * sinks[sink_index][0]
                total_transports += 1
                source_remaining -= take
                sink_remaining -= take
                if source_remaining == 0:
                    source_index += 1
                    if source_index < len(sources):
                        source_remaining = sources[source_index][1]
                if sink_remaining == 0:
                    sink_index += 1
                    if sink_index < len(sinks):
                        sink_remaining = sinks[sink_index][1]
            if source_index == len(sources):
                source_limited += 1
            else:
                sink_limited += 1
            total += pole_powers[emitted] * value
        image.append(total)
    return image, {
        "transport_segments": total_transports,
        "source_limited_cells": source_limited,
        "sink_limited_cells": sink_limited,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, default=DEFAULT_ANCHORS)
    parser.add_argument("--name", default=DEFAULT_NAME)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    self_check()

    anchor = load_anchor(args.anchors, args.name)
    fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))

    # Power iteration only proposes a positive test vector.  Soundness comes
    # from the exact outward check below, not from convergence of this loop.
    diagnostic_histograms = block_histograms(fugacities)
    diagnostic_caps = exact_point_caps(fugacities)
    diagnostic_kernel = SharedDriveStratifiedKernel(
        diagnostic_histograms,
        diagnostic_caps,
        split_caps,
        anchor.pole,
    )
    _eigenvalue, _domination, values, _worst = witness(
        diagnostic_kernel, args.iterations
    )

    histograms_upper = block_histograms_outward(fugacities)
    caps_upper = exact_point_caps_outward(fugacities)
    if not np.all(histograms_upper >= diagnostic_histograms):
        raise SystemExit("hard-face drive inner: outward histogram self-check failed")
    if not np.all(caps_upper >= diagnostic_caps):
        raise SystemExit("hard-face drive inner: outward cap self-check failed")

    image, transport_stats = outward_transport_image(
        values,
        histograms_upper,
        caps_upper,
        split_caps,
        anchor.pole,
    )
    value_fractions = [exact_float(value) for value in values]
    ratios = [entry / value for entry, value in zip(image, value_fractions)]
    eigenvalue_upper = max(ratios)
    worst_state = ratios.index(eigenvalue_upper)
    domination = max(Fraction(1, 1) / value for value in value_fractions)
    domination_state = max(
        range(BITS + 1), key=lambda state: Fraction(1, 1) / value_fractions[state]
    )

    log_lambda = log2_fraction(eigenvalue_upper)
    log_domination = log2_fraction(domination)
    log_initial = log2_fraction(value_fractions[0])
    mgf = log_domination + log_lambda.times_int(INNER_BLOCKS) + log_initial
    maximum_local_inflation = max(
        float(histograms_upper[index] / diagnostic_histograms[index])
        for index in zip(*np.nonzero(diagnostic_histograms))
    )
    maximum_cap_inflation = max(
        float(caps_upper[index] / diagnostic_caps[index])
        for index in zip(*np.nonzero(diagnostic_caps))
    )

    report = {
        "status": "OUTWARD_CERTIFIED_SHARED_DRIVE_INNER_COLLATZ",
        "anchor": anchor.name,
        "arithmetic": (
            "exact dyadic fugacities/test vector; one-ulp upward local DP; "
            "exact Fraction transportation; outward Decimal log2"
        ),
        "iterations_used_only_to_propose_vector": args.iterations,
        "outward_float_invariant": (
            "positive dyadic multiplication must remain positive in "
            "binary64; any underflow aborts certification"
        ),
        "lambda_log2_interval": [str(log_lambda.lo), str(log_lambda.hi)],
        "domination_log2_interval": [str(log_domination.lo), str(log_domination.hi)],
        "initial_witness_log2_interval": [str(log_initial.lo), str(log_initial.hi)],
        "inner_mgf_log2_interval": [str(mgf.lo), str(mgf.hi)],
        "inner_mgf_log2_upper": float(mgf.hi),
        "worst_state": worst_state,
        "domination_state": domination_state,
        "maximum_histogram_outward_ratio": maximum_local_inflation,
        "maximum_cap_outward_ratio": maximum_cap_inflation,
        **transport_stats,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
