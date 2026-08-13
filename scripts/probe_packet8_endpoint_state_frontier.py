#!/usr/bin/env python3
"""65-state endpoint frontier with the structured outer count.

For every requested exact density of ``0xff`` packets among zero packets:

1. use the fast byte trellis to choose the best scalar OFF/LIVE pole and
   fugacity from a fixed pole grid;
2. evaluate that same valid Cauchy point with the stronger 65-state exact-
   accumulator/robust-split operator; and
3. add the universal degree-8-Finner constant-packet outer count.

The outer count and all local combinatorics are rigorous reductions, but the
stored universal exponent and state-operator arithmetic are binary64
diagnostics.  This scans a frontier; it is not yet an interval certificate.
"""

from __future__ import annotations

import argparse
import math

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_constant_value_scalar import (
    D,
    INNER_BLOCKS,
    PACKET_SLOTS,
    log2_binomial,
    optimize_two_state_x,
)
from probe_packet8_repeated_value_state_transfer import (
    RobustRepeatedKernel,
    build_histograms,
    witness,
)
from probe_packet8_weight_transfer import build_split_caps, load_exact
from scan_packet8_repeated_values import trellis_coefficients


OUTER_LOG2_UPPER = 4096.001


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--denominator", type=int, default=32)
    parser.add_argument(
        "--output-poles", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98"
    )
    parser.add_argument(
        "--numerators",
        help="comma-separated density numerators to scan (default: every interior row)",
    )
    parser.add_argument(
        "--robust-all-poles",
        action="store_true",
        help="evaluate the 65-state bound at every pole instead of scalar screening",
    )
    parser.add_argument("--iterations", type=int, default=64)
    args = parser.parse_args()
    if args.denominator <= 1 or PACKET_SLOTS % args.denominator:
        raise SystemExit("endpoint state frontier: denominator must divide packet slots")
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("endpoint state frontier: poles must lie in (0,1)")
    numerators = (
        [int(value) for value in args.numerators.split(",")]
        if args.numerators
        else list(range(1, args.denominator))
    )
    if any(not 0 < numerator < args.denominator for numerator in numerators):
        raise SystemExit("endpoint state frontier: numerators must be interior")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    histograms = build_histograms(0xFF)
    coefficient_tables = {
        pole: trellis_coefficients(0xFF, pole) for pole in poles
    }

    rows = []
    print("packet-8 endpoint 65-state frontier")
    for numerator in numerators:
        packets = PACKET_SLOTS * numerator // args.denominator
        denominator_log2 = log2_binomial(PACKET_SLOTS, packets)
        scalar_rows = []
        for pole in poles:
            cauchy, x, _sequence, _turnoff_state, _turnoff, _live_state = (
                optimize_two_state_x(coefficient_tables[pole], packets, 8)
            )
            scalar_probability = (
                cauchy - denominator_log2 - D * math.log2(pole)
            )
            scalar_rows.append((scalar_probability, pole, x))
        candidates = scalar_rows if args.robust_all_poles else [min(scalar_rows)]
        robust_best = None
        for _scalar_probability, pole, x in candidates:
            kernel = RobustRepeatedKernel(histograms, split_caps, pole, x)
            eigenvalue, domination, values, worst_state = witness(
                kernel, args.iterations
            )
            mgf_log2 = (
                math.log2(domination)
                + INNER_BLOCKS * math.log2(eigenvalue)
                + math.log2(float(values[0]))
            )
            robust_probability = (
                mgf_log2
                - packets * math.log2(x)
                - denominator_log2
                - D * math.log2(pole)
            )
            robust_row = (robust_probability, pole, x, worst_state)
            if robust_best is None or robust_row < robust_best:
                robust_best = robust_row
        assert robust_best is not None
        robust_probability, pole, x, worst_state = robust_best
        combined = OUTER_LOG2_UPPER + robust_probability
        rows.append((combined, numerator, robust_probability, pole, x, worst_state))
        print(
            f"density={numerator}/{args.denominator} "
            f"robust_inner_log2={robust_probability:.6f} "
            f"combined_log2={combined:.6f} pole={pole:.6f} x={x:.9f} "
            f"worst_state={worst_state}"
        )
    worst = max(rows)
    print(
        f"worst_combined_log2={worst[0]:.6f} "
        f"worst_density={worst[1]}/{args.denominator}"
    )
    print("status=DIAGNOSTIC_ENDPOINT_65_STATE_OUTER_COMBINED_FRONTIER")


if __name__ == "__main__":
    main()
