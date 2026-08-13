#!/usr/bin/env python3
"""Probe fixed positive witnesses at every nonzero packet-weight vertex.

The profile simplex has nine packet-weight classes.  Class zero alone is the
zero message and is excluded, while the eight remaining vertices are useful
anchors for an adaptive cover.  A zero fugacity for an inactive class proves
only the vertex itself; this probe instead assigns every inactive class a
strictly positive floor so each returned witness remains finite on a genuine
neighborhood of that vertex.

The three-band outer inequality and robust 65-state inner relaxation are
rigorous.  Optimization, Perron iteration, and the printed exponents use
binary64 arithmetic and remain diagnostic until outward hardened.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_packet8_repeated_value_state_transfer import (
    EXACT_SLICES,
    SPECTRUM,
    build_split_caps,
    load_ebch128_spectrum,
    load_exact,
    witness,
)
from probe_packet8_small_alphabet import log2_multinomial
from probe_packet8_three_band_linear_bl2 import optimize_outer
from probe_packet8_three_band_profile_enumerator import N
from probe_packet8_weight_profile_scalar import CLASSES
from probe_packet8_weight_profile_state_transfer import (
    D,
    INNER_BLOCKS,
    RobustProfileKernel,
    block_histograms,
)


PACKETS = N // 8


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-poles", default="0.05,0.075,0.1,0.15,0.2,0.3,0.4"
    )
    parser.add_argument("--inactive-floor", type=float, default=1e-12)
    parser.add_argument("--outer-log-fugacity-bound", type=float, default=40.0)
    parser.add_argument("--iterations", type=int, default=64)
    args = parser.parse_args()
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("pure-class frontier: invalid output pole")
    if not 0.0 < args.inactive_floor < 1.0:
        raise SystemExit("pure-class frontier: inactive floor must lie in (0,1)")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))

    print("packet-8 pure packet-weight-class frontier")
    print(
        f"packets={PACKETS} inactive_floor={args.inactive_floor:.12g} "
        f"iterations={args.iterations}"
    )
    for active_weight in range(1, 9):
        profile = [0] * 9
        profile[active_weight] = PACKETS
        normalization_log2 = log2_multinomial(profile) + (
            PACKETS * math.log2(CLASSES[active_weight])
        )

        outer_log, outer_result = optimize_outer(
            profile,
            optimize_band_coefficients=True,
            log_fugacity_bound=args.outer_log_fugacity_bound,
        )
        outer_log2 = outer_log / math.log(2.0)
        band1 = float(outer_result.x[8])

        fugacities = np.full(9, args.inactive_floor, dtype=np.float64)
        fugacities[active_weight] = 1.0
        histograms = block_histograms(fugacities)
        rows = []
        for pole in poles:
            kernel = RobustProfileKernel(
                histograms,
                split_caps,
                pole,
                1.0,
            )
            eigenvalue, domination, values, worst_state = witness(
                kernel, args.iterations
            )
            inner_log2 = (
                math.log2(domination)
                + INNER_BLOCKS * math.log2(eigenvalue)
                + math.log2(float(values[0]))
                - normalization_log2
                - D * math.log2(pole)
            )
            combined = outer_log2 + inner_log2
            rows.append(
                (
                    combined,
                    pole,
                    inner_log2,
                    worst_state,
                    eigenvalue,
                    domination,
                )
            )
        best = min(rows, key=lambda row: row[0])
        print(
            f"weight={active_weight} outer_log2={outer_log2:.12f} "
            f"band_coefficients={1.0-band1:.12f},{band1:.12f},{band1:.12f} "
            f"best_pole={best[1]:.12f} inner_log2={best[2]:.12f} "
            f"combined_log2={best[0]:.12f} worst_state={best[3]}"
        )
    print("status=DIAGNOSTIC_BINARY64_POSITIVE_VERTEX_WITNESSES")


if __name__ == "__main__":
    main()
