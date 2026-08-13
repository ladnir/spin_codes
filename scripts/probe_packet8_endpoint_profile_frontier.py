#!/usr/bin/env python3
"""Diagnostic {weight-0,weight-8} profile frontier.

For a grid of exact all-ones-packet densities, combine:

* the three-band Finner outer profile coefficient bound; and
* the repeated-0xff OFF/LIVE conditional transfer optimized over a pole grid.

This isolates the endpoint phase left open by the intermediate profile bound.
Both ingredients use binary64 optimization and the outer model is still the
unpunctured data code, so the frontier is diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_packet8_constant_value_scalar import D, PACKET_SLOTS, log2_binomial, optimize_two_state_x
from probe_packet8_three_band_profile_enumerator import optimize_outer
from scan_packet8_repeated_values import trellis_coefficients


def inner_probability(packets: int, poles: list[float]):
    denominator = log2_binomial(PACKET_SLOTS, packets)
    best = None
    for pole in poles:
        coefficients = trellis_coefficients(0xFF, pole)
        cauchy, x, _sequence, turnoff_state, _turnoff, live_state = (
            optimize_two_state_x(coefficients, packets, 8)
        )
        probability = cauchy - denominator - D * math.log2(pole)
        row = (probability, pole, x, turnoff_state, live_state)
        if best is None or row < best:
            best = row
    assert best is not None
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--denominator", type=int, default=16)
    parser.add_argument("--output-poles", default="0.4,0.5,0.6,0.7,0.8,0.9,0.95")
    args = parser.parse_args()
    if args.denominator <= 1 or PACKET_SLOTS % args.denominator:
        raise SystemExit("endpoint frontier: denominator must divide packet slots")
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("endpoint frontier: poles must lie in (0,1)")

    print("packet-8 endpoint-profile frontier")
    rows = []
    for numerator in range(1, args.denominator):
        packets = PACKET_SLOTS * numerator // args.denominator
        profile = [PACKET_SLOTS - packets, 0, 0, 0, 0, 0, 0, 0, packets]
        coefficient_log, result, _components = optimize_outer(profile)
        outer_log2 = coefficient_log / math.log(2.0)
        probability, pole, x, turnoff_state, live_state = inner_probability(
            packets, poles
        )
        combined = outer_log2 + probability
        row = (combined, numerator, outer_log2, probability, pole, x)
        rows.append(row)
        print(
            f"density={numerator}/{args.denominator} outer_log2={outer_log2:.6f} "
            f"inner_log2={probability:.6f} combined_log2={combined:.6f} "
            f"pole={pole:.6f} x={x:.9f} turnoff_state={turnoff_state} "
            f"live_state={live_state} optimizer_success={result.success}"
        )
    worst = max(rows)
    print(
        f"worst_combined_log2={worst[0]:.6f} "
        f"worst_density={worst[1]}/{args.denominator}"
    )
    print("status=DIAGNOSTIC_ENDPOINT_PROFILE_FRONTIER")


if __name__ == "__main__":
    main()
