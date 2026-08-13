#!/usr/bin/env python3
"""Direct fugacity tuning of the robust 65-state profile transfer.

The base 65-state probe evaluates fugacities optimized for the coarser
OFF/LIVE transfer.  This script instead minimizes the robust kernel's own
Cauchy objective at one fixed output pole by deterministic coordinate search.
It is an optimization diagnostic; the exact-local kernel and split-cap
relaxation are unchanged.
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
from probe_packet8_weight_profile_scalar import CLASSES, optimize_profile, parse_profile
from probe_packet8_weight_profile_state_transfer import (
    D,
    INNER_BLOCKS,
    RobustProfileKernel,
    block_histograms,
)
from probe_packet8_small_alphabet import log2_multinomial


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--coordinate-iterations", type=int, default=12)
    parser.add_argument("--witness-iterations", type=int, default=48)
    parser.add_argument(
        "--inactive-fugacity-floor",
        type=float,
        default=0.0,
        help="fixed positive fugacity assigned to inactive nonzero classes",
    )
    parser.add_argument(
        "--start-fugacities",
        help="optional comma-separated nine-value warm start with class zero fixed to one",
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"tuned state transfer: {error}") from error
    if not 0.0 < args.output_pole < 1.0:
        raise SystemExit("tuned state transfer: output pole must lie in (0,1)")
    if not 0.0 <= args.inactive_fugacity_floor < 1.0:
        raise SystemExit("tuned state transfer: invalid inactive fugacity floor")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    normalization = log2_multinomial(profile) + sum(
        count * math.log2(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )
    _scalar, _sequence, _turnoff, _live, initial = optimize_profile(
        profile, args.output_pole
    )
    if args.start_fugacities is not None:
        initial = np.array(
            [float(value) for value in args.start_fugacities.split(",")],
            dtype=np.float64,
        )
        if (
            len(initial) != 9
            or abs(float(initial[0]) - 1.0) > 1e-12
            or np.any(initial <= 0.0)
        ):
            raise SystemExit("tuned state transfer: invalid warm-start fugacities")
    active = [weight for weight in range(1, 9) if profile[weight]]

    evaluations = 0

    def objective(point: np.ndarray):
        nonlocal evaluations
        evaluations += 1
        fugacities = np.full(
            9, args.inactive_fugacity_floor, dtype=np.float64
        )
        fugacities[0] = 1.0
        fugacities[active] = np.exp(point)
        histograms = block_histograms(fugacities)
        kernel = RobustProfileKernel(
            histograms,
            split_caps,
            args.output_pole,
            float(np.max(fugacities)) ** 8,
        )
        eigenvalue, domination, values, worst_state = witness(
            kernel, args.witness_iterations
        )
        mgf = (
            math.log2(domination)
            + INNER_BLOCKS * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        charge = sum(
            profile[weight] * math.log2(float(fugacities[weight]))
            for weight in active
        )
        probability = (
            mgf
            - charge
            - normalization
            - D * math.log2(args.output_pole)
        )
        return probability, fugacities, eigenvalue, domination, worst_state

    point = np.log(initial[active])
    initial_row = objective(point)
    best_row = initial_row
    for _pass in range(args.passes):
        old = point.copy()
        for coordinate in range(len(point)):
            low = -10.0
            high = 2.0
            for _ in range(args.coordinate_iterations):
                left = (2.0 * low + high) / 3.0
                right = (low + 2.0 * high) / 3.0
                left_point = point.copy()
                right_point = point.copy()
                left_point[coordinate] = left
                right_point[coordinate] = right
                left_value = objective(left_point)[0]
                right_value = objective(right_point)[0]
                if left_value <= right_value:
                    high = right
                else:
                    low = left
            point[coordinate] = (low + high) / 2.0
            candidate = objective(point)
            if candidate[0] < best_row[0]:
                best_row = candidate
        if float(np.max(np.abs(point - old))) < 1e-5:
            break

    final_row = objective(point)
    if best_row[0] < final_row[0]:
        final_row = best_row
    probability, fugacities, eigenvalue, domination, worst_state = final_row
    print("packet-weight-profile tuned 65-state transfer probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"pole={args.output_pole:.12f} evaluations={evaluations}")
    print(f"initial_probability_log2={initial_row[0]:.12f}")
    print(f"tuned_probability_log2={probability:.12f}")
    print(f"improvement_bits={initial_row[0]-probability:.12f}")
    print(f"lambda_log2={math.log2(eigenvalue):.12f}")
    print(f"domination_log2={math.log2(domination):.12f}")
    print(f"worst_state={worst_state}")
    print("fugacities=" + ",".join(f"{value:.10g}" for value in fugacities))
    print("status=DIAGNOSTIC_DIRECT_ROBUST_FUGACITY_TUNING")


if __name__ == "__main__":
    main()
