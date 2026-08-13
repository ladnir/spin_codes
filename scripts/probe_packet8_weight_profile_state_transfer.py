#!/usr/bin/env python3
"""65-state robust transfer for an exact packet-weight profile.

For fixed class fugacities, enumerate one eight-bit packet jointly over its
concrete byte value and the selected incoming-state bits.  An exact packet DP
records selected-state count, incoming/outgoing accumulator parity, and
emitted weight.  Convolving eight packets gives, for every incoming state
weight, the exact fugacity-weighted accumulator histogram for one 64-bit
recursive block.

The rigorous systematic EBCH split caps then define the same monotone greedy
65-state operator used by the repeated-value probe.  Fugacities are seeded by
the OFF/LIVE multivariate optimization.  Local combinatorics are exact before
binary64 evaluation; the pole/fugacity choice, greedy curves, and nonlinear
Perron witness remain diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_constant_value_scalar import BITS, D, INNER_BLOCKS
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_small_alphabet import log2_multinomial
from probe_packet8_weight_profile_scalar import (
    CLASSES,
    optimize_profile,
    parse_profile,
)
from probe_packet8_weight_transfer import build_split_caps, load_exact


PACKET_BITS = 8


def one_packet_transfer(fugacities: np.ndarray) -> np.ndarray:
    transfer = np.zeros((2, 2, PACKET_BITS + 1, PACKET_BITS + 1), dtype=np.float64)
    for value in range(1 << PACKET_BITS):
        value_weight = float(fugacities[value.bit_count()])
        if value_weight == 0.0:
            continue
        for selected_mask in range(1 << PACKET_BITS):
            selected = selected_mask.bit_count()
            for incoming_parity in range(2):
                parity = incoming_parity
                emitted = 0
                for position in range(PACKET_BITS):
                    parity ^= (value >> position & 1) ^ (
                        selected_mask >> position & 1
                    )
                    emitted += parity
                transfer[
                    incoming_parity, parity, selected, emitted
                ] += value_weight
    return transfer


def block_histograms(fugacities: np.ndarray) -> np.ndarray:
    packet = one_packet_transfer(fugacities)
    entries = [
        (incoming, outgoing, selected, emitted, float(packet[incoming, outgoing, selected, emitted]))
        for incoming in range(2)
        for outgoing in range(2)
        for selected in range(PACKET_BITS + 1)
        for emitted in range(PACKET_BITS + 1)
        if packet[incoming, outgoing, selected, emitted]
    ]
    dp = np.zeros((2, BITS + 1, BITS + 1), dtype=np.float64)
    dp[0, 0, 0] = 1.0
    maximum = 0
    for _packet_index in range(8):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1]
            next_dp[
                outgoing,
                selected : selected + maximum + 1,
                emitted : emitted + maximum + 1,
            ] += coefficient * source
        dp = next_dp
        maximum += PACKET_BITS
    histograms = dp.sum(axis=0)
    # For each state weight, summing over output weights gives C(64,s) times
    # the eighth power of the one-packet concrete-value fugacity mass.
    packet_mass = sum(
        classes * float(fugacity)
        for classes, fugacity in zip(CLASSES, fugacities)
    )
    for state in range(BITS + 1):
        expected = math.comb(BITS, state) * packet_mass**8
        actual = float(np.sum(histograms[state]))
        if abs(actual - expected) > 2e-12 * max(1.0, expected):
            raise SystemExit("weight-profile state transfer: histogram mass mismatch")
    return histograms


class RobustProfileKernel:
    def __init__(self, histograms, split_caps, pole: float, maximum_value_weight: float):
        self.histograms = histograms
        self.split_caps = split_caps
        self.pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
        self.denominators = np.array(
            [float(math.comb(BITS, state)) for state in range(BITS + 1)]
        )
        self.maximum_value_weight = maximum_value_weight

    def apply(self, values: np.ndarray) -> np.ndarray:
        moments = np.zeros(BITS + 1, dtype=np.float64)
        order = np.argsort(values)[::-1]
        scores = values[order]
        for emitted in range(BITS + 1):
            # After averaging the uniform incoming state support, the weight
            # attached to any one emitted 64-bit value is at most the maximum
            # fugacity weight of a concrete current input.  Scale the exact
            # split capacities by that point-mass cap before greedy filling.
            counts = self.histograms[:, emitted] / self.denominators
            if not np.any(counts):
                continue
            capacities = np.array(
                [
                    self.maximum_value_weight
                    * float(self.split_caps[emitted][state])
                    for state in order
                ]
            )
            positive = capacities > 0.0
            active_capacities = capacities[positive]
            active_scores = scores[positive]
            cumulative = np.cumsum(active_capacities)
            weighted = np.cumsum(active_capacities * active_scores)
            if len(cumulative) == 0:
                raise SystemExit(
                    "weight-profile state transfer: positive mass has no split capacity"
                )
            capacity_slack = 2e-12 * max(1.0, float(cumulative[-1]))
            if float(np.max(counts)) > float(cumulative[-1]) + capacity_slack:
                raise SystemExit(
                    "weight-profile state transfer: scaled split capacities miss mass"
                )
            indices = np.searchsorted(cumulative, counts, side="left")
            indices = np.minimum(indices, len(cumulative) - 1)
            prior_indices = np.maximum(indices - 1, 0)
            prior_mass = np.where(indices == 0, 0.0, cumulative[prior_indices])
            prior_weighted = np.where(indices == 0, 0.0, weighted[prior_indices])
            curve = prior_weighted + (counts - prior_mass) * active_scores[indices]
            moments += self.pole_powers[emitted] * curve
        return moments


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output-poles", default="0.5,0.6,0.7,0.8,0.9")
    parser.add_argument("--iterations", type=int, default=64)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"weight-profile state transfer: {error}") from error
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("weight-profile state transfer: invalid pole")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    normalization_log2 = log2_multinomial(profile) + sum(
        count * math.log2(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )
    rows = []
    print("packet-weight-profile 65-state robust transfer probe")
    print(
        f"profile={','.join(map(str, profile))} "
        f"normalization_log2={normalization_log2:.9f}"
    )
    for pole in poles:
        _scalar, _sequence, _turnoff, _live, fugacities = optimize_profile(
            profile, pole
        )
        histograms = block_histograms(fugacities)
        maximum_value_weight = float(np.max(fugacities)) ** 8
        kernel = RobustProfileKernel(
            histograms, split_caps, pole, maximum_value_weight
        )
        eigenvalue, domination, values, worst_state = witness(kernel, args.iterations)
        mgf_log2 = (
            math.log2(domination)
            + INNER_BLOCKS * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        cauchy_charge = sum(
            count * math.log2(float(fugacity))
            for count, fugacity in zip(profile, fugacities)
            if count
        )
        probability = (
            mgf_log2 - cauchy_charge - normalization_log2 - D * math.log2(pole)
        )
        row = (probability, pole, fugacities, worst_state, eigenvalue, domination)
        rows.append(row)
        print(
            f"pole={pole:.6f} probability_log2={probability:.6f} "
            f"lambda_log2={math.log2(eigenvalue):.12f} "
            f"domination_log2={math.log2(domination):.12f} "
            f"worst_state={worst_state} "
            f"fugacities={','.join(f'{value:.8g}' for value in fugacities)}"
        )
    best = min(rows)
    print(
        f"best_probability_log2={best[0]:.6f} best_pole={best[1]:.6f} "
        f"best_worst_state={best[3]}"
    )
    print("status=DIAGNOSTIC_EXACT_LOCAL_PROFILE_65_STATE_BINARY64")


if __name__ == "__main__":
    main()
