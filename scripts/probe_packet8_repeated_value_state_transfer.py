#!/usr/bin/env python3
"""65-state robust transfer for one repeated concrete packet value.

For every incoming state weight and every subset of the eight current packet
slots, ``accumulator_histograms`` gives the exact number of permuted state
supports producing each emitted weight.  For a fixed emitted weight ``y``, the
systematic EBCH split table gives exact or rigorous caps on the number of
``y``-weight values producing each next-state weight ``q``.  Greedily assigning
the exact accumulator count to the largest ``z^y v_q`` scores under those caps
defines a monotone homogeneous 65-state upper operator.

This retains state-weight reachability that the OFF/LIVE probe discards.  The
local counts and split caps are rigorous; binary64 curve evaluation, nonlinear
power iteration, and the selected pole/fugacity make the result diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_constant_value_scalar import (
    BITS,
    D,
    INNER_BLOCKS,
    PACKET_SLOTS,
    accumulator_histograms,
    log2_binomial,
)
from probe_packet8_weight_transfer import build_split_caps, load_exact


PACKETS_PER_INNER = 8


def build_histograms(packet_value: int) -> np.ndarray:
    histograms = np.zeros(
        (BITS + 1, 1 << PACKETS_PER_INNER, BITS + 1), dtype=np.uint64
    )
    for mask in range(1 << PACKETS_PER_INNER):
        current = 0
        for slot in range(PACKETS_PER_INNER):
            if mask >> slot & 1:
                current |= packet_value << (8 * slot)
        histograms[:, mask] = accumulator_histograms(current)
    return histograms


class RobustRepeatedKernel:
    def __init__(
        self,
        histograms: np.ndarray,
        split_caps,
        pole: float,
        fugacity: float,
        *,
        fixed_all_active: bool = False,
    ):
        self.histograms = histograms
        self.split_caps = split_caps
        self.pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
        if fixed_all_active:
            self.mask_weights = np.zeros(256, dtype=np.float64)
            self.mask_weights[-1] = 1.0
        else:
            packet_counts = np.array([mask.bit_count() for mask in range(256)])
            self.mask_weights = fugacity**packet_counts
        self.denominators = np.array(
            [float(math.comb(BITS, state)) for state in range(BITS + 1)]
        )

    def apply(self, values: np.ndarray) -> np.ndarray:
        moments = np.zeros((BITS + 1, 256), dtype=np.float64)
        for emitted in range(BITS + 1):
            counts = self.histograms[:, :, emitted].astype(np.float64)
            if not np.any(counts):
                continue
            order = np.argsort(values)[::-1]
            scores = values[order]
            capacities = np.array(
                [float(self.split_caps[emitted][state]) for state in order]
            )
            positive = capacities > 0.0
            capacities = capacities[positive]
            scores = scores[positive]
            cumulative = np.cumsum(capacities)
            weighted = np.cumsum(capacities * scores)
            flat = counts.ravel()
            indices = np.searchsorted(cumulative, flat, side="left")
            indices = np.minimum(indices, len(cumulative) - 1)
            prior_indices = np.maximum(indices - 1, 0)
            prior_mass = np.where(indices == 0, 0.0, cumulative[prior_indices])
            prior_weighted = np.where(indices == 0, 0.0, weighted[prior_indices])
            curve = prior_weighted + (flat - prior_mass) * scores[indices]
            moments += self.pole_powers[emitted] * curve.reshape(counts.shape)
        return (moments @ self.mask_weights) / self.denominators


def witness(kernel: RobustRepeatedKernel, iterations: int):
    values = np.ones(BITS + 1, dtype=np.float64)
    for _ in range(iterations):
        image = kernel.apply(values)
        scale = float(np.max(image))
        values = image / scale
    image = kernel.apply(values)
    ratios = image / values
    eigenvalue = float(np.max(ratios))
    domination = float(np.max(1.0 / values))
    return eigenvalue, domination, values, int(np.argmax(ratios))


def best_witness(kernel: RobustRepeatedKernel, iterations: int, blocks: int):
    """Return the strongest finite-block bound on one power trajectory.

    Every positive trajectory vector is a valid Collatz test vector.  The
    finite-block objective is not monotone in the iteration number, so the
    last vector need not be the strongest one encountered.
    """

    if iterations < 0 or blocks < 1:
        raise ValueError("best witness requires nonnegative iterations and blocks")
    values = np.ones(BITS + 1, dtype=np.float64)
    best = None
    for iteration in range(iterations + 1):
        image = kernel.apply(values)
        if np.any(values <= 0.0) or np.any(image <= 0.0):
            raise FloatingPointError("Collatz trajectory left the positive cone")
        ratios = image / values
        eigenvalue = float(np.max(ratios))
        domination = float(np.max(1.0 / values))
        score_log2 = (
            math.log2(domination)
            + blocks * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        candidate = (
            score_log2,
            iteration,
            eigenvalue,
            domination,
            values.copy(),
            int(np.argmax(ratios)),
        )
        if best is None or candidate[:2] < best[:2]:
            best = candidate
        if iteration != iterations:
            values = image / float(np.max(image))
    assert best is not None
    score_log2, iteration, eigenvalue, domination, values, worst_state = best
    return eigenvalue, domination, values, worst_state, iteration, score_log2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-value", type=lambda value: int(value, 0), default=0xFF)
    parser.add_argument("--active-numerator", type=int, default=21)
    parser.add_argument("--active-denominator", type=int, default=128)
    parser.add_argument("--output-pole", type=float, default=0.8)
    parser.add_argument("--fugacity", type=float, default=0.559040145)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument("--fixed-all-active", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.packet_value <= 0xFF:
        raise SystemExit("repeated state transfer: packet value must be nonzero")
    if not 0.0 < args.output_pole < 1.0 or args.fugacity <= 0.0:
        raise SystemExit("repeated state transfer: invalid pole or fugacity")
    packets = (
        PACKET_SLOTS
        if args.fixed_all_active
        else PACKET_SLOTS * args.active_numerator // args.active_denominator
    )
    if (
        not args.fixed_all_active
        and packets * args.active_denominator
        != PACKET_SLOTS * args.active_numerator
    ):
        raise SystemExit("repeated state transfer: nonintegral packet count")

    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    exact = load_exact(EXACT_SLICES)
    split_caps = build_split_caps(spectrum, exact)
    histograms = build_histograms(args.packet_value)
    kernel = RobustRepeatedKernel(
        histograms,
        split_caps,
        args.output_pole,
        args.fugacity,
        fixed_all_active=args.fixed_all_active,
    )
    eigenvalue, domination, values, worst_state = witness(kernel, args.iterations)
    mgf_log2 = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    if args.fixed_all_active:
        probability_log2 = mgf_log2 - D * math.log2(args.output_pole)
    else:
        probability_log2 = (
            mgf_log2
            - packets * math.log2(args.fugacity)
            - log2_binomial(PACKET_SLOTS, packets)
            - D * math.log2(args.output_pole)
        )
    print("repeated-packet 65-state robust transfer probe")
    print(
        f"packet_value=0x{args.packet_value:02x} active_packets={packets} "
        f"pole={args.output_pole} fugacity={args.fugacity} "
        f"fixed_all_active={args.fixed_all_active}"
    )
    print(f"lambda_log2={math.log2(eigenvalue):.12f}")
    print(f"domination_log2={math.log2(domination):.12f}")
    print(f"worst_state={worst_state}")
    print(f"conditional_probability_log2={probability_log2:.9f}")
    print("status=DIAGNOSTIC_EXACT_ACCUMULATOR_ROBUST_SPLIT_BINARY64")


if __name__ == "__main__":
    main()
