#!/usr/bin/env python3
"""Packet-weight-aware robust transfer probe for eight-bit Riffle atoms.

This probe does not replace the two EBCH [128,64,22] local codes and does not
assume that eight shuffled packets form a uniform 64-bit Hamming slice.

For a fixed 64-bit current input U of weight u and an independently permuted
state of weight s, the drive W=U+S has an exact hypergeometric weight-shell
decomposition.  On a shell of drive weight t:

* the exact accumulator table caps the number of drives producing emitted
  systematic-half weight y;
* the exact/rigorously capped systematic EBCH split table caps the number of
  y-weight values producing next-state weight q.

Greedily filling the highest score cells under the minimum of those two caps
is a rigorous relaxation for every fixed U of weight u.  The current script
uses binary64 Perron iteration and iid input-weight laws only to diagnose
whether this stronger local operator has useful dense-regime margin.  A final
proof must replace the iid law by coefficient extraction for the exact global
packet profile and harden the selected Perron witnesses with outward or exact
arithmetic.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_nosinger_bch_kernel import accumulator_weight_entries
from analyze_systematic_group_kernel import B, EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum


CHAIN_BLOCKS = 32768
TOTAL_LENGTH = 1 << 21
DISTANCE = 9 * TOTAL_LENGTH // 100


def load_exact(path: Path) -> dict[int, Counter[int]]:
    result: dict[int, Counter[int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            result.setdefault(int(row["input_weight"]), Counter())[
                int(row["state_weight"])
            ] += int(row["count"])
    for weight, histogram in result.items():
        if sum(histogram.values()) != math.comb(B, weight):
            raise SystemExit(f"packet8 transfer: exact split row {weight} has wrong mass")
    return result


def build_split_caps(
    spectrum: dict[int, int], exact: dict[int, Counter[int]]
) -> tuple[tuple[int, ...], ...]:
    caps: list[tuple[int, ...]] = []
    for emitted in range(B + 1):
        row = [0] * (B + 1)
        if emitted == 0:
            row[0] = 1
        elif emitted in exact:
            for state, count in exact[emitted].items():
                row[state] = count
        else:
            for state in range(1, B + 1):
                row[state] = min(
                    spectrum.get(emitted + state, 0), math.comb(B, state)
                )
        if sum(row) < math.comb(B, emitted):
            raise SystemExit(f"packet8 transfer: split caps do not cover row {emitted}")
        caps.append(tuple(row))
    return tuple(caps)


def build_accumulator_counts() -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []
    for drive in range(B + 1):
        row = [0] * (B + 1)
        if drive == 0:
            row[0] = 1
        else:
            for emitted, count in accumulator_weight_entries(drive):
                row[emitted] = count
        if sum(row) != math.comb(B, drive):
            raise SystemExit(f"packet8 transfer: accumulator row {drive} has wrong mass")
        rows.append(tuple(row))
    return tuple(rows)


def build_xor_shells() -> tuple[tuple[tuple[tuple[int, int], ...], ...], ...]:
    all_rows = []
    for state in range(B + 1):
        input_rows = []
        denominator = math.comb(B, state)
        for current in range(B + 1):
            shells = []
            for intersection in range(
                max(0, state + current - B), min(state, current) + 1
            ):
                drive = state + current - 2 * intersection
                count = math.comb(current, intersection) * math.comb(
                    B - current, state - intersection
                )
                if count:
                    shells.append((drive, count))
            if sum(count for _drive, count in shells) != denominator:
                raise SystemExit("packet8 transfer: XOR shell mass mismatch")
            input_rows.append(tuple(shells))
        all_rows.append(tuple(input_rows))
    return tuple(all_rows)


class RobustPacketWeightKernel:
    def __init__(self, output_pole: float, spectrum, exact) -> None:
        self.output_powers = np.array(
            [output_pole**weight for weight in range(B + 1)], dtype=np.float64
        )
        self.split_caps = build_split_caps(spectrum, exact)
        self.accumulator = build_accumulator_counts()
        self.shells = build_xor_shells()
        self.denominators = np.array(
            [float(math.comb(B, state)) for state in range(B + 1)]
        )

    def _shell_curves(self, values: np.ndarray):
        scores = [
            (self.output_powers[emitted] * values[state], emitted, state)
            for emitted in range(B + 1)
            for state in range(B + 1)
        ]
        scores.sort(reverse=True)
        curves = []
        for drive in range(B + 1):
            masses = [0]
            weighted = [0.0]
            interval_scores = []
            for score, emitted, state in scores:
                cap = min(
                    self.accumulator[drive][emitted],
                    self.split_caps[emitted][state],
                )
                if cap:
                    masses.append(masses[-1] + cap)
                    weighted.append(weighted[-1] + cap * score)
                    interval_scores.append(score)
            if masses[-1] < math.comb(B, drive):
                raise SystemExit(f"packet8 transfer: joint caps do not cover drive {drive}")
            curves.append((masses, weighted, interval_scores))
        return curves

    @staticmethod
    def _curve_value(curve, mass: int) -> float:
        masses, weighted, interval_scores = curve
        index = bisect.bisect_left(masses, mass)
        if index == 0:
            return 0.0
        if masses[index] == mass:
            return weighted[index]
        previous = index - 1
        return weighted[previous] + (mass - masses[previous]) * interval_scores[previous]

    def apply(self, values: np.ndarray, input_law: np.ndarray) -> np.ndarray:
        curves = self._shell_curves(values)
        image = np.zeros(B + 1, dtype=np.float64)
        for state in range(B + 1):
            total = 0.0
            for current, probability in enumerate(input_law):
                if probability == 0.0:
                    continue
                row = 0.0
                for drive, shell_count in self.shells[state][current]:
                    row += self._curve_value(curves[drive], shell_count)
                total += probability * row
            image[state] = total / self.denominators[state]
        return image


def binomial_input_law(theta: float) -> np.ndarray:
    return np.array(
        [
            math.comb(B, weight) * theta**weight * (1.0 - theta) ** (B - weight)
            for weight in range(B + 1)
        ],
        dtype=np.float64,
    )


def fixed_packet_weight_law(active_density: float, packet_weight: int) -> np.ndarray:
    law = np.zeros(B + 1, dtype=np.float64)
    for active_packets in range(9):
        law[active_packets * packet_weight] = (
            math.comb(8, active_packets)
            * active_density**active_packets
            * (1.0 - active_density) ** (8 - active_packets)
        )
    return law


def binary_entropy(theta: float) -> float:
    if theta in (0.0, 1.0):
        return 0.0
    return -theta * math.log2(theta) - (1.0 - theta) * math.log2(1.0 - theta)


def run(kernel, law: np.ndarray, label: str, entropy_theta: float | None, iterations: int, output_pole: float) -> None:
    values = np.ones(B + 1, dtype=np.float64)
    for _ in range(iterations):
        image = kernel.apply(values, law)
        scale = float(np.max(image))
        if not math.isfinite(scale) or scale <= 0:
            raise SystemExit("packet8 transfer: invalid power-iteration scale")
        values = image / scale
    image = kernel.apply(values, law)
    ratios = image / values
    eigenvalue = float(np.max(ratios))
    domination = float(np.max(1.0 / values))
    mgf_log2 = (
        math.log2(domination)
        + CHAIN_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    probability_log2 = mgf_log2 - DISTANCE * math.log2(output_pole)
    random_rate_half_spectrum = (
        float("nan")
        if entropy_theta is None
        else TOTAL_LENGTH * (binary_entropy(entropy_theta) - 0.5)
    )
    sum_diagnostic = (
        float("nan")
        if entropy_theta is None
        else probability_log2 + max(0.0, random_rate_half_spectrum)
    )
    print(
        f"{label} lambda_log2={math.log2(eigenvalue):.12f} "
        f"domination_log2={math.log2(domination):.12f} "
        f"low_output_probability_log2={probability_log2:.6f} "
        f"random_rate_half_spectrum_log2={random_rate_half_spectrum:.6f} "
        f"sum_diagnostic={sum_diagnostic:.6f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-pole", type=float, default=0.997)
    parser.add_argument("--theta", default="0.5")
    parser.add_argument("--packet-active-density", type=float)
    parser.add_argument("--packet-weight", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=48)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    args = parser.parse_args()
    if not 0.0 < args.output_pole < 1.0:
        raise SystemExit("packet8 transfer: output pole must lie in (0,1)")
    thetas = [float(value) for value in args.theta.split(",")]
    if any(not 0.0 < theta < 1.0 for theta in thetas):
        raise SystemExit("packet8 transfer: theta must lie in (0,1)")
    if args.packet_active_density is not None and not 0.0 < args.packet_active_density < 1.0:
        raise SystemExit("packet8 transfer: packet active density must lie in (0,1)")
    if not 1 <= args.packet_weight <= 8:
        raise SystemExit("packet8 transfer: packet weight must lie in 1..8")
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(args.spectrum)
        if weight and count
    }
    exact = load_exact(args.exact_slices)
    kernel = RobustPacketWeightKernel(args.output_pole, spectrum, exact)
    print("eight-bit packet-weight robust transfer probe")
    print(
        f"inner=[128,64,22] outer=[128,64,22] B={CHAIN_BLOCKS} "
        f"d={DISTANCE} output_pole={args.output_pole} iterations={args.iterations}"
    )
    for theta in thetas:
        run(
            kernel,
            binomial_input_law(theta),
            f"theta={theta:.8f}",
            theta,
            args.iterations,
            args.output_pole,
        )
    if args.packet_active_density is not None:
        run(
            kernel,
            fixed_packet_weight_law(args.packet_active_density, args.packet_weight),
            (
                f"packet_active_density={args.packet_active_density:.8f} "
                f"packet_weight={args.packet_weight}"
            ),
            None,
            args.iterations,
            args.output_pole,
        )
    print("status=DIAGNOSTIC_IID_PROFILE_BINARY64_PERRON")


if __name__ == "__main__":
    main()
