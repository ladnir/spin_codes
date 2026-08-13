#!/usr/bin/env python3
"""Floating pole scan for the robust systematic one-episode envelope.

The theorem certificate will use rational poles and exact/outward arithmetic.
This probe identifies whether several output poles materially improve the
short-episode interval that a single pole at 2333/2373 leaves trivial.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter

import numpy as np

from analyze_systematic_group_kernel import B, EXACT_SLICES, SPECTRUM
from analyze_nosinger_bch_kernel import accumulator_weight_entries
from certificate_spectra import load_ebch128_spectrum


N = 1 << 21
DISTANCE = 9 * N // 100


def load_rows() -> dict[int, Counter[int]]:
    rows: dict[int, Counter[int]] = {}
    with EXACT_SLICES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(int(row["input_weight"]), Counter())[
                int(row["state_weight"])
            ] += int(row["count"])
    return rows


def maximize_row(
    t: int,
    scores: np.ndarray,
    spectrum: dict[int, int],
    rows: dict[int, Counter[int]],
) -> float:
    population = math.comb(B, t)
    if t in rows:
        return sum(count * scores[q] for q, count in rows[t].items()) / population
    remaining = population
    total = 0.0
    for q in np.argsort(-scores):
        q = int(q)
        if q == 0:
            continue
        take = min(remaining, spectrum.get(t + q, 0), math.comb(B, q))
        total += take * scores[q]
        remaining -= take
        if not remaining:
            return total / population
    raise RuntimeError(f"caps do not cover row {t}")


def envelope(
    pole: float,
    spectrum: dict[int, int],
    rows: dict[int, Counter[int]],
    burnin: int,
    accumulate: bool,
) -> tuple[float, float]:
    values = np.ones(B + 1)
    maxima: list[float] = []
    powers = np.array([pole**t for t in range(B + 1)])
    mixer = np.zeros((B + 1, B + 1))
    mixer[0, 0] = 1.0
    for source in range(1, B + 1):
        if accumulate:
            denominator = math.comb(B, source)
            for target, count in accumulator_weight_entries(source):
                mixer[source, target] = count / denominator
        else:
            mixer[source, source] = 1.0

    def apply(values: np.ndarray) -> np.ndarray:
        split = np.zeros(B + 1)
        for t in range(1, B + 1):
            split[t] = powers[t] * maximize_row(t, values, spectrum, rows)
        return mixer @ split

    for _ in range(burnin):
        image = apply(values)
        values = image
        maxima.append(float(np.max(values)))
    image = apply(values)
    rho = float(np.max(image[1:] / values[1:]))
    prefactor = max(value / rho**step for step, value in enumerate(maxima))
    return rho, prefactor


def geometric_envelope(
    pole: float,
    spectrum: dict[int, int],
    rows: dict[int, Counter[int]],
    accumulate: bool,
    iterations: int = 160,
) -> tuple[float, float]:
    """Floating Collatz envelope from a normalized positive test vector."""

    powers = np.array([pole**t for t in range(B + 1)])
    mixer = np.zeros((B + 1, B + 1))
    for source in range(1, B + 1):
        if accumulate:
            denominator = math.comb(B, source)
            for target, count in accumulator_weight_entries(source):
                mixer[source, target] = count / denominator
        else:
            mixer[source, source] = 1.0

    def apply(values: np.ndarray) -> np.ndarray:
        split = np.zeros(B + 1)
        for t in range(1, B + 1):
            split[t] = powers[t] * maximize_row(t, values, spectrum, rows)
        return mixer @ split

    values = np.ones(B + 1)
    values[0] = 0.0
    for _ in range(iterations):
        values = apply(values)
        values /= np.max(values)
    image = apply(values)
    rho = float(np.max(image[1:] / values[1:]))
    prefactor = float(np.max(1 / values[1:]) * np.max(image))
    return rho, prefactor


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--poles",
        default="0.50,0.60,0.70,0.80,0.86,0.90,0.93,0.95,0.97,0.98,0.9831437008,0.99,0.995",
    )
    parser.add_argument(
        "--geometric",
        action="store_true",
        help="also report a generic normalized Collatz-vector envelope",
    )
    parser.add_argument("--burnin", type=int, default=12)
    parser.add_argument(
        "--accumulate",
        action="store_true",
        help="insert one random-permute-plus-accumulate round before E_sys",
    )
    args = parser.parse_args()
    poles = [float(value) for value in args.poles.split(",")]
    if any(not 0 < pole < 1 for pole in poles):
        raise SystemExit("multipole probe: poles must lie in (0,1)")
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    rows = load_rows()
    lengths = (2950, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7419)
    results = []
    for pole in poles:
        rho, prefactor = envelope(
            pole, spectrum, rows, args.burnin, args.accumulate
        )
        bounds = [
            min(
                0.0,
                -DISTANCE * math.log2(pole)
                + math.log2(prefactor)
                + (length - 1) * math.log2(rho),
            )
            for length in lengths
        ]
        results.append((pole, rho, prefactor, bounds))
        print(
            f"mixer={'PA1' if args.accumulate else 'identity'} "
            f"pole={pole:.10f} rho_log2={math.log2(rho):.12f} "
            f"prefactor_log2={math.log2(prefactor):.12f} "
            f"bounds=" + ",".join(f"{value:.6f}" for value in bounds)
        )
        if args.geometric:
            geometric_rho, geometric_prefactor = geometric_envelope(
                pole, spectrum, rows, args.accumulate
            )
            print(
                f"geometric_rho_log2={math.log2(geometric_rho):.12f} "
                f"geometric_prefactor_log2={math.log2(geometric_prefactor):.12f}"
            )
    best = [max(-math.inf, min(row[3][index] for row in results)) for index in range(len(lengths))]
    print("lengths=" + ",".join(map(str, lengths)))
    print("multipole_best=" + ",".join(f"{value:.6f}" for value in best))
    print("status=DIAGNOSTIC_ONLY_NEEDS_RATIONAL_OR_OUTWARD_CERTIFICATES")


if __name__ == "__main__":
    main()
