#!/usr/bin/env python3
"""Nonlinear Perron probe for the systematic striped group-chain inner.

For a fixed outer word, mark a physical 64-coordinate group active iff it is
nonzero.  A uniform group permutation makes its g active positions a uniform
g-subset of the B chain positions.  Local lane and incoming-state
permutations make their supports uniform conditional on weight.

This probe builds two robust 65-state operators:

* Z: one zero-input group;
* A: one nonzero-input group, maximized over its unknown weight 1..64.

The systematic split kernel uses exact C[t,q] rows where available and the
rigorous ordinary-spectrum caps elsewhere.  For positive x,

    H_x(v) = Z(v) + x A(v)

upper-bounds the active-position generating function.  A positive vector v
with H_x(v)<=lambda*v gives a nonlinear Perron certificate for B steps.
Dividing the x-Cauchy coefficient bound by C(B,g) yields the MGF averaged over
a uniform g-subset.  This is currently a long-double diagnostic; the local
row caps and combinatorics are rigorous, but outward arithmetic remains to be
added after the pole schedule is fixed.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import B, EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum


CHAIN_BLOCKS = 32768
TOTAL_LENGTH = 1 << 21
DISTANCE = 9 * TOTAL_LENGTH // 100


def log2_binomial(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2)


def load_exact(path: Path) -> dict[int, Counter[int]]:
    result: dict[int, Counter[int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            t = int(row["input_weight"])
            q = int(row["state_weight"])
            result.setdefault(t, Counter())[q] += int(row["count"])
    for t, histogram in result.items():
        if sum(histogram.values()) != math.comb(B, t):
            raise SystemExit(f"group operator: exact row t={t} has wrong mass")
    return result


def build_xor_probabilities() -> np.ndarray:
    probabilities = np.zeros((B + 1, B, B + 1), dtype=np.longdouble)
    for state_weight in range(B + 1):
        for input_weight in range(1, B + 1):
            denominator = math.comb(B, input_weight)
            for intersection in range(
                max(0, state_weight + input_weight - B),
                min(state_weight, input_weight) + 1,
            ):
                xor_weight = state_weight + input_weight - 2 * intersection
                probabilities[state_weight, input_weight - 1, xor_weight] = (
                    np.longdouble(math.comb(state_weight, intersection))
                    * math.comb(B - state_weight, input_weight - intersection)
                    / denominator
                )
    if float(np.max(np.abs(np.sum(probabilities, axis=2) - 1))) > 1e-15:
        raise SystemExit("group operator: XOR rows do not sum to one")
    return probabilities


class RobustKernel:
    def __init__(
        self,
        *,
        output_pole: float,
        spectrum: dict[int, int],
        exact: dict[int, Counter[int]],
    ) -> None:
        self.output_powers = np.array(
            [np.longdouble(output_pole) ** t for t in range(B + 1)],
            dtype=np.longdouble,
        )
        self.spectrum = spectrum
        self.exact = exact
        self.populations = [math.comb(B, t) for t in range(B + 1)]
        self.xor_probabilities = build_xor_probabilities()

    def split(self, value: np.ndarray) -> np.ndarray:
        result = np.zeros(B + 1, dtype=np.longdouble)
        result[0] = value[0]
        order = np.argsort(-value)
        for t in range(1, B + 1):
            if t in self.exact:
                result[t] = self.output_powers[t] * sum(
                    np.longdouble(count) * value[q]
                    for q, count in self.exact[t].items()
                ) / self.populations[t]
                continue
            remaining = self.populations[t]
            total = np.longdouble(0)
            for q in order:
                if q == 0:
                    continue
                take = min(
                    remaining,
                    self.spectrum.get(t + int(q), 0),
                    self.populations[int(q)],
                )
                total += np.longdouble(take) * value[q]
                remaining -= take
                if not remaining:
                    break
            if remaining:
                raise SystemExit(f"group operator: caps do not cover t={t}")
            result[t] = self.output_powers[t] * total / self.populations[t]
        return result

    def apply(self, value: np.ndarray, active_pole: float) -> np.ndarray:
        split = self.split(value)
        zero = split.copy()
        # For each incoming state weight, maximize over the nonzero input
        # group weight.  Matrix multiplication averages over the exact
        # hypergeometric XOR-weight law.
        active_by_weight = self.xor_probabilities @ split
        active = np.max(active_by_weight, axis=1)
        return zero + np.longdouble(active_pole) * active


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-pole", type=float, default=0.99)
    parser.add_argument("--active-pole", type=float, required=True)
    parser.add_argument("--active-groups", type=str, default="22,24,32,44,64,128,256,500")
    parser.add_argument("--iterations", type=int, default=160)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    args = parser.parse_args()
    if not 0 < args.output_pole < 1:
        raise SystemExit("group operator: output pole must lie in (0,1)")
    if args.active_pole <= 0:
        raise SystemExit("group operator: active pole must be positive")
    groups = sorted({int(value) for value in args.active_groups.split(",")})
    if not groups or groups[0] < 1 or groups[-1] > CHAIN_BLOCKS:
        raise SystemExit("group operator: invalid active group list")

    loaded = load_ebch128_spectrum(args.spectrum)
    spectrum = {weight: count for weight, count in loaded if weight and count}
    exact = load_exact(args.exact_slices)
    if any(
        count > spectrum.get(t + q, 0)
        for t, histogram in exact.items()
        for q, count in histogram.items()
    ):
        raise SystemExit("group operator: exact split row exceeds BCH spectrum")
    kernel = RobustKernel(
        output_pole=args.output_pole, spectrum=spectrum, exact=exact
    )

    value = np.ones(B + 1, dtype=np.longdouble)
    for _ in range(args.iterations):
        next_value = kernel.apply(value, args.active_pole)
        next_value /= np.max(next_value)
        value = next_value
    image = kernel.apply(value, args.active_pole)
    ratios = image / value
    lam = float(np.max(ratios))
    domination = float(np.max(1 / value))
    sequence_log2 = (
        math.log2(domination)
        + CHAIN_BLOCKS * math.log2(lam)
        + math.log2(float(value[0]))
    )

    print("systematic striped group-chain nonlinear Perron probe")
    print(
        f"B={CHAIN_BLOCKS} d={DISTANCE} output_pole={args.output_pole} "
        f"active_pole={args.active_pole} iterations={args.iterations}"
    )
    print(f"lambda_log2={math.log2(lam):.12f}")
    print(f"domination_log2={math.log2(domination):.12f}")
    print(f"sequence_generating_log2_upper={sequence_log2:.12f}")
    for g in groups:
        mgf = (
            sequence_log2
            - g * math.log2(args.active_pole)
            - log2_binomial(CHAIN_BLOCKS, g)
        )
        probability = mgf - DISTANCE * math.log2(args.output_pole)
        print(
            f"g={g} averaged_mgf_log2_upper={mgf:.12f} "
            f"low_output_probability_log2_upper={probability:.12f}"
        )
    print("status=DIAGNOSTIC_ONLY_NEEDS_EXACT_OR_OUTWARD_ARITHMETIC")


if __name__ == "__main__":
    main()
