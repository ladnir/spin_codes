#!/usr/bin/env python3
"""Bulk exponent diagnostic for full-length designed BA-3.

For normalized weights x, y, and z after the constituent, first accumulator,
and second accumulator, respectively, the expected-spectrum exponent is

    a_C(x) + p_acc(x,y) + p_acc(y,z).

The constituent exponent a_C is obtained from the saddle point of W_C(r)^m.
The accumulator exponent is the Stirling limit of its exact IOWE.  This script
is a bulk diagnostic, not a finite-N interval certificate; the sparse finite-N
evaluator covers the constant-weight boundary separately.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, differential_evolution
from scipy.special import logsumexp


LN2 = math.log(2.0)


def read_spectrum(path: Path, local_length: int, local_dimension: int) -> np.ndarray:
    values = np.zeros(local_length + 1, dtype=np.float64)
    integer_sum = 0
    zero = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            weight = int(row["weight"])
            count = int(row["count"])
            values[weight] += count
            integer_sum += count
            if weight == 0:
                zero += count
    if zero != 1 or integer_sum != 1 << local_dimension:
        raise ValueError("invalid constituent spectrum")
    return values


def binary_entropy(value: float) -> float:
    if value <= 0.0 or value >= 1.0:
        return 0.0 if value in (0.0, 1.0) else -math.inf
    return -value * math.log2(value) - (1.0 - value) * math.log2(1.0 - value)


class ConstituentExponent:
    def __init__(self, spectrum: np.ndarray):
        self.spectrum = spectrum
        self.length = len(spectrum) - 1
        self.weights = np.arange(self.length + 1, dtype=np.float64)
        self.support = spectrum > 0
        self.log_counts = np.full(spectrum.shape, -np.inf)
        self.log_counts[self.support] = np.log(spectrum[self.support])
        self.minimum = int(np.flatnonzero(self.support)[0])
        self.maximum = int(np.flatnonzero(self.support)[-1])

    def mean_at_log2_r(self, log2_r: float) -> float:
        logs = self.log_counts + self.weights * log2_r * LN2
        normalizer = logsumexp(logs[self.support])
        probabilities = np.zeros(self.spectrum.shape)
        probabilities[self.support] = np.exp(logs[self.support] - normalizer)
        return float(np.dot(self.weights, probabilities))

    def __call__(self, normalized_weight: float) -> tuple[float, float]:
        target_mean = self.length * normalized_weight
        if target_mean <= self.minimum:
            if target_mean == 0.0 and self.spectrum[0] == 1:
                return 0.0, -math.inf
        if target_mean < 0.0 or target_mean > self.maximum:
            return -math.inf, math.nan

        def equation(log2_r: float) -> float:
            return self.mean_at_log2_r(log2_r) - target_mean

        log2_r = brentq(equation, -80.0, 80.0)
        logs = self.log_counts + self.weights * log2_r * LN2
        log2_w = float(logsumexp(logs[self.support]) / LN2)
        exponent = log2_w / self.length - normalized_weight * log2_r
        return exponent, log2_r


def accumulator_exponent(input_weight: float, output_weight: float) -> float:
    if input_weight <= 0.0 or input_weight >= 1.0:
        return -math.inf
    half = input_weight / 2.0
    if output_weight < half or output_weight > 1.0 - half:
        return -math.inf
    left_denominator = 1.0 - output_weight
    right_denominator = output_weight
    if left_denominator <= 0.0 or right_denominator <= 0.0:
        return -math.inf
    left = half / left_denominator
    right = half / right_denominator
    if not (0.0 <= left <= 1.0 and 0.0 <= right <= 1.0):
        return -math.inf
    return (
        left_denominator * binary_entropy(left)
        + right_denominator * binary_entropy(right)
        - binary_entropy(input_weight)
    )


def optimize(
    constituent: ConstituentExponent,
    relative_distance: float,
    x_min: float,
    seed: int,
) -> dict[str, float]:
    x_max = min(1.0, 4.0 * relative_distance)
    if x_min >= x_max:
        raise ValueError("x_min leaves an empty bulk domain")

    def decode(point: np.ndarray) -> tuple[float, float, float]:
        x, y_fraction, z_fraction = point
        y_low = x / 2.0
        y_high = 2.0 * relative_distance
        y = y_low + y_fraction * (y_high - y_low)
        z_low = y / 2.0
        z = z_low + z_fraction * (relative_distance - z_low)
        return float(x), float(y), float(z)

    def objective(point: np.ndarray) -> float:
        x, y, z = decode(point)
        outer, _ = constituent(x)
        first = accumulator_exponent(x, y)
        second = accumulator_exponent(y, z)
        value = outer + first + second
        return 1e6 if not math.isfinite(value) else -value

    result = differential_evolution(
        objective,
        bounds=((x_min, x_max), (0.0, 1.0), (0.0, 1.0)),
        seed=seed,
        tol=1e-11,
        polish=True,
        updating="immediate",
        workers=1,
    )
    x, y, z = decode(result.x)
    outer, log2_r = constituent(x)
    first = accumulator_exponent(x, y)
    second = accumulator_exponent(y, z)
    return {
        "max_exponent_per_output_bit": outer + first + second,
        "outer_relative_weight": x,
        "first_accumulator_relative_weight": y,
        "final_relative_weight": z,
        "constituent_exponent": outer,
        "first_accumulator_exponent": first,
        "second_accumulator_exponent": second,
        "constituent_saddle_log2_r": log2_r,
        "optimizer_success": bool(result.success),
        "optimizer_fun": float(result.fun),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--code-name", required=True)
    parser.add_argument("--local-length", type=int, required=True)
    parser.add_argument("--local-dimension", type=int, required=True)
    parser.add_argument("--length", type=int, default=1 << 21)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--x-min", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    spectrum = read_spectrum(args.spectrum, args.local_length, args.local_dimension)
    constituent = ConstituentExponent(spectrum)
    optimum = optimize(constituent, args.relative_distance, args.x_min, args.seed)
    result = {
        "schema": "full-designed-ba3-bulk-exponent-v1",
        "construction": {
            "name": f"DBA-3({args.code_name})",
            "length": args.length,
            "dimension": args.length // 2,
            "local_length": args.local_length,
            "local_dimension": args.local_dimension,
        },
        "domain": {
            "relative_distance": args.relative_distance,
            "outer_relative_weight_min": args.x_min,
            "outer_relative_weight_max": min(1.0, 4.0 * args.relative_distance),
        },
        "method": {
            "constituent": "coefficient saddle point of W_C(r)^(N/local_length)",
            "accumulators": "Stirling exponent of the exact accumulator IOWE",
            "status": "bulk diagnostic; finite-N interval certification remains separate",
        },
        "result": {
            **optimum,
            "scaled_exponent_at_length": optimum["max_exponent_per_output_bit"] * args.length,
        },
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
