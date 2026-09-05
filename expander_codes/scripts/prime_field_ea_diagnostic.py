#!/usr/bin/env python3
"""Asymptotic diagnostics for labeled prime-field Expand--Accumulate.

This script evaluates only the proved Bernoulli asymptotic conditions.  It does
not claim a finite-length failure probability.
"""

from __future__ import annotations

import argparse
import math

from scipy.optimize import minimize_scalar


def binary_entropy(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log(x) - (1.0 - x) * math.log1p(-x)


def shell_entropy(prime: int, delta: float) -> float:
    return binary_entropy(delta) + delta * math.log(prime - 1)


def gv_distance(prime: int, rate: float) -> float:
    """Return the smaller p-ary GV root in relative symbol weight."""
    target = (1.0 - rate) * math.log(prime)
    lo = 0.0
    hi = (prime - 1.0) / prime
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if shell_entropy(prime, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sparse_rate(prime: int, delta: float) -> float:
    return (
        math.sqrt(1.0 - delta)
        - math.sqrt(delta / (prime - 1.0))
    ) ** 2


def regular_shell_distribution(prime: int, length: int, steps: int) -> list[float]:
    """Weight law of a labeled unit-vector walk in one regular region."""
    if length <= 0 or steps < 0:
        raise ValueError("length must be positive and steps nonnegative")
    s = prime - 1.0
    current = [0.0] * (length + 1)
    current[0] = 1.0
    for _ in range(steps):
        following = [0.0] * (length + 1)
        for weight, probability in enumerate(current):
            if probability == 0.0:
                continue
            if weight < length:
                following[weight + 1] += probability * (length - weight) / length
            if weight:
                following[weight] += probability * weight / length * (prime - 2) / s
                following[weight - 1] += probability * weight / length / s
        current = following
    return current


def spectral_radius(prime: int, nonzero_probability: float, z: float) -> float:
    s = prime - 1.0
    q = nonzero_probability
    trace = 1.0 - q + (1.0 - q / s) * z
    difference = 1.0 - q - (1.0 - q / s) * z
    discriminant = difference * difference + 4.0 * q * q * z / s
    return 0.5 * (trace + math.sqrt(discriminant))


def tail_rate(prime: int, delta: float, nonzero_probability: float) -> float:
    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return -(
            delta * log_z
            - math.log(spectral_radius(prime, nonzero_probability, z))
        )

    optimum = minimize_scalar(
        objective,
        bounds=(-max(80.0, math.log(prime) + 20.0), 0.0),
        method="bounded",
        options={"xatol": 1e-13, "maxiter": 300},
    )
    return -float(optimum.fun)


def report(prime: int, rate: float, length: int) -> str:
    delta = gv_distance(prime, rate)
    information = sparse_rate(prime, delta)
    c_min = 1.0 / information
    degree = c_min * math.log(length)
    return (
        f"p={prime:>6d}  GV={delta:.9f}  "
        f"I_sparse={information:.9f}  C_min={c_min:.6f}  "
        f"C_min*ln(n)={degree:.3f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rate", type=float, default=0.5)
    parser.add_argument("--length", type=int, default=2**21)
    parser.add_argument(
        "--primes", type=int, nargs="+", default=[2, 3, 5, 7, 17, 257]
    )
    args = parser.parse_args()
    if not 0.0 < args.rate < 1.0:
        raise ValueError("rate must lie in (0,1)")
    if args.length <= 1:
        raise ValueError("length must exceed one")
    for prime in args.primes:
        if prime < 2:
            raise ValueError("all primes must be at least two")
        print(report(prime, args.rate, args.length))


if __name__ == "__main__":
    main()
