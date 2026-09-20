#!/usr/bin/env python3
"""Asymptotic rate calculations for binary wrapped Expand--Convolute."""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar


def binary_entropy(x: float) -> float:
    """Return binary entropy in nats."""
    if not 0.0 <= x <= 1.0:
        raise ValueError("entropy argument must lie in [0,1]")
    if x == 0.0 or x == 1.0:
        return 0.0
    return -x * math.log(x) - (1.0 - x) * math.log1p(-x)


def binary_gv_distance(rate: float) -> float:
    """Return the relative distance solving h(delta)=(1-rate) ln 2."""
    if not 0.0 < rate < 1.0:
        raise ValueError("rate must lie in (0,1)")
    target = (1.0 - rate) * math.log(2.0)
    return float(brentq(lambda delta: binary_entropy(delta) - target, 0.0, 0.5))


def sparse_rate(memory: int, relative_distance: float) -> float:
    """Return the wrapped-EC sparse exponent I^w_{memory,delta}."""
    if memory < 1:
        raise ValueError("memory must be positive")
    if not 0.0 < relative_distance < 0.5:
        raise ValueError("relative_distance must lie in (0,1/2)")
    correction = math.ldexp(relative_distance, 1 - memory)
    return (
        math.sqrt(1.0 - 2.0 * relative_distance + correction)
        - math.sqrt(correction)
    ) ** 2


def effective_eigenvalue(memory: int, theta: float) -> float:
    """Return the dominant eigenvalue of the sparse effective generator."""
    if memory < 1:
        raise ValueError("memory must be positive")
    if theta < 0.0:
        raise ValueError("theta must be nonnegative")
    denominator = float((1 << memory) - 1)
    return (
        -(1.0 + 1.0 / denominator + theta * (1 << (memory - 1)) / denominator)
        + math.sqrt(
            (
                1.0
                - 1.0 / denominator
                - theta * (1 << (memory - 1)) / denominator
            )
            ** 2
            + 4.0 / denominator
        )
    ) / 2.0


def optimized_sparse_rate(memory: int, relative_distance: float) -> float:
    """Numerically optimize the two-class sparse exponent."""
    result = minimize_scalar(
        lambda theta: relative_distance * theta
        + effective_eigenvalue(memory, theta),
        bounds=(0.0, 64.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    if not result.success:
        raise RuntimeError("sparse optimization failed")
    return -float(result.fun)


def wrapping_transition_matrix(q: float, z: float, memory: int) -> np.ndarray:
    """Return the Bernoulli-input wrapping transfer matrix N_{memory,q}(z)."""
    if memory < 1:
        raise ValueError("memory must be positive")
    if not 0.0 <= q <= 0.5:
        raise ValueError("q must lie in [0,1/2]")
    if not 0.0 < z <= 1.0:
        raise ValueError("z must lie in (0,1]")
    matrix = np.zeros((memory + 1, memory + 1), dtype=np.float64)
    for state in range(memory - 1):
        matrix[state, 0] = z / 2.0
        matrix[state, state + 1] = 0.5
    matrix[memory - 1, 0] = (1.0 - q) * z
    matrix[memory - 1, memory] = q
    matrix[memory, 0] = q * z
    matrix[memory, memory] = 1.0 - q
    return matrix


def spectral_radius(q: float, z: float, memory: int) -> float:
    """Return the Perron root of the nonnegative wrapping transfer matrix."""
    eigenvalues = np.linalg.eigvals(wrapping_transition_matrix(q, z, memory))
    return float(np.max(np.abs(eigenvalues)))


def tail_rate(q: float, memory: int, relative_distance: float) -> float:
    """Numerically evaluate J^w_{memory,delta}(q)."""
    if not 0.0 < q <= 0.5:
        raise ValueError("q must lie in (0,1/2]")
    result = minimize_scalar(
        lambda log_z: (
            math.log(spectral_radius(q, math.exp(log_z), memory))
            - relative_distance * log_z
        ),
        bounds=(-64.0, 0.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if not result.success:
        raise RuntimeError("tail-rate optimization failed")
    return -float(result.fun)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory", type=int, default=21)
    parser.add_argument("--relative-distance", type=float, default=0.05)
    parser.add_argument("--rate", type=float, default=0.2)
    args = parser.parse_args()

    exponent = sparse_rate(args.memory, args.relative_distance)
    gv_distance = binary_gv_distance(args.rate)
    print(f"sparse exponent: {exponent:.12f}")
    print(f"density threshold C: {1.0 / exponent:.12f}")
    print(f"binary GV distance: {gv_distance:.12f}")
    print(
        "GV margin: "
        f"{(1.0 - args.rate) * math.log(2.0) - binary_entropy(args.relative_distance):.12f}"
    )


if __name__ == "__main__":
    main()
