#!/usr/bin/env python3
"""Search simple sufficient constants for the Random SPIN sparse bound.

This is a binary64 design-space search, not a proof certificate.  Pass a
reported rational candidate to certify_random_spin_linear_exponent.py for
outward-rounded verification.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math


def rational(value: str) -> Fraction:
    return Fraction(value)


def evaluate(
    delta: Fraction,
    eta: Fraction,
    chi: Fraction,
    gamma: Fraction,
) -> dict[str, object] | None:
    delta_f = float(delta)
    eta_f = float(eta)
    chi_f = float(chi)
    gamma_f = float(gamma)
    a = eta_f / (delta_f * (chi_f + 1.0))
    if chi_f * a >= 1.0:
        return None
    b = (1.0 + a) / 2.0
    memory_margin = gamma_f * math.log2(1.0 / b) - 1.0
    if memory_margin <= 0.0:
        return None
    epsilon = (
        chi_f * chi_f * a / ((chi_f + 1.0) * (1.0 - chi_f * a))
        + eta_f / (2.0 * (1.0 - eta_f))
    )
    rho = delta_f * (chi_f + 1.0) * math.exp(epsilon)
    if rho >= math.sqrt(2.0) - 1.0:
        return None
    decay = math.log2(math.sqrt(2.0) / (1.0 + rho))
    polynomial_exponent = 2.0 * gamma_f + 1.5
    minimum_block_constant = math.floor(polynomial_exponent / decay) + 1
    block_margin = minimum_block_constant * decay - polynomial_exponent
    return {
        "delta": str(delta),
        "sparse_cutoff": str(eta),
        "interior_parameter": str(chi),
        "memory_gamma": str(gamma),
        "minimum_integer_block_constant": minimum_block_constant,
        "rho": rho,
        "outer_decay_log2": decay,
        "memory_margin": memory_margin,
        "block_margin": block_margin,
    }


def search(delta: Fraction) -> dict[str, object]:
    candidates: list[dict[str, object]] = []
    for cutoff_denominator in range(1000, 5001, 250):
        eta = Fraction(1, cutoff_denominator)
        for chi_hundredths in range(101, 151):
            chi = Fraction(chi_hundredths, 100)
            a = float(eta / (delta * (chi + 1)))
            b = (1.0 + a) / 2.0
            threshold = 1.0 / math.log2(1.0 / b)
            gamma_hundredths = math.floor(100.0 * threshold) + 2
            gamma = Fraction(gamma_hundredths, 100)
            candidate = evaluate(delta, eta, chi, gamma)
            if candidate is not None:
                candidates.append(candidate)
    if not candidates:
        raise RuntimeError("the search grid contains no sufficient candidate")
    return min(
        candidates,
        key=lambda item: (
            int(item["minimum_integer_block_constant"]),
            -float(item["block_margin"]),
            -float(item["memory_margin"]),
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta", default="5501/50000")
    parser.add_argument("--sparse-cutoff", default="1/5000")
    parser.add_argument("--interior-parameter", default="101/100")
    parser.add_argument("--memory-gamma", default="51/50")
    args = parser.parse_args()

    delta = rational(args.delta)
    requested = evaluate(
        delta,
        rational(args.sparse_cutoff),
        rational(args.interior_parameter),
        rational(args.memory_gamma),
    )
    if requested is None:
        raise SystemExit("the requested constants do not satisfy the sparse conditions")
    payload = {
        "schema": "random-spin-asymptotic-constant-search-v1",
        "arithmetic": "binary64 diagnostic; certify the selected candidate separately",
        "requested_candidate": requested,
        "best_grid_candidate": search(delta),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
