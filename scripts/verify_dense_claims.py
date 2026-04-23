#!/usr/bin/env python3
"""Reproduce the explicit dense+dense constant checks used in the paper.

This script evaluates the explicit asymptotic rate functions appearing in:
  - outerDense.tex
  - innerDense.tex
  - integration.tex

It is intentionally self-contained and uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass


ETA_CRIT = 1.0 - 2.0 ** -0.5
OUTER_LOW_SLOPE = math.log2(1.0 + math.sqrt(2.0))


def h2(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        if x == 0.0 or x == 1.0:
            return 0.0
        raise ValueError(f"x must lie in [0,1], got {x}")
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def outer_exponent(eta: float) -> float:
    if eta <= 0.0 or eta >= 1.0:
        raise ValueError(f"eta must lie in (0,1), got {eta}")
    if eta <= ETA_CRIT:
        return OUTER_LOW_SLOPE * eta
    return h2(eta) - 0.5


def psi_run(eta: float, theta: float) -> float:
    if not (0.0 < eta < 1.0):
        raise ValueError(f"eta must lie in (0,1), got {eta}")
    if not (0.0 < theta < 1.0 - eta):
        raise ValueError(f"theta must lie in (0,1-eta), got theta={theta}, eta={eta}")
    alpha = theta * eta / (1.0 - eta)
    if not (0.0 < alpha < 1.0):
        raise ValueError(f"alpha out of range: {alpha}")
    return h2(eta) - eta * h2(theta) - (1.0 - eta) * h2(alpha)


def e_bin(delta: float, xi: float) -> float:
    if not (0.0 < delta < 0.5):
        raise ValueError(f"delta must lie in (0,1/2), got {delta}")
    if xi <= 0.0:
        raise ValueError(f"xi must be positive, got {xi}")
    q = 1.0 / (2.0 + xi)
    return (2.0 + xi) * delta * (1.0 - h2(q))


def inner_exponent(eta: float, theta: float, delta: float, xi: float) -> float:
    return min(psi_run(eta, theta), e_bin(delta, xi))


@dataclass
class GapResult:
    worst_gap: float
    worst_eta: float
    outer_at_worst: float
    inner_at_worst: float
    bin_at_worst: float
    run_at_worst: float
    samples: int


def scan_gap(
    *,
    eta_lo: float,
    eta_hi: float,
    theta: float,
    delta: float,
    xi: float,
    step: float,
) -> GapResult:
    worst_gap = -float("inf")
    worst_eta = None
    outer_val = inner_val = bin_val = run_val = None
    samples = 0

    i = 0
    while True:
        eta = eta_lo + i * step
        if eta > eta_hi + 1e-18:
            break
        if eta >= 1.0 - theta:
            break
        out = outer_exponent(eta)
        run = psi_run(eta, theta)
        bin_term = e_bin(delta, xi)
        inn = min(run, bin_term)
        gap = out - inn
        if gap > worst_gap:
            worst_gap = gap
            worst_eta = eta
            outer_val = out
            inner_val = inn
            bin_val = bin_term
            run_val = run
        i += 1
        samples += 1

    if worst_eta is None:
        raise RuntimeError("No sample points evaluated")

    return GapResult(
        worst_gap=worst_gap,
        worst_eta=worst_eta,
        outer_at_worst=outer_val,
        inner_at_worst=inner_val,
        bin_at_worst=bin_val,
        run_at_worst=run_val,
        samples=samples,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delta", type=float, default=0.12)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-lo", type=float, default=ETA_CRIT)
    parser.add_argument("--eta-hi", type=float, default=0.99)
    parser.add_argument("--step", type=float, default=1e-4)
    args = parser.parse_args()

    result = scan_gap(
        eta_lo=args.eta_lo,
        eta_hi=args.eta_hi,
        theta=args.theta,
        delta=args.delta,
        xi=args.xi,
        step=args.step,
    )

    print("Dense+dense explicit-constant verification")
    print(f"eta_crit = {ETA_CRIT:.12f}")
    print(f"outer low-weight slope = log2(1+sqrt(2)) = {OUTER_LOW_SLOPE:.12f}")
    print(f"delta = {args.delta}")
    print(f"theta = {args.theta}")
    print(f"xi = {args.xi}")
    print(f"E_bin(delta, xi) = {e_bin(args.delta, args.xi):.12f}")
    print(f"window = [{args.eta_lo:.12f}, {args.eta_hi:.12f}]")
    print(f"step = {args.step}")
    print(f"samples = {result.samples}")
    print()
    print(f"worst gap = {result.worst_gap:.12f}")
    print(f"worst eta = {result.worst_eta:.12f}")
    print(f"outer(eta*) = {result.outer_at_worst:.12f}")
    print(f"run(eta*)   = {result.run_at_worst:.12f}")
    print(f"inner(eta*) = {result.inner_at_worst:.12f}")
    print(f"bin(eta*)   = {result.bin_at_worst:.12f}")
    print()
    if result.worst_gap < 0.0:
        print("STATUS: VERIFIED ON THE SAMPLED GRID (gap stays negative).")
        return 0

    print("STATUS: FAILED (gap became nonnegative).")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
