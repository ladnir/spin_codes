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
    length_frac = min(1.0, (2.0 + xi) * delta)
    threshold_frac = delta / length_frac
    if threshold_frac >= 0.5:
        return 0.0
    return length_frac * (1.0 - h2(threshold_frac))


def inner_exponent(eta: float, theta: float, delta: float, xi: float) -> float:
    return min(psi_run(eta, theta), e_bin(delta, xi))


def first_start_binomial_exponent(eta: float, delta: float, off_budget: float = 0.0) -> tuple[float, float]:
    """Integrated first-start placement plus suffix-binomial exponent.

    For a first one near time alpha*n, the exact hypergeometric placement cost
    contributes H2(eta) - (1-alpha) H2(eta/(1-alpha)). If at most
    off_budget*n later times are OFF, then the available ON-time suffix is
    1-alpha-off_budget, contributing the binomial tail exponent
    (1-alpha-off_budget) * (1 - H2(delta/(1-alpha-off_budget))).
    The returned value is the minimum over alpha. This is a diagnostic exponent,
    not the full inner theorem.
    """
    best = float("inf")
    best_alpha = 0.0
    def value(alpha: float) -> float:
        placement_suffix = 1.0 - alpha
        on_suffix = placement_suffix - off_budget
        if on_suffix <= 0.0:
            return float("inf")
        x = delta / on_suffix
        if x >= 0.5:
            return float("inf")
        if eta > placement_suffix:
            return float("inf")
        placement = h2(eta) - placement_suffix * h2(eta / placement_suffix)
        bin_exp = on_suffix * (1.0 - h2(x))
        return placement + bin_exp

    lo = 0.0
    hi = min(max(0.0, 1.0 - off_budget - 2.0 * delta), max(0.0, 1.0 - eta))
    if hi <= 0.0:
        return value(0.0), 0.0

    golden = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - golden * (hi - lo)
    d = lo + golden * (hi - lo)
    fc = value(c)
    fd = value(d)
    for _ in range(80):
        if fc > fd:
            lo = c
            c = d
            fc = fd
            d = lo + golden * (hi - lo)
            fd = value(d)
        else:
            hi = d
            d = c
            fd = fc
            c = hi - golden * (hi - lo)
            fc = value(c)
    best_alpha = (lo + hi) / 2.0
    best = value(best_alpha)
    return best, best_alpha


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
    parser.add_argument("--scan-delta", action="store_true", help="Scan for the largest sampled delta with negative gap.")
    parser.add_argument(
        "--diagnostic-first-start",
        action="store_true",
        help="Use the diagnostic first-start-plus-suffix-binomial exponent instead of the current inner exponent.",
    )
    parser.add_argument(
        "--off-budget",
        type=float,
        default=0.0,
        help="Reserve this fraction of the post-start suffix for OFF/restart time in the first-start diagnostic.",
    )
    parser.add_argument("--delta-min", type=float, default=0.05)
    parser.add_argument("--delta-max", type=float, default=0.13)
    parser.add_argument("--delta-step", type=float, default=0.0005)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-lo", type=float, default=ETA_CRIT)
    parser.add_argument("--eta-hi", type=float, default=0.99)
    parser.add_argument("--step", type=float, default=1e-4)
    args = parser.parse_args()

    if args.scan_delta:
        best = None
        first_fail = None
        i = 0
        while True:
            delta = args.delta_min + i * args.delta_step
            if delta > args.delta_max + 1e-18:
                break
            result = scan_gap(
                eta_lo=args.eta_lo,
                eta_hi=args.eta_hi,
                theta=args.theta,
                delta=delta,
                xi=args.xi,
                step=args.step,
            )
            if args.diagnostic_first_start:
                worst = -float("inf")
                worst_eta = args.eta_lo
                worst_alpha = 0.0
                eta_i = 0
                while True:
                    eta = args.eta_lo + eta_i * args.step
                    if eta > args.eta_hi + 1e-18 or eta >= 1.0 - args.theta:
                        break
                    inn, alpha = first_start_binomial_exponent(eta, delta, args.off_budget)
                    gap = outer_exponent(eta) - inn
                    if gap > worst:
                        worst = gap
                        worst_eta = eta
                        worst_alpha = alpha
                    eta_i += 1
                result = GapResult(worst, worst_eta, 0.0, 0.0, 0.0, 0.0, 0)
            ok = result.worst_gap < 0.0
            if ok:
                best = (delta, result)
            elif best is not None and first_fail is None:
                first_fail = (delta, result)
                break
            i += 1

        print("Dense+dense delta scan")
        print(f"theta = {args.theta}")
        print(f"xi = {args.xi}")
        print(f"mode = {'first-start diagnostic' if args.diagnostic_first_start else 'current inner exponent'}")
        if args.diagnostic_first_start:
            print(f"off budget = {args.off_budget}")
        print(f"window = [{args.eta_lo:.12f}, {args.eta_hi:.12f}]")
        print(f"eta step = {args.step}")
        print(f"delta grid = [{args.delta_min}, {args.delta_max}] step {args.delta_step}")
        if best is None:
            print("no passing delta found")
            return 1
        delta, result = best
        print()
        print(f"largest passing grid delta = {delta:.6f}")
        print(f"  worst gap = {result.worst_gap:.12f}")
        print(f"  worst eta = {result.worst_eta:.12f}")
        print(f"  E_bin = {e_bin(delta, args.xi):.12f}")
        if first_fail is not None:
            delta_fail, result_fail = first_fail
            print(f"first failing grid delta = {delta_fail:.6f}")
            print(f"  worst gap = {result_fail.worst_gap:.12f}")
            print(f"  worst eta = {result_fail.worst_eta:.12f}")
            print(f"  E_bin = {e_bin(delta_fail, args.xi):.12f}")
        return 0

    result = scan_gap(
        eta_lo=args.eta_lo,
        eta_hi=args.eta_hi,
        theta=args.theta,
        delta=args.delta,
        xi=args.xi,
        step=args.step,
    )
    if args.diagnostic_first_start:
        worst = -float("inf")
        worst_eta = args.eta_lo
        worst_alpha = 0.0
        i = 0
        while True:
            eta = args.eta_lo + i * args.step
            if eta > args.eta_hi + 1e-18 or eta >= 1.0 - args.theta:
                break
            inn, alpha = first_start_binomial_exponent(eta, args.delta, args.off_budget)
            gap = outer_exponent(eta) - inn
            if gap > worst:
                worst = gap
                worst_eta = eta
                worst_alpha = alpha
                result = GapResult(
                    worst_gap=gap,
                    worst_eta=eta,
                    outer_at_worst=outer_exponent(eta),
                    inner_at_worst=inn,
                    bin_at_worst=inn,
                    run_at_worst=0.0,
                    samples=0,
                )
            i += 1

    print("Dense+dense explicit-constant verification")
    print(f"eta_crit = {ETA_CRIT:.12f}")
    print(f"outer low-weight slope = log2(1+sqrt(2)) = {OUTER_LOW_SLOPE:.12f}")
    print(f"delta = {args.delta}")
    print(f"theta = {args.theta}")
    print(f"xi = {args.xi}")
    print(f"mode = {'first-start diagnostic' if args.diagnostic_first_start else 'current inner exponent'}")
    if args.diagnostic_first_start:
        print(f"off budget = {args.off_budget}")
    print(f"episode length fraction = {min(1.0, (2.0 + args.xi) * args.delta):.12f}")
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
