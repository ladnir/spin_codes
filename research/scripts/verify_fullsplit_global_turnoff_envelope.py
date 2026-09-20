#!/usr/bin/env python3
"""Verify a global one-step termination envelope for current full-split wrappers.

The tiny-prefix certificate uses a sharper finite-prefix exact-cancellation
atom.  The theorem-facing all-episode wrappers need a single atom that is safe
for every larger H currently covered by the late-window certificate.

For n=bT and a fixed state pattern of weight q, exact cancellation has

    C(n-b, H-q) / C(n, H)
      = (H)_q (n-H)_{b-q} / (n)_b.

With alpha=H/n and n>=b*T_min,

    (H)_q (n-H)_{b-q} / (n)_b
      <= correction * alpha^q (1-alpha)^{b-q},

where correction=(1-(b-1)/(b*T_min))^{-b}.  Averaging over Q gives a
Bernstein polynomial.  Since

    sum_q Pr[Q=q] alpha^q(1-alpha)^{b-q}
      = sum_q (Pr[Q=q]/C(b,q)) B_{q,b}(alpha),

it is enough to check that the largest nonzero Bernstein coefficient is the
q=b coefficient.  This script verifies that finite coefficient inequality and
prints the resulting global p_term upper bound.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from analyze_fullsplit_exact_cancellation import load_spectrum, split_q_law
from dense_largek_eval import log2add


ROOT = Path(__file__).resolve().parent


def check_close(name: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise SystemExit(
            f"{name}: expected {expected:.12f}, got {actual:.12f}; "
            f"delta={actual - expected:.6g}"
        )


def check_leq(name: str, actual: float, upper: float, tol: float) -> None:
    if actual > upper + tol:
        raise SystemExit(
            f"{name}: expected <= {upper:.12f}, got {actual:.12f}; "
            f"excess={actual - upper:.6g}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--t-min", type=int, default=5949)
    parser.add_argument("--expected-p0-log2", type=float, default=-64.00457432724919)
    parser.add_argument("--expected-eta-envelope-log2", type=float, default=-62.98700592417877)
    parser.add_argument("--expected-pterm-envelope-log2", type=float, default=-62.407875819071386)
    parser.add_argument(
        "--pterm-upper-log2",
        type=float,
        default=-62.4078758,
        help="Conservative rounded global atom for theorem-facing wrappers.",
    )
    parser.add_argument("--tolerance", type=float, default=5e-10)
    args = parser.parse_args()

    b = args.block_bits
    if args.t_min <= 1:
        raise SystemExit("--t-min must be larger than 1")

    pq = split_q_law(load_spectrum(args.spectrum), b)
    p0_log2 = math.log2(pq[0])
    q_full_log2 = math.log2(pq[b])

    max_coeff_log2 = float("-inf")
    max_coeff_q = -1
    for q in range(1, b + 1):
        if pq[q] <= 0.0:
            continue
        coeff = math.log2(pq[q]) - math.log2(math.comb(b, q))
        if coeff > max_coeff_log2:
            max_coeff_log2 = coeff
            max_coeff_q = q

    # (n)_b >= n^b * (1-(b-1)/n)^b, and n>=b*t_min.
    n_min = b * args.t_min
    correction_log2 = -b * math.log2(1.0 - (b - 1) / n_min)
    eta_envelope_log2 = q_full_log2 + correction_log2
    pterm_envelope_log2 = log2add(p0_log2, eta_envelope_log2)

    check_close("p0_log2", p0_log2, args.expected_p0_log2, args.tolerance)
    check_close(
        "eta_envelope_log2",
        eta_envelope_log2,
        args.expected_eta_envelope_log2,
        args.tolerance,
    )
    check_close(
        "pterm_envelope_log2",
        pterm_envelope_log2,
        args.expected_pterm_envelope_log2,
        args.tolerance,
    )
    check_leq("pterm_envelope_log2", pterm_envelope_log2, args.pterm_upper_log2, args.tolerance)
    if max_coeff_q != b:
        raise SystemExit(f"Bernstein coefficient maximum at q={max_coeff_q}, expected q={b}")
    check_leq("max_bernstein_coeff_log2", max_coeff_log2, q_full_log2, args.tolerance)

    print("Full-split global termination envelope verifier")
    print(f"spectrum,{args.spectrum}")
    print(f"block_bits,{b}")
    print(f"T_min,{args.t_min}")
    print(f"self_turnoff_p0_log2,{p0_log2:.14f}")
    print(f"full_state_q_log2,{q_full_log2:.14f}")
    print(f"max_bernstein_coeff_q,{max_coeff_q}")
    print(f"max_bernstein_coeff_log2,{max_coeff_log2:.14f}")
    print(f"finite_population_correction_log2,{correction_log2:.14f}")
    print(f"eta_envelope_log2,{eta_envelope_log2:.14f}")
    print(f"computed_pterm_envelope_log2,{pterm_envelope_log2:.14f}")
    print(f"certified_global_pterm_upper_log2,{args.pterm_upper_log2:.14f}")
    print("status,PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
