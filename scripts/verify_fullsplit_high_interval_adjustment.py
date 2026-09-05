#!/usr/bin/env python3
"""Verify the high-interval p_term adjustment used by the finite ledger.

The raw high-interval rows were generated with the sharper finite-prefix
termination atom.  The theorem-facing all-episode wrapper now uses the global
Bernstein atom.  This helper verifies that the ledger's additive allowance
covers the switch.

For the e>=1 branch, the one-termination term increases by

    global_pterm_log2 - prefix_pterm_log2.

The fixed-pole e>=2 multiplier contributes an extra factor
log2(1+R).  Uniformly over the current high-window certificate,
H<=bT and T<=T_max, so with p=2^global_pterm_log2 and 1-M<=1,

    R <= sum_{k>=1} C(H,k) C(T,k) p^k
      <= exp(H T p)-1
      <= exp(b T_max^2 p)-1.

Thus the extra log2 factor is at most b*T_max^2*p/ln(2), which is
far below one microbit for the current parameters.
"""

from __future__ import annotations

import argparse
import math

from verify_fullsplit_finite_ledger import HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--t-max", type=int, default=17948)
    parser.add_argument("--prefix-pterm-log2", type=float, default=-63.8926492)
    parser.add_argument("--global-pterm-log2", type=float, default=-62.4078758)
    parser.add_argument("--ledger-adjustment-bits", type=float, default=HIGH_INTERVAL_TURNOFF_ADJUSTMENT_BITS)
    args = parser.parse_args()

    if args.block_bits <= 0 or args.t_max <= 0:
        raise SystemExit("block size and T_max must be positive")

    delta = args.global_pterm_log2 - args.prefix_pterm_log2
    y = args.block_bits * args.t_max * args.t_max * (2.0 ** args.global_pterm_log2)
    tail_factor_bits = y / math.log(2.0)
    required = delta + tail_factor_bits
    slack = args.ledger_adjustment_bits - required
    if slack < 0.0:
        raise SystemExit(
            f"high interval adjustment too small: required {required:.12f}, "
            f"ledger has {args.ledger_adjustment_bits:.12f}"
        )

    print("Full-split high-interval p_term adjustment verifier")
    print(f"prefix_pterm_log2,{args.prefix_pterm_log2:.12f}")
    print(f"global_pterm_log2,{args.global_pterm_log2:.12f}")
    print(f"delta_bits,{delta:.12f}")
    print(f"tail_multiplier_extra_bits_bound,{tail_factor_bits:.12g}")
    print(f"required_adjustment_bits,{required:.12f}")
    print(f"ledger_adjustment_bits,{args.ledger_adjustment_bits:.12f}")
    print(f"adjustment_slack_bits,{slack:.12g}")
    print("status,PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
