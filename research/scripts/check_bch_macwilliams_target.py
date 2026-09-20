#!/usr/bin/env python3
"""Print the finite MacWilliams correction target for BCH spectrum envelopes.

For a parent [n,k_parent] code C and a desired [n,k_sub] subcode, MacWilliams
gives

    A_h(C) = 2^{-(n-k_parent)} (binom(n,h) + S_h),

where S_h is the nonzero-dual Krawtchouk correction.  If the block certificate
can tolerate B bits of loss relative to the random-like [n,k_sub] spectrum,
then the parent spectrum is sufficient whenever

    1 + S_h / binom(n,h) <= 2^{B - (k_parent-k_sub)}.

This script reports that correction budget.  It is a theorem-target calculator,
not a proof of the BCH spectrum bound.
"""

from __future__ import annotations

import argparse
import math


def parse_weights(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            vals = [int(x) for x in part.split(":")]
            if len(vals) == 2:
                lo, hi = vals
                step = 1
            elif len(vals) == 3:
                lo, hi, step = vals
            else:
                raise ValueError(f"bad range {part!r}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=511)
    parser.add_argument("--parent-dim", type=int, default=259)
    parser.add_argument("--subcode-dim", type=int, default=256)
    parser.add_argument("--slack-bits", type=float, default=56.75)
    parser.add_argument("--weights", default="61,70,89,100,150,220")
    args = parser.parse_args()

    dim_gap = args.parent_dim - args.subcode_dim
    envelope_factor = args.slack_bits - dim_gap
    correction = math.log2(2.0**envelope_factor - 1.0) if envelope_factor > -40.0 else float("-inf")

    print("BCH MacWilliams correction target")
    print(f"n={args.n}, parent_dim={args.parent_dim}, subcode_dim={args.subcode_dim}")
    print(f"dim_gap={dim_gap}, local_slack_budget={args.slack_bits:.6f}")
    print(f"required log2(1 + S_h/binom(n,h)) <= {envelope_factor:.6f}")
    print(f"equivalently log2(S_h/binom(n,h)) <= {correction:.6f}")
    print("h,log2_allowed_correction_over_binom")
    for h in parse_weights(args.weights):
        print(f"{h},{correction:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
