#!/usr/bin/env python3
"""Audit locator-polynomial constraints at the BCH minimum-weight boundary.

At the boundary h=Delta for an odd designed distance, the conditions
p_1=...=p_{h-1}=0 imply by Newton identities that the locator polynomial has
only even elementary coefficients below the top degree:

    f(X) = X^h + e_2 X^{h-2} + e_4 X^{h-4} + ... + e_{h-1} X + e_h.

When p_h is free and nonzero, e_h=p_h.  For the [511,259,>=61] parent,
gcd(61,511)=1, so scaling normalizes e_h=1.  This script prints the resulting
locator family and the first post-boundary power-sum constraints that would
cut it down further.
"""

from __future__ import annotations

import argparse
import math

from bch_candidate_params import cyclotomic_coset


def defining_set(n: int, delta: int) -> set[int]:
    out: set[int] = set()
    for e in range(1, delta):
        out.update(cyclotomic_coset(e, n))
    return out


def free_power_sums(n: int, delta: int, limit: int) -> list[int]:
    zeros = defining_set(n, delta)
    return [i for i in range(1, limit + 1) if (i % n) not in zeros]


def locator_terms(h: int, normalized: bool) -> list[tuple[str, int]]:
    terms: list[tuple[str, int]] = [("1", h)]
    for j in range(2, h, 2):
        terms.append((f"e_{j}", h - j))
    terms.append(("1" if normalized else f"e_{h}", 0))
    return terms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=9)
    parser.add_argument("--delta", type=int, default=61)
    parser.add_argument("--h", type=int, default=None)
    parser.add_argument("--post-limit", type=int, default=140)
    args = parser.parse_args()

    n = (1 << args.m) - 1
    h = args.h if args.h is not None else args.delta
    gcd_hn = math.gcd(h, n)
    normalized = gcd_hn == 1
    free_evens = list(range(2, h, 2))
    terms = locator_terms(h, normalized)
    free_powers = free_power_sums(n, args.delta, args.post_limit)
    post_boundary = [p for p in free_powers if p > h]

    print("BCH boundary locator audit")
    print(f"m={args.m}, n={n}, delta={args.delta}, h={h}")
    print(f"gcd(h,n)={gcd_hn}")
    print(f"normalized_constant={int(normalized)}")
    print(f"free_even_coefficients={len(free_evens)}")
    print(f"free_even_e_j={','.join(str(j) for j in free_evens)}")
    print("locator_terms= " + " + ".join(f"{c}*X^{d}" for c, d in terms))
    print("derivative_terms= " + " + ".join(f"{c}*X^{d-1}" for c, d in terms if d % 2 == 1))
    print(f"first_free_power_sums={','.join(str(p) for p in free_powers[:20])}")
    print(f"first_post_boundary_free_power_sums={','.join(str(p) for p in post_boundary[:20])}")
    print(f"num_post_boundary_free_power_sums_to_{args.post_limit}={len(post_boundary)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
