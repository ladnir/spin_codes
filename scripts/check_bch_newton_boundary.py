#!/usr/bin/env python3
"""Audit the Newton-identity boundary data for primitive BCH codes.

For a primitive narrow-sense BCH code of length n=2^m-1 and designed distance
Delta, a weight-h support T satisfies p_i(T)=sum_{x in T} x^i=0 for every
exponent i in the defining cyclotomic closure of {1,...,Delta-1}.

At the minimum-weight edge h=Delta, Newton's identities make the first free
power sum p_h especially important.  If gcd(h,n)=1 and p_h is nonzero, scaling
the locators by beta can normalize p_h to one.  This is the same normalization
used in the locator-polynomial BCH minimum-distance literature.
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


def coset_reps(n: int, delta: int) -> list[int]:
    seen: set[int] = set()
    reps: list[int] = []
    for e in range(1, delta):
        c = cyclotomic_coset(e, n)
        if not any(x in seen for x in c):
            reps.append(e)
            seen.update(c)
    return reps


def odd_elementary_free_indices(h: int) -> list[int]:
    # If p_1,...,p_{h-1}=0, Newton gives e_j=0 for odd j<h.
    return [j for j in range(1, h) if j % 2 == 0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=9)
    parser.add_argument("--delta", type=int, default=61)
    parser.add_argument("--h", type=int, default=None)
    parser.add_argument("--show-free-powers-to", type=int, default=140)
    args = parser.parse_args()

    n = (1 << args.m) - 1
    h = args.h if args.h is not None else args.delta
    zeros = defining_set(n, args.delta)
    reps = coset_reps(n, args.delta)
    free_powers = [i for i in range(1, args.show_free_powers_to + 1) if (i % n) not in zeros]
    first_free = free_powers[0] if free_powers else None
    print("BCH Newton boundary audit")
    print(f"m={args.m}, n={n}, delta={args.delta}, h={h}")
    print(f"redundancy={len(zeros)}, dimension={n - len(zeros)}")
    print(f"coset_reps={','.join(str(x) for x in reps)}")
    print(f"num_cosets={len(reps)}")
    print(f"first_free_power={first_free}")
    print(f"p_h_forced_zero={(h % n) in zeros}")
    print(f"gcd(h,n)={math.gcd(h,n)}")
    scaling_kernel = math.gcd(h, n)
    print(f"p_h_nonzero_at_weight_h={h % 2 == 1}")
    print(f"p_h_scaling_image_size={n // scaling_kernel}")
    print(f"p_h_scaling_classes={scaling_kernel}")
    print(f"p_h_scaling_normalizable_to_one={scaling_kernel == 1}")
    evens = odd_elementary_free_indices(h)
    print(f"if p_1..p_{h-1}=0, odd e_j below h vanish; free even e_j count={len(evens)}")
    print(f"free_even_e_j_prefix={','.join(str(x) for x in evens[:20])}")
    print(f"free_power_sums_to_{args.show_free_powers_to}={','.join(str(x) for x in free_powers)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
