#!/usr/bin/env python3
"""List primitive narrow-sense binary BCH candidates near rate 1/2.

For length n=2^m-1, the generator degree is the size of the union of binary
cyclotomic cosets meeting {1,...,delta-1}.  The dimension is n-degree.  This
script reports the largest odd designed distance with dimension at least a
target message size, defaulting to b=2^(m-1).
"""

from __future__ import annotations

import argparse


def cyclotomic_coset(e: int, n: int) -> list[int]:
    out: list[int] = []
    x = e % n
    while x not in out:
        out.append(x)
        x = (2 * x) % n
    return out


def bch_dimension(m: int, designed_delta: int) -> tuple[int, int, int]:
    n = (1 << m) - 1
    seen: set[int] = set()
    for e in range(1, designed_delta):
        seen.update(cyclotomic_coset(e, n))
    degree = len(seen)
    return n, n - degree, degree


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-min", type=int, default=7)
    parser.add_argument("--m-max", type=int, default=11)
    parser.add_argument("--delta-max", type=int, default=255)
    parser.add_argument(
        "--target",
        type=int,
        default=None,
        help="Target dimension. Default for each m is 2^(m-1).",
    )
    args = parser.parse_args()

    print("m,n,target_dim,best_designed_delta,dimension,generator_degree")
    for m in range(args.m_min, args.m_max + 1):
        n = (1 << m) - 1
        target = args.target if args.target is not None else (1 << (m - 1))
        best: tuple[int, int, int] | None = None
        for delta in range(3, args.delta_max + 1, 2):
            _, dim, degree = bch_dimension(m, delta)
            if dim >= target:
                best = (delta, dim, degree)
        if best is None:
            print(f"{m},{n},{target},,,")
        else:
            delta, dim, degree = best
            print(f"{m},{n},{target},{delta},{dim},{degree}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
