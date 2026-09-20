#!/usr/bin/env python3
"""Positive-shell diagnostics for exact-row wrapped Expand--Convolute.

The exact-row XOR walk determines the intermediate Hamming shell.  Conditional
on shell weight ``u``, the wrapped convolution is bounded with the exact
bivariate transfer matrix.  This is a floating-point parameter diagnostic;
the selected markers must be rechecked with interval arithmetic before a
paper theorem calls the result certified.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from exact_row_bounds import exact_row_shell_step, log_binomial
from expander_bounds import log_slice_tail_bound


LOG_TWO = math.log(2.0)


@dataclass(frozen=True)
class ECMessageWeightTerm:
    message_weight: int
    fixed_message_log2_tail: float
    first_moment_log2_term: float
    dominant_shell: int


def diagnostic_terms(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    memory: int,
    max_message_weight: int,
    wrapping: bool = True,
) -> list[ECMessageWeightTerm]:
    """Bound exact-row EC first-moment terms through one message weight."""
    if not (1 <= max_message_weight <= k):
        raise ValueError("max_message_weight must lie in [1,k]")
    shell_distribution = np.array([0.0], dtype=np.float64)
    slice_cache: dict[int, float] = {0: 0.0}
    terms: list[ECMessageWeightTerm] = []

    for r in range(1, max_message_weight + 1):
        shell_distribution = exact_row_shell_step(shell_distribution, n, row_weight)
        shells = np.flatnonzero(np.isfinite(shell_distribution))
        for shell in shells:
            u = int(shell)
            if u not in slice_cache:
                slice_cache[u] = log_slice_tail_bound(
                    input_weight=u,
                    length=n,
                    cutoff=cutoff,
                    memory=memory,
                    wrapping=wrapping,
                )[0]
        contributions = np.asarray(
            [shell_distribution[u] + slice_cache[int(u)] for u in shells]
        )
        fixed_tail = float(logsumexp(contributions))
        dominant_shell = int(shells[int(np.argmax(contributions))])
        first_moment = log_binomial(k, r) + fixed_tail
        terms.append(
            ECMessageWeightTerm(
                message_weight=r,
                fixed_message_log2_tail=fixed_tail / LOG_TWO,
                first_moment_log2_term=first_moment / LOG_TWO,
                dominant_shell=dominant_shell,
            )
        )
    return terms


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, required=True)
    parser.add_argument("--rate-denominator", type=int, default=5)
    parser.add_argument("--relative-cutoff", type=float, default=0.05)
    parser.add_argument("--row-weight", type=int, required=True)
    parser.add_argument("--memory", type=int, default=21)
    parser.add_argument("--max-message-weight", type=int, default=8)
    parser.add_argument("--nonwrapping", action="store_true")
    args = parser.parse_args()

    k = 1 << args.k_log2
    n = args.rate_denominator * k
    cutoff = math.floor(args.relative_cutoff * n)
    terms = diagnostic_terms(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=args.row_weight,
        memory=args.memory,
        max_message_weight=args.max_message_weight,
        wrapping=not args.nonwrapping,
    )
    print("r\tlog2 fixed tail\tlog2 first-moment term\tdominant shell")
    for term in terms:
        print(
            f"{term.message_weight}\t{term.fixed_message_log2_tail:.6f}\t"
            f"{term.first_moment_log2_term:.6f}\t{term.dominant_shell}"
        )


if __name__ == "__main__":
    main()
