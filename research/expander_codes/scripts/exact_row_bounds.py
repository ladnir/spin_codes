#!/usr/bin/env python3
"""Positive-coefficient diagnostics for exact-row expand--accumulate codes.

Each selected expander row is uniform on the Hamming sphere of radius ``d``.
The XOR of the selected rows is therefore exchangeable, and its Hamming
weight evolves as a one-dimensional Markov chain. Conditional on shell
weight ``u``, the accumulator low-weight event is a hypergeometric tail.

This module deliberately separates fast floating-point diagnostics from a
future rigorous certificate verifier. No alternating Krawtchouk sums occur
in the diagnostic path.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.stats import hypergeom


LOG_TWO = math.log(2.0)


def log_binomial(n: int, k: int) -> float:
    """Return log(binomial(n,k)), with -inf outside the natural range."""
    if k < 0 or k > n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def accumulator_slice_logtail(n: int, cutoff: int, shell_weight: int) -> float:
    """Return log Pr[accumulator weight <= cutoff | input shell weight u].

    For a uniformly random ``u``-subset of input positions, the event is that
    at least ceil(u/2) selected positions lie among the first ``cutoff``
    positions. Thus it is a positive hypergeometric upper tail.
    """
    u = shell_weight
    if not (0 <= cutoff <= n and 0 <= u <= n):
        raise ValueError("cutoff and shell weight must lie in [0,n]")
    if u == 0:
        return 0.0
    threshold = (u + 1) // 2
    if threshold > min(u, cutoff):
        return -math.inf
    value = float(hypergeom.logsf(threshold - 1, n, cutoff, u))
    if math.isfinite(value):
        return value

    # Some SciPy builds form sf by subtracting a cdf. The direct positive
    # log-sum is slower but preserves very small tails needed in diagnostics.
    upper = min(u, cutoff)
    selected = np.arange(threshold, upper + 1, dtype=np.int64)
    terms = (
        gammaln(cutoff + 1)
        - gammaln(selected + 1)
        - gammaln(cutoff - selected + 1)
        + gammaln(n - cutoff + 1)
        - gammaln(u - selected + 1)
        - gammaln(n - cutoff - u + selected + 1)
        - log_binomial(n, u)
    )
    return float(logsumexp(terms))


def accumulator_slice_logtails(n: int, cutoff: int, max_weight: int) -> np.ndarray:
    """Precompute accumulator slice log tails for shells 0..max_weight."""
    if max_weight > n:
        raise ValueError("max_weight cannot exceed n")
    result = np.empty(max_weight + 1, dtype=np.float64)
    for u in range(max_weight + 1):
        result[u] = accumulator_slice_logtail(n, cutoff, u)
    return result


def exact_row_shell_step(
    log_distribution: np.ndarray, n: int, row_weight: int
) -> np.ndarray:
    """Advance the exact-row XOR shell chain by one row in log space."""
    d = row_weight
    if not (0 <= d <= n):
        raise ValueError("row_weight must lie in [0,n]")
    old_max = len(log_distribution) - 1
    new_max = min(n, old_max + d)
    result = np.full(new_max + 1, -math.inf, dtype=np.float64)

    for u in np.flatnonzero(np.isfinite(log_distribution)):
        j_lo = max(0, d - (n - int(u)))
        j_hi = min(d, int(u))
        overlaps = np.arange(j_lo, j_hi + 1, dtype=np.int64)
        next_shells = int(u) + d - 2 * overlaps
        log_transitions = hypergeom.logpmf(overlaps, n, int(u), d)
        np.logaddexp.at(
            result,
            next_shells,
            log_distribution[u] + log_transitions,
        )
    return result


def exact_row_shell_logdistribution(n: int, row_weight: int, rows: int) -> np.ndarray:
    """Return log probabilities for the XOR weight after ``rows`` rows."""
    if rows < 0:
        raise ValueError("rows must be nonnegative")
    distribution = np.array([0.0], dtype=np.float64)
    for _ in range(rows):
        distribution = exact_row_shell_step(distribution, n, row_weight)
    return distribution


def fixed_message_logtail(
    n: int,
    cutoff: int,
    log_shell_distribution: np.ndarray,
    slice_logtails: np.ndarray | None = None,
) -> float:
    """Return the log low-output probability for one fixed message."""
    max_shell = len(log_shell_distribution) - 1
    if slice_logtails is None:
        slice_logtails = accumulator_slice_logtails(n, cutoff, max_shell)
    if len(slice_logtails) <= max_shell:
        raise ValueError("slice_logtails does not cover every shell")
    return float(logsumexp(log_shell_distribution + slice_logtails[: max_shell + 1]))


@dataclass(frozen=True)
class MessageWeightTerm:
    message_weight: int
    fixed_message_log2_tail: float
    first_moment_log2_term: float
    most_likely_contributing_shell: int


def diagnostic_terms(
    *, k: int, n: int, cutoff: int, row_weight: int, max_message_weight: int
) -> list[MessageWeightTerm]:
    """Compute exact-row EA first-moment terms up to a message-weight cutoff."""
    if not (1 <= max_message_weight <= k):
        raise ValueError("max_message_weight must lie in [1,k]")
    max_shell = min(n, row_weight * max_message_weight)
    slice_tails = accumulator_slice_logtails(n, cutoff, max_shell)
    shell_distribution = np.array([0.0], dtype=np.float64)
    terms: list[MessageWeightTerm] = []

    for r in range(1, max_message_weight + 1):
        shell_distribution = exact_row_shell_step(shell_distribution, n, row_weight)
        contributions = shell_distribution + slice_tails[: len(shell_distribution)]
        fixed_logtail = float(logsumexp(contributions))
        shell = int(np.argmax(contributions))
        term = log_binomial(k, r) + fixed_logtail
        terms.append(
            MessageWeightTerm(
                message_weight=r,
                fixed_message_log2_tail=fixed_logtail / LOG_TWO,
                first_moment_log2_term=term / LOG_TWO,
                most_likely_contributing_shell=shell,
            )
        )
    return terms


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, required=True)
    parser.add_argument("--rate-denominator", type=int, default=5)
    parser.add_argument("--relative-cutoff", type=float, default=0.05)
    parser.add_argument("--row-weight", type=int, required=True)
    parser.add_argument("--max-message-weight", type=int, default=64)
    args = parser.parse_args()

    k = 1 << args.k_log2
    n = args.rate_denominator * k
    cutoff = math.floor(args.relative_cutoff * n)
    terms = diagnostic_terms(
        k=k,
        n=n,
        cutoff=cutoff,
        row_weight=args.row_weight,
        max_message_weight=args.max_message_weight,
    )
    print("r\tlog2 fixed tail\tlog2 first-moment term\tdominant shell")
    for term in terms:
        print(
            f"{term.message_weight}\t{term.fixed_message_log2_tail:.6f}\t"
            f"{term.first_moment_log2_term:.6f}\t"
            f"{term.most_likely_contributing_shell}"
        )


if __name__ == "__main__":
    main()
