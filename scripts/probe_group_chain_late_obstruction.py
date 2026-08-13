#!/usr/bin/env python3
"""Exact late-suffix obstruction for contiguous group-chain Riffle layouts.

The calculation is a lower bound on construction failure, not an upper-bound
certificate.  It uses only the maximum possible 64 output bits per group and
the guaranteed 40-dimensional kernel of a codimension-24 graph restricted to
one 64-coordinate message block.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
B = 32768
D = 188743
T = D // 64
UNPUNCTURED_DATA_BLOCKS = 16384 - 256
GRAPH_CODIMENSION = 24


def no_complete_pair_count(pair_count: int, singles: int, chosen: int) -> int:
    """Subsets choosing at most one member from each distinguished pair."""

    total = 0
    for from_pairs in range(max(0, chosen - singles), min(pair_count, chosen) + 1):
        total += (
            math.comb(pair_count, from_pairs)
            * (1 << from_pairs)
            * math.comb(singles, chosen - from_pairs)
        )
    return total


def load_spectrum() -> dict[int, int]:
    result: dict[int, int] = {}
    for raw in (ROOT / "EBCH128_64.wd").read_text(encoding="utf-8").splitlines():
        fields = raw.split()
        if len(fields) == 2 and fields[0].isdigit():
            result[int(fields[0])] = int(fields[1])
    if sum(result.values()) != 1 << 64:
        raise SystemExit("late obstruction: invalid EBCH spectrum mass")
    return result


def main() -> None:
    suffix_capacity = 64 * T
    if suffix_capacity > D or 64 * (T + 1) <= D:
        raise SystemExit("late obstruction: suffix threshold mismatch")

    parity_in_suffix = Fraction(math.comb(T, 2), math.comb(B, 2))
    remaining_positions = B - 2
    remaining_chosen = T - 2
    singles = remaining_positions - 2 * UNPUNCTURED_DATA_BLOCKS
    no_pair = Fraction(
        no_complete_pair_count(
            UNPUNCTURED_DATA_BLOCKS, singles, remaining_chosen
        ),
        math.comb(remaining_positions, remaining_chosen),
    )
    failure_lower = parity_in_suffix * (1 - no_pair)

    # Diagnostic for the striped repair: expected number of the same
    # one-data-block, graph-zero words whose entire w+2 group support is late.
    # This is also a sufficient-event lower bound, not a proof upper bound.
    spectrum = load_spectrum()
    striped_late_expectation = sum(
        Fraction(
            UNPUNCTURED_DATA_BLOCKS
            * count
            * math.comb(T, weight + 2),
            (1 << GRAPH_CODIMENSION) * math.comb(B, weight + 2),
        )
        for weight, count in spectrum.items()
        if weight and weight + 2 <= T
    )

    print("group-chain universal late-suffix obstruction")
    print(f"B={B} d={D} T=floor(d/64)={T} suffix_capacity={suffix_capacity}")
    print(f"unpunctured_disjoint_data_pairs={UNPUNCTURED_DATA_BLOCKS}")
    print(f"parity_pair_late_log2={math.log2(float(parity_in_suffix)):.12f}")
    print(f"conditional_no_data_pair_late_log2={math.log2(float(no_pair)):.12f}")
    print(f"contiguous_failure_lower_log2={math.log2(float(failure_lower)):.12f}")
    print(
        "striped_same_family_all_late_expectation_log2="
        f"{math.log2(float(striped_late_expectation)):.12f}"
    )


if __name__ == "__main__":
    main()
