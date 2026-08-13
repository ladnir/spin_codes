#!/usr/bin/env python3
"""Combine the current disjoint three-band outer failure families.

The entries are outputs of the exact-combinatorial/floating-arithmetic probes
named below.  This small script makes the current *partial* margin and the
contribution of every included family explicit.  It is not a valid final
upper bound yet: decorated-star extremality is unproved for sizes 9..64,
connected sizes 65..308 are absent, disconnected composition is absent, and
the rows other than the exact size-7/8 frontier still need outward versions.
"""

from __future__ import annotations

import math


ROWS = (
    (
        "all_collision_free_PA1",
        -44.922832756341,
        "multipole PA1 collision-free outer-weight diagnostic",
    ),
    (
        "one_band_stars_s_ge_2",
        -59.327072025850,
        "probe_fixed_three_band_star.py",
    ),
    (
        "mixed_two_blocks",
        -56.292258901063,
        "probe_bch_three_band_cell_cap_motifs.py --blocks 2",
    ),
    (
        "mixed_three_blocks",
        -63.934907029017,
        "probe_bch_three_band_cell_cap_motifs.py --blocks 3",
    ),
    (
        "mixed_four_blocks",
        -67.631081527906,
        "probe_bch_three_band_cell_cap_motifs.py --blocks 4",
    ),
    (
        "mixed_five_blocks",
        -65.783158930521,
        "probe_bch_three_band_cell_cap_motifs.py --blocks 5",
    ),
    (
        "mixed_six_blocks",
        -62.620455832711,
        "probe_bch_three_band_cell_cap_motifs.py --blocks 6",
    ),
    (
        "mixed_connected_tail_s_7_to_64",
        -55.299068676898,
        "probe_three_band_decorated_star_tail.py",
    ),
)


def logaddexp2(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(sum(2 ** (value - maximum) for value in values))


def main() -> None:
    total = logaddexp2([value for _name, value, _source in ROWS])
    leader = max(value for _name, value, _source in ROWS)
    print("three-band partial outer diagnostic ledger")
    for name, value, source in ROWS:
        print(
            f"family={name} log2={value:.12f} "
            f"relative_to_leader={value - leader:.12f} source={source}"
        )
    print(f"combined_log2={total:.12f}")
    print(f"margin_below_2^-40={-40.0 - total:.12f}")
    print("disconnected_component_exponential_correction_not_yet_added")
    print("connected_sizes_65_to_308_not_yet_added")
    print("status=INCOMPLETE_DIAGNOSTIC_NOT_A_FAILURE_PROBABILITY_BOUND")


if __name__ == "__main__":
    main()
