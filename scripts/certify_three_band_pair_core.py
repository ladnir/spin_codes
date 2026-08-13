#!/usr/bin/env python3
"""Exact/outward certificate for the fixed-band pair-collision core.

If every collision tile has size two, the collision pattern in each band is a
matching on the active blocks.  Dropping connectivity and the prohibition on
using the same block pair in two bands leaves three independent matchings.

At Chernoff pole q=1/20, let F_S be the one-codeword cell-cap factor when the
block is paired in precisely the bands in S.  This verifier proves the simple
dyadic envelope

    F_S <= 2^(|S|-57).

Consequently a matching edge has activity 2^2=4 after extracting the common
per-vertex factor 2^-57.  The weighted number J_s of matchings satisfies

    J_0 = J_1 = 1,
    J_s = J_(s-1) + 4 (s-1) J_(s-2).

For a connected s-block pattern, the affine tile map has at most
2^14 * 64^(s-1) placements.  The three independent matching relaxation,
exact-length puncture adjustment, and graph incremental-MGF factor therefore
give the exact rational bound used below.  Finally,
(0.181/(1/20))^106 = (181/50)^106 is charged once for the support-106
Chernoff cutoff.
"""

from __future__ import annotations

import math
from fractions import Fraction

from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from three_band_exact_moments import (
    factor_table_upper,
    greedy_factor_bound_upper,
)
from certify_three_band_exact_length import (
    BAND_ZERO_COLUMNS,
    graph_factor_upper,
    puncture_selection_adjustment,
)


BLOCKS = 1 << 14
LANES = 1 << 6
SUPPORT_CUTOFF = 106
FIRST_CORE_SIZE = 7
LAST_CORE_SIZE = 323
ROOT_BITS = 160


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    pole = Fraction(1, 20)
    tables = {
        (band, cluster_size): factor_table_upper(
            BAND_SIZES[band], cluster_size, pole, root_bits=ROOT_BITS
        )
        for band in range(3)
        for cluster_size in (1, 2)
    }

    print("exact three-band pair-collision core certificate")
    print(f"pole={pole} root_bits={ROOT_BITS}")
    for mask in range(8):
        degree = mask.bit_count()
        factors = tuple(
            tables[band, 2 if mask & (1 << band) else 1]
            for band in range(3)
        )
        value = greedy_factor_bound_upper(
            states, caps, diagonal_mass, by_total, factors
        )
        puncture_multiplier = tuple(
            1 + (1 / pole - 1) * Fraction(weight, BAND_ZERO_COLUMNS)
            for weight in range(BAND_ZERO_COLUMNS + 1)
        )
        punctured_factors = (
            tuple(
                factors[0][weight] * puncture_multiplier[weight]
                for weight in range(BAND_ZERO_COLUMNS + 1)
            ),
            factors[1],
            factors[2],
        )
        punctured_value = greedy_factor_bound_upper(
            states, caps, diagonal_mass, by_total, punctured_factors
        )
        value *= puncture_selection_adjustment(value, punctured_value)
        envelope = Fraction(1, 1 << (57 - degree))
        if value > envelope:
            raise SystemExit(
                f"pair core: mask {mask:03b} exceeds 2^({degree}-57)"
            )
        print(
            f"mask={mask:03b} degree={degree} "
            f"factor_log2_upper={log2_fraction(value):.12f} "
            f"envelope_log2={degree - 57}"
        )

    weighted_matchings = [1, 1]
    for size in range(2, LAST_CORE_SIZE + 1):
        weighted_matchings.append(
            weighted_matchings[-1]
            + 4 * (size - 1) * weighted_matchings[-2]
        )

    terms: list[Fraction] = []
    for size in range(FIRST_CORE_SIZE, LAST_CORE_SIZE + 1):
        numerator = (
            BLOCKS
            * LANES ** (size - 1)
            * weighted_matchings[size] ** 3
        )
        denominator = math.factorial(size) * (1 << (57 * size))
        terms.append(
            Fraction(numerator, denominator) * graph_factor_upper(pole)
        )

    tilted_total = sum(terms, Fraction(0))
    cutoff_conversion = Fraction(181, 50) ** SUPPORT_CUTOFF
    target_total = tilted_total * cutoff_conversion
    if target_total > Fraction(1, 1 << 149):
        raise SystemExit("pair core: corrected total exceeds 2^-149")

    print(
        f"size_{FIRST_CORE_SIZE}_term_log2="
        f"{log2_fraction(terms[0]):.12f}"
    )
    print(f"tilted_total_log2={log2_fraction(tilted_total):.12f}")
    print(
        f"cutoff_conversion_log2={log2_fraction(cutoff_conversion):.12f}"
    )
    print(f"target_total_log2={log2_fraction(target_total):.12f}")
    print("target_total_le_2^-149=PASS")
    print("status=EXACT_COUNTS_AND_RATIONAL_OUTWARD_FACTORS")


if __name__ == "__main__":
    main()
