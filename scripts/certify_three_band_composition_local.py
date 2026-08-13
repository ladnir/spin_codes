#!/usr/bin/env python3
"""Exact/outward local inequalities for the three-band composition lemma.

This verifier covers two reusable rows at Chernoff pole 1/20.

* A high tile of degree r whose blocks are pessimistically paired in both
  external bands has one-codeword factor F_b(r,2,2).  After choosing a parent,
  anchor, root band, and remaining lanes, every such child batch has activity
  at most 2^-112, uniformly for r=3..64 and all root bands.

* When a degree-three batch attaches through a later degree-three root tile,
  the final shared tile has degree four.  The exact static Holder types are
  (3,4,2) for the earlier anchor, two copies of (3,2,2), and three copies of
  (2,4,2), up to band permutation.  Including root/lane choices, the one
  support-106 conversion, and the one graph-codimension charge gives at most
  ``(3/4) 2^-39`` after the exact-length puncture adjustment.

All BCH factors use exact rational diagonal moments, 160-bit upward dyadic
roots, and the exact cell-cap greedy relaxation.  Printed logarithms are only
diagnostic; every PASS comparison is exact rational arithmetic.
"""

from __future__ import annotations

import math
from fractions import Fraction

from outward_log2 import log2_fraction
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
    worst_puncture_block_adjustment,
)


MAXIMUM_ACTIVE_BLOCKS = 323
ROOT_BITS = 160
POLE = Fraction(1, 20)
CUTOFF_CONVERSION = Fraction(181, 50) ** 106


def main() -> None:
    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    tables = {
        (band, degree): factor_table_upper(
            BAND_SIZES[band], degree, POLE, root_bits=ROOT_BITS
        )
        for band in range(3)
        for degree in range(1, 65)
    }
    factors: dict[tuple[int, int, int], Fraction] = {}
    punctured_factors: dict[tuple[int, int, int], Fraction] = {}
    adjusted_factors: dict[tuple[int, int, int], Fraction] = {}

    def factor(kind: tuple[int, int, int]) -> Fraction:
        if kind not in factors:
            factors[kind] = greedy_factor_bound_upper(
                states,
                caps,
                diagonal_mass,
                by_total,
                tuple(tables[band, kind[band]] for band in range(3)),
            )
        return factors[kind]

    def punctured_factor(kind: tuple[int, int, int]) -> Fraction:
        if kind not in punctured_factors:
            puncture_multiplier = tuple(
                1
                + (1 / POLE - 1) * Fraction(weight, BAND_ZERO_COLUMNS)
                for weight in range(BAND_ZERO_COLUMNS + 1)
            )
            punctured_tables = (
                tuple(
                    tables[0, kind[0]][weight] * puncture_multiplier[weight]
                    for weight in range(BAND_ZERO_COLUMNS + 1)
                ),
                tables[1, kind[1]],
                tables[2, kind[2]],
            )
            punctured_factors[kind] = greedy_factor_bound_upper(
                states,
                caps,
                diagonal_mass,
                by_total,
                punctured_tables,
            )
        return punctured_factors[kind]

    def adjusted_factor(kind: tuple[int, int, int]) -> Fraction:
        if kind not in adjusted_factors:
            full = factor(kind)
            punctured = punctured_factor(kind)
            adjusted_factors[kind] = full * puncture_selection_adjustment(
                full, punctured
            )
        return adjusted_factors[kind]

    worst_child = Fraction(0)
    worst_child_kind = (0, 0)
    worst_root = Fraction(0)
    worst_root_kind = (0, 0)
    worst_root_degree_at_least_five = Fraction(0)
    worst_root_degree_at_least_five_kind = (0, 0)
    for root_band in range(3):
        for degree in range(3, 65):
            kind = [2, 2, 2]
            kind[root_band] = degree
            value = (
                adjusted_factor(tuple(kind)) ** degree
                * MAXIMUM_ACTIVE_BLOCKS
                * 378
                * math.comb(63, degree - 1)
            )
            if value > worst_child:
                worst_child = value
                worst_child_kind = (root_band, degree)
            root_value = adjusted_factor(tuple(kind)) ** degree * math.comb(64, degree)
            if root_value > worst_root:
                worst_root = root_value
                worst_root_kind = (root_band, degree)
            if degree >= 5 and root_value > worst_root_degree_at_least_five:
                worst_root_degree_at_least_five = root_value
                worst_root_degree_at_least_five_kind = (root_band, degree)
    if worst_child > Fraction(1, 1 << 112):
        raise SystemExit("composition local: decorated child exceeds 2^-112")
    if worst_root > Fraction(1, 1 << 125):
        raise SystemExit("composition local: decorated root exceeds 2^-125")
    if worst_root_degree_at_least_five > Fraction(1, 1 << 165):
        raise SystemExit("composition local: degree>=5 root exceeds 2^-165")

    worst_adjacent = Fraction(0)
    worst_adjacent_bands = (0, 0)
    for earlier_band in range(3):
        for later_band in range(3):
            if earlier_band == later_band:
                continue
            earlier_other = [2, 2, 2]
            earlier_other[earlier_band] = 3
            earlier_anchor = earlier_other.copy()
            earlier_anchor[later_band] = 4
            later = [2, 2, 2]
            later[later_band] = 4
            raw = (
                adjusted_factor(tuple(earlier_other)) ** 2
                * adjusted_factor(tuple(earlier_anchor))
                * adjusted_factor(tuple(later)) ** 3
            )
            # Choose the later root tile and its three lanes; select its
            # attachment lane/band; then choose the other two earlier-root
            # lanes around the determined earlier anchor.
            placement = (
                768
                * math.comb(64, 3)
                * (3 * 126)
                * math.comb(63, 2)
            )
            value = raw * placement * CUTOFF_CONVERSION * graph_factor_upper(POLE)
            if value > worst_adjacent:
                worst_adjacent = value
                worst_adjacent_bands = (earlier_band, later_band)
    adjacent_envelope = Fraction(3, 4 * (1 << 39))
    if worst_adjacent > adjacent_envelope:
        raise SystemExit("composition local: adjacent 3/3 base exceeds (3/4) 2^-39")

    worst_one_residual = Fraction(0)
    worst_one_residual_kind = (0, 0)
    for root_band in range(3):
        external = [band for band in range(3) if band != root_band]
        for degree in range(5, 65):
            ordinary = [1, 1, 1]
            ordinary[root_band] = degree
            anchors = []
            for external_band in external:
                anchor = ordinary.copy()
                anchor[external_band] = 2
                anchors.append(tuple(anchor))
            residual = [1, 1, 1]
            residual[external[0]] = 2
            residual[external[1]] = 2
            # This family is far from the final bottleneck.  Use the
            # pointwise one-coordinate puncture adjustment instead of solving
            # a second cell-cap problem for hundreds of degree/kind pairs.
            raw = (
                factor(tuple(ordinary)) ** (degree - 2)
                * factor(anchors[0])
                * factor(anchors[1])
                * factor(tuple(residual))
                * worst_puncture_block_adjustment(POLE) ** (degree + 1)
            )
            placement = (
                768
                * math.comb(64, degree)
                * math.comb(degree, 2)
                * 4
                * 63
            )
            value = raw * placement * CUTOFF_CONVERSION * graph_factor_upper(POLE)
            if value > worst_one_residual:
                worst_one_residual = value
                worst_one_residual_kind = (root_band, degree)
    if worst_one_residual > Fraction(1, 1 << 50):
        raise SystemExit("composition local: degree>=5 one-residual base exceeds 2^-50")

    print("exact three-band composition local certificate")
    print(f"pole={POLE} root_bits={ROOT_BITS}")
    print(
        f"worst_decorated_child_band_degree={worst_child_kind} "
        f"log2_upper={log2_fraction(worst_child).hi}"
    )
    print("worst_decorated_child_le_2^-112=PASS")
    print(
        f"worst_decorated_root_band_degree={worst_root_kind} "
        f"log2_upper={log2_fraction(worst_root).hi}"
    )
    print("worst_decorated_root_le_2^-125=PASS")
    print(
        "worst_decorated_root_degree_at_least_five_band_degree="
        f"{worst_root_degree_at_least_five_kind} "
        f"log2_upper={log2_fraction(worst_root_degree_at_least_five).hi}"
    )
    print("worst_decorated_root_degree_at_least_five_le_2^-165=PASS")
    print(
        f"worst_adjacent_3_3_bands={worst_adjacent_bands} "
        f"log2_upper={log2_fraction(worst_adjacent).hi}"
    )
    print("worst_adjacent_3_3_le_(3/4)*2^-39=PASS")
    print(
        f"worst_one_residual_band_degree={worst_one_residual_kind} "
        f"log2_upper={log2_fraction(worst_one_residual).hi}"
    )
    print("worst_one_residual_degree_at_least_five_le_2^-50=PASS")
    print("status=EXACT_RATIONAL_OUTWARD_LOCAL_COMPOSITION")


if __name__ == "__main__":
    main()
