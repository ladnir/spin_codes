#!/usr/bin/env python3
"""Bell-overcounted decorated-star envelope for the size-7..64 frontier.

The degree-type relaxation identifies the extremal mixed connected profile on
``s >= 7`` blocks as one ``(s-1)``-cluster, one external block, and one spoke
from that block in each of the other bands.  This script evaluates that
profile with the solver-free BCH cell-cap factors and then multiplies by
``Bell(s)^3``, counting every triple of set partitions, compatible or not.

The numerical envelope uses tilt 0.075 at s=7 and 0.05 from s=8 through the
maximum single-tile degree 64.  This is a frontier, not the entire connected
tail: a connected pattern can span several tiles.  Minimum BCH weight 22 and
physical group capacity 64 imply only ``s <= floor(64*106/22) = 308`` under
the support-106 cutoff.  Sizes 65..308 therefore require a separate envelope.

The displayed rows are also conditional on the decorated-star extremal lemma
checked at individual sizes by ``probe_three_band_tail_degree_mip.py``.
"""

from __future__ import annotations

import math

from probe_bch_three_band_cell_cap_motifs import (
    BAND_SIZES,
    cell_caps,
    factor_table,
    greedy_factor_bound,
    logaddexp2,
)
from probe_bch_three_band_transport_lp import (
    exact_rows,
    remove_trivial_codewords,
    triples,
)


BLOCKS = 16384
LANES = 64
GRAPH_CODIMENSION = 24
TARGET_POLE = 0.181
SUPPORT_CUTOFF = 106
MAXIMUM_CLUSTER = 64


def bell_numbers(maximum: int) -> list[int]:
    values = [1]
    for size in range(1, maximum + 1):
        values.append(
            sum(math.comb(size - 1, prefix) * values[prefix] for prefix in range(size))
        )
    return values


def main() -> None:
    states = triples()
    exact = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, exact)
    bells = bell_numbers(MAXIMUM_CLUSTER)
    factor_cache = {}
    value_cache = {}

    def factor(tilt: float, band: int, degree: int):
        key = (tilt, band, degree)
        if key not in factor_cache:
            factor_cache[key] = factor_table(BAND_SIZES[band], degree, tilt)
        return factor_cache[key]

    def value(tilt: float, kind: tuple[int, int, int]) -> float:
        key = (tilt, kind)
        if key not in value_cache:
            value_cache[key] = greedy_factor_bound(
                states,
                caps,
                diagonal_mass,
                by_total,
                tuple(factor(tilt, band, kind[band]) for band in range(3)),
            )
        return value_cache[key]

    rows = []
    total = -math.inf
    for size in range(7, MAXIMUM_CLUSTER + 1):
        tilt = 0.075 if size == 7 else 0.05
        dense_values = []
        for dense_band in range(3):
            other_bands = [band for band in range(3) if band != dense_band]
            external = [1, 1, 1]
            external[other_bands[0]] = 2
            external[other_bands[1]] = 2
            first_arm = [1, 1, 1]
            first_arm[dense_band] = size - 1
            first_arm[other_bands[0]] = 2
            second_arm = [1, 1, 1]
            second_arm[dense_band] = size - 1
            second_arm[other_bands[1]] = 2
            core = [1, 1, 1]
            core[dense_band] = size - 1
            dense_values.append(
                math.log2(BLOCKS)
                + (size - 1) * math.log2(LANES)
                - math.lgamma(size + 1) / math.log(2)
                - GRAPH_CODIMENSION
                + math.log2(value(tilt, tuple(external)))
                + math.log2(value(tilt, tuple(first_arm)))
                + math.log2(value(tilt, tuple(second_arm)))
                + (size - 3) * math.log2(value(tilt, tuple(core)))
                + SUPPORT_CUTOFF * math.log2(TARGET_POLE / tilt)
            )
        extremal = max(dense_values)
        bound = extremal + 3 * math.log2(bells[size])
        rows.append((bound, size, tilt, extremal))
        total = logaddexp2(total, bound)

    print("three-band decorated-star size-7..64 frontier")
    print(
        f"sizes=7..{MAXIMUM_CLUSTER} target_pole={TARGET_POLE:.12f} "
        f"support_cutoff={SUPPORT_CUTOFF}"
    )
    print(f"tail_log2={total:.12f}")
    for bound, size, tilt, extremal in sorted(rows, reverse=True)[:12]:
        print(
            f"size={size} tilt={tilt:.12f} extremal_individual_log2="
            f"{extremal:.12f} bell_overcounted_log2={bound:.12f}"
        )
    print(
        "status=CONDITIONAL_FRONTIER_ONLY_MISSING_EXTREMAL_LEMMA_AND_SIZES_65_308"
    )


if __name__ == "__main__":
    main()
