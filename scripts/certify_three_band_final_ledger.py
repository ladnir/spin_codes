#!/usr/bin/env python3
"""Exact/outward support-at-most-106 Riffle failure ledger.

Components of sizes one through six are charged directly at the target pole
181/1000.  The size-two-through-six activities deliberately use the broader
all-pattern envelopes, so treating them again as connected activities only
overcounts.  Every remaining connected component is in one of three common
tail-pole families: a size-at-least-six star, the degree-two pair core, or the
high-cluster composition.

For a word of total group support at most 106 containing a tail-pole component,
the conversion
``(181/50)^106`` and the graph factor are charged exactly once globally.
Collections of components are bounded by the positive geometric series
``x/(1-x)``, which is coarser than the exponential formula and hence safe.
All arithmetic after the already outward-certified local factors is exact
rational arithmetic.
"""

from __future__ import annotations

from fractions import Fraction

from certify_three_band_exact_length import (
    BAND_ZERO_COLUMNS,
    graph_factor_upper,
    puncture_selection_adjustment,
)
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


BLOCKS = 1 << 14
TARGET_POLE = Fraction(181, 1000)
TAIL_POLE = Fraction(1, 20)
ROOT_BITS = 160
SUPPORT_CUTOFF = 106

# Post-graph dyadic envelopes proved by
# certify_three_band_small_motifs.py.  They include all compatible patterns,
# not only connected ones, which is a useful safe overcount here.
SMALL_POST_GRAPH = sum(
    (Fraction(1, 1 << bits) for bits in (56, 62, 63, 60, 53)),
    Fraction(0),
)

# Post-graph, post-cutoff connected tail envelopes proved respectively by
# certify_three_band_composition_ledger.py,
# certify_three_band_star.py, and certify_three_band_pair_core.py.
TAIL_POST_GRAPH = (
    Fraction(1, 1 << 41)
    + Fraction(1, 1 << 53)
    + Fraction(1, 1 << 149)
)


def singleton_activity() -> Fraction:
    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    tables = tuple(
        factor_table_upper(
            BAND_SIZES[band], 1, TARGET_POLE, root_bits=ROOT_BITS
        )
        for band in range(3)
    )
    full = greedy_factor_bound_upper(
        states, caps, diagonal_mass, by_total, tables
    )
    puncture_multiplier = tuple(
        1
        + (1 / TARGET_POLE - 1) * Fraction(weight, BAND_ZERO_COLUMNS)
        for weight in range(BAND_ZERO_COLUMNS + 1)
    )
    punctured_tables = (
        tuple(
            tables[0][weight] * puncture_multiplier[weight]
            for weight in range(BAND_ZERO_COLUMNS + 1)
        ),
        tables[1],
        tables[2],
    )
    punctured = greedy_factor_bound_upper(
        states, caps, diagonal_mass, by_total, punctured_tables
    )
    adjusted = full * puncture_selection_adjustment(full, punctured)
    return BLOCKS * adjusted


def main() -> None:
    target_graph = graph_factor_upper(TARGET_POLE)
    tail_graph = graph_factor_upper(TAIL_POLE)
    cutoff_conversion = (TARGET_POLE / TAIL_POLE) ** SUPPORT_CUTOFF

    singleton = singleton_activity()
    low_component_activity = singleton + SMALL_POST_GRAPH / target_graph
    if low_component_activity >= 1:
        raise SystemExit("low-support ledger: low component activity does not contract")
    # exp(x)-1 <= x/(1-x) and exp(x) <= 1/(1-x) for 0<=x<1.
    low_only = (
        target_graph
        * low_component_activity
        / (1 - low_component_activity)
    )

    tail_component_activity = (
        TAIL_POST_GRAPH / (tail_graph * cutoff_conversion)
    )
    if tail_component_activity >= 1:
        raise SystemExit("low-support ledger: tail component activity does not contract")
    with_tail = (
        target_graph
        * cutoff_conversion
        * tail_component_activity
        / (1 - tail_component_activity)
        / (1 - low_component_activity)
    )

    low_support_failure = low_only + with_tail
    print("support-at-most-106 three-band connected/disconnected ledger")
    print(f"singleton_pre_graph_log2={log2_fraction(singleton).hi}")
    print(
        "small_sizes_2_6_pre_graph_envelope_log2="
        f"{log2_fraction(SMALL_POST_GRAPH / target_graph).hi}"
    )
    print(f"low_component_activity_log2={log2_fraction(low_component_activity).hi}")
    print(f"low_only_post_graph_log2={log2_fraction(low_only).hi}")
    print(
        "tail_component_activity_at_1_20_log2="
        f"{log2_fraction(tail_component_activity).hi}"
    )
    print(f"with_tail_post_graph_log2={log2_fraction(with_tail).hi}")
    print(
        "low_support_failure_log2_upper="
        f"{log2_fraction(low_support_failure).hi}"
    )
    if low_support_failure > Fraction(1, 1 << 40):
        raise SystemExit("low-support ledger: failure exceeds 2^-40")
    print("low_support_failure_le_2^-40=PASS")
    print("status=EXACT_RATIONAL_LOW_SUPPORT_COMPONENT_LEDGER")


if __name__ == "__main__":
    main()
