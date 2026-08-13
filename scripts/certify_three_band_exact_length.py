#!/usr/bin/env python3
"""Exact bookkeeping certificate for the power-of-two Riffle layout.

The exact-length construction uses the following zero-hot-path-cost sampling
rule.  Choose 128 of the 256 band-zero tiles uniformly.  In every chosen tile
choose one of its 64 data lanes uniformly and puncture a uniformly chosen
band-zero coordinate of that data word.  The 128 graph-word coordinates are
then mapped by a uniform bijection to the freed positions.

Consequences used by the proof are:

* graph coordinates occupy distinct band-zero tiles;
* every data block is punctured with marginal probability 1/128;
* conditioned on graph weight w, its occupied holes use a uniform w-subset
  of the 256 band-zero tiles and an independent uniform column in each tile;
* for data group support h, Maclaurin's inequality bounds the graph incremental
  moment by ``(q + (1-q) h/(256*42))**w``; and
* a punctured active block loses at most one group, so 21 rather than 22 is
  the deterministic per-block support floor.  The support-106 component range
  is consequently at most floor(64*106/21)=323 blocks.

All comparisons in this file use exact integers and rationals.  Logarithms
are printed only as diagnostics.
"""

from __future__ import annotations

import csv
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GRAPH_SPECTRUM = ROOT / "ebch128_graph24_spectrum.csv"
TILES = 256
LANES = 64
BAND_ZERO_COLUMNS = 42
PUNCTURES = 128
GRAPH_DIMENSION = 24
SUPPORT_CUTOFF = 106
DATA_DISTANCE = 22
PUNCTURED_DISTANCE = DATA_DISTANCE - 1
MAXIMUM_COMPONENT_BLOCKS = LANES * SUPPORT_CUTOFF // PUNCTURED_DISTANCE


def load_graph_spectrum() -> tuple[int, ...]:
    spectrum = [0] * 129
    with GRAPH_SPECTRUM.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            spectrum[int(row["weight"])] = int(row["count"])
    if sum(spectrum) != 1 << GRAPH_DIMENSION:
        raise SystemExit("exact length: graph spectrum cardinality mismatch")
    if spectrum[0] != 1:
        raise SystemExit("exact length: graph code does not have one zero word")
    if next(weight for weight, count in enumerate(spectrum) if weight and count) != 22:
        raise SystemExit("exact length: graph minimum distance changed")
    return tuple(spectrum)


def graph_factor_upper(
    pole: Fraction,
    *,
    data_support: int = SUPPORT_CUTOFF,
    spectrum: tuple[int, ...] | None = None,
) -> Fraction:
    """Return the normalized graph incremental-MGF upper bound.

    For occupied-group counts ``u_t`` in the 256 band-zero tiles, selecting a
    graph hole in tile t has incremental factor

        q + (1-q) u_t/42.

    Averaging products over a uniform w-subset of tiles and applying
    Maclaurin gives the w-th power of the arithmetic mean.  Replacing the
    band-zero occupied support by the total support is pessimistic.
    """

    if not 0 < pole < 1:
        raise ValueError("pole must lie in (0,1)")
    if not 0 <= data_support <= TILES * BAND_ZERO_COLUMNS:
        raise ValueError("invalid data support")
    if spectrum is None:
        spectrum = load_graph_spectrum()
    effective_pole = pole + (1 - pole) * Fraction(
        data_support, TILES * BAND_ZERO_COLUMNS
    )
    return sum(
        (count * effective_pole**weight for weight, count in enumerate(spectrum)),
        Fraction(0),
    ) / (1 << GRAPH_DIMENSION)


def puncture_selection_adjustment(
    full_factor: Fraction, punctured_factor: Fraction
) -> Fraction:
    """Turn one-block full/punctured bounds into a composable average bound.

    In a chosen band-zero tile only one of 64 lanes is punctured, and 128 of
    256 tiles are sampled without replacement.  If block i has factor ratio
    r_i, Maclaurin followed by ``1+sum x_i <= product(1+x_i)`` gives the
    per-block multiplier below.  Multiplying these adjustments over active
    blocks is therefore a valid upper bound despite the without-replacement
    dependence.
    """

    if full_factor <= 0 or punctured_factor < full_factor:
        raise ValueError("invalid full/punctured factor pair")
    ratio = punctured_factor / full_factor
    return (1 + (ratio - 1) / (TILES * LANES)) ** PUNCTURES


def worst_puncture_block_adjustment(pole: Fraction) -> Fraction:
    """Composable pointwise bound when no BCH-weight refinement is used."""

    if not 0 < pole < 1:
        raise ValueError("pole must lie in (0,1)")
    return (1 + (1 / pole - 1) / (TILES * LANES)) ** PUNCTURES


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    spectrum = load_graph_spectrum()
    # The graph factor is 2^-24 times a ratio whose nonzero-word excess is
    # tiny.  These dyadic envelopes cover the two poles used by the final
    # connected ledger and the target-pole diagnostics.
    checks = (
        (Fraction(181, 1000), 48),
        (Fraction(1, 20), 85),
    )
    print("exact three-band power-of-two layout certificate")
    print(
        f"tiles={TILES} lanes={LANES} band_zero_columns={BAND_ZERO_COLUMNS} "
        f"punctures={PUNCTURES}"
    )
    print("puncture_marginal=1/128")
    print("graph_holes_use_distinct_band_zero_tiles=BY_CONSTRUCTION")
    for pole, excess_bits in checks:
        factor = graph_factor_upper(pole, spectrum=spectrum)
        ratio = factor * (1 << GRAPH_DIMENSION)
        envelope = 1 + Fraction(1, 1 << excess_bits)
        if ratio > envelope:
            raise SystemExit(
                f"exact length: pole {pole} graph excess exceeds 2^-{excess_bits}"
            )
        print(
            f"pole={pole} graph_factor_log2_upper={log2_fraction(factor):.15f} "
            f"relative_excess_log2={log2_fraction(ratio - 1):.12f} "
            f"relative_excess_le_2^-{excess_bits}=PASS"
        )
    if MAXIMUM_COMPONENT_BLOCKS != 323:
        raise SystemExit("exact length: component-size cap changed")
    print(f"maximum_component_blocks={MAXIMUM_COMPONENT_BLOCKS} PASS")
    print("status=EXACT_RATIONAL_POWER_OF_TWO_BOOKKEEPING")


if __name__ == "__main__":
    main()
