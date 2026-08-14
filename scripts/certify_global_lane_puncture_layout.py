#!/usr/bin/env python3
"""Exact certificate for the global-lane puncture rule.

The three-band layout sends block (t,l) to tile labels

    (t, t + 9*l, t + 20*l) mod 256.

Sample one global lane L uniformly. Choose 128 distinct band-zero tiles
uniformly, and puncture block (t,L) in every chosen tile t. Choose the
punctured band-zero coordinate independently and uniformly from 42
coordinates. Finally, map the 128 graph coordinates by a uniform bijection
to the 128 holes.

This file exposes the sampler and certifies, with integer or rational
arithmetic:

* every fixed lane is a perfect matching in every band;
* every data block has puncture marginal 1/128;
* every data coordinate in band zero has marginal 1/(128*42);
* holes occupy distinct band-zero tiles and physical groups; and
* the retained lane can be regrouped into complete EBCH codewords in the
  one-conditioned-row exact-graph proof.

The global lane creates correlation between puncture events in different
tiles. This certificate claims uniform marginals, not independent lanes.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from fractions import Fraction


TILES = 256
LANES = 64
SLOPES = (0, 9, 20)
BAND_ZERO_COLUMNS = 42
HOLES = 128
BANDS = len(SLOPES)


@dataclass(frozen=True)
class PunctureHole:
    """One freed physical lane position."""

    band_zero_tile: int
    data_base: int
    data_lane: int
    band_zero_coordinate: int


@dataclass(frozen=True)
class GlobalLanePunctureSample:
    """One sample from the construction-facing puncture interface."""

    retained_lane: int
    holes: tuple[PunctureHole, ...]
    graph_coordinate_to_hole: tuple[int, ...]


def tile_label(base: int, lane: int, band: int) -> int:
    """Return the physical tile label for one codeword band."""

    if not 0 <= base < TILES:
        raise ValueError("base tile is out of range")
    if not 0 <= lane < LANES:
        raise ValueError("lane is out of range")
    if not 0 <= band < BANDS:
        raise ValueError("band is out of range")
    return (base + SLOPES[band] * lane) % TILES


def retained_base_for_tile(lane: int, band: int, tile: int) -> int:
    """Invert the fixed-lane matching in one band."""

    if not 0 <= lane < LANES:
        raise ValueError("lane is out of range")
    if not 0 <= band < BANDS:
        raise ValueError("band is out of range")
    if not 0 <= tile < TILES:
        raise ValueError("tile is out of range")
    return (tile - SLOPES[band] * lane) % TILES


def sample_global_lane_punctures(
    rng: random.Random,
) -> GlobalLanePunctureSample:
    """Sample the global lane, holes, coordinates, and graph bijection."""

    lane = rng.randrange(LANES)
    selected_tiles = sorted(rng.sample(range(TILES), HOLES))
    holes = tuple(
        PunctureHole(
            band_zero_tile=tile,
            data_base=tile,
            data_lane=lane,
            band_zero_coordinate=rng.randrange(BAND_ZERO_COLUMNS),
        )
        for tile in selected_tiles
    )
    graph_coordinate_to_hole = list(range(HOLES))
    rng.shuffle(graph_coordinate_to_hole)
    sample = GlobalLanePunctureSample(
        retained_lane=lane,
        holes=holes,
        graph_coordinate_to_hole=tuple(graph_coordinate_to_hole),
    )
    validate_sample(sample)
    return sample


def validate_sample(sample: GlobalLanePunctureSample) -> None:
    """Validate one sampler output without probabilistic assumptions."""

    if not 0 <= sample.retained_lane < LANES:
        raise AssertionError("sampled lane is out of range")
    if len(sample.holes) != HOLES:
        raise AssertionError("sample has the wrong hole count")
    if sorted(sample.graph_coordinate_to_hole) != list(range(HOLES)):
        raise AssertionError("graph-to-hole map is not a bijection")

    tiles: set[int] = set()
    physical_groups: set[tuple[int, int]] = set()
    for hole in sample.holes:
        if hole.data_lane != sample.retained_lane:
            raise AssertionError("hole does not use the global lane")
        if hole.data_base != hole.band_zero_tile:
            raise AssertionError("band-zero base/tile identity changed")
        if not 0 <= hole.band_zero_tile < TILES:
            raise AssertionError("hole tile is out of range")
        if not 0 <= hole.band_zero_coordinate < BAND_ZERO_COLUMNS:
            raise AssertionError("hole coordinate is out of range")
        if tile_label(hole.data_base, hole.data_lane, 0) != hole.band_zero_tile:
            raise AssertionError("punctured block is not retained in its hole tile")
        tiles.add(hole.band_zero_tile)
        physical_groups.add(
            (hole.band_zero_tile, hole.band_zero_coordinate)
        )
    if len(tiles) != HOLES:
        raise AssertionError("holes do not occupy distinct band-zero tiles")
    if len(physical_groups) != HOLES:
        raise AssertionError("holes do not occupy distinct physical groups")


def certify_perfect_matchings() -> None:
    """Exhaustively certify all 64 lane matchings and their inverses."""

    expected_tiles = set(range(TILES))
    for lane in range(LANES):
        for band in range(BANDS):
            labels = {
                tile_label(base, lane, band) for base in range(TILES)
            }
            if labels != expected_tiles:
                raise AssertionError(
                    f"lane {lane}, band {band} is not a perfect matching"
                )
            for tile in range(TILES):
                base = retained_base_for_tile(lane, band, tile)
                if tile_label(base, lane, band) != tile:
                    raise AssertionError("fixed-lane inverse map failed")


def certify_one_row_regrouping() -> None:
    """Certify the global reindexing used by the conditioned-row proof."""

    for lane in range(LANES):
        factors_by_physical_tile = {
            (
                retained_base_for_tile(lane, band, tile),
                band,
            )
            for band in range(BANDS)
            for tile in range(TILES)
        }
        factors_by_complete_codeword = {
            (base, band)
            for base in range(TILES)
            for band in range(BANDS)
        }
        if factors_by_physical_tile != factors_by_complete_codeword:
            raise AssertionError("retained factors do not regroup by codeword")

        # Every possible chosen band-zero tile punctures the retained row.
        for tile in range(TILES):
            base = retained_base_for_tile(lane, 0, tile)
            if base != tile or tile_label(base, lane, 0) != tile:
                raise AssertionError("hole does not puncture the retained row")


def certify_exact_marginals() -> None:
    """Check the sampler probabilities with exact rational arithmetic."""

    lane_probability = Fraction(1, LANES)
    tile_probability = Fraction(HOLES, TILES)
    block_probability = lane_probability * tile_probability
    coordinate_probability = block_probability / BAND_ZERO_COLUMNS
    if block_probability != Fraction(1, 128):
        raise AssertionError("per-block puncture marginal changed")
    if coordinate_probability != Fraction(1, 128 * BAND_ZERO_COLUMNS):
        raise AssertionError("per-coordinate puncture marginal changed")


def certify_graph_hole_distribution() -> None:
    """Check that graph support maps to a uniform subset of all 256 tiles."""

    # For a fixed graph support of size w, first choose 128 hole tiles and
    # then inject the w occupied graph coordinates uniformly into those holes.
    # The image probability of a fixed w-subset of all tiles is
    #
    # C(256-w,128-w) / C(256,128) / C(128,w) = 1 / C(256,w).
    for weight in range(HOLES + 1):
        containing = math.comb(TILES - weight, HOLES - weight)
        sampled = math.comb(TILES, HOLES)
        inside = math.comb(HOLES, weight)
        probability = Fraction(containing, sampled * inside)
        expected = Fraction(1, math.comb(TILES, weight))
        if probability != expected:
            raise AssertionError(
                f"graph-hole subset distribution failed at weight {weight}"
            )


def main() -> None:
    certify_perfect_matchings()
    certify_one_row_regrouping()
    certify_exact_marginals()
    certify_graph_hole_distribution()

    # Exercise the public sampler interface. Exact claims above do not rely
    # on this particular pseudorandom sample.
    sample = sample_global_lane_punctures(random.Random(0))
    print("global-lane puncture layout certificate")
    print(f"tiles={TILES} lanes={LANES} slopes={SLOPES}")
    print(f"holes={HOLES} band_zero_columns={BAND_ZERO_COLUMNS}")
    print("fixed_lane_perfect_matching_all_bands=PASS")
    print("one_conditioned_row_global_regrouping=PASS")
    print("per_block_puncture_marginal=1/128 PASS")
    print(f"per_band_zero_coordinate_marginal=1/{128 * BAND_ZERO_COLUMNS} PASS")
    print("distinct_band_zero_tiles_and_groups=PASS")
    print("graph_support_uniform_subset_of_256_tiles=PASS")
    print("lane_events_independent=NO_BY_CONSTRUCTION")
    print(f"sampler_interface_smoke_lane={sample.retained_lane} PASS")
    print("status=EXACT_GLOBAL_LANE_PUNCTURE_LAYOUT_CERTIFICATE")


if __name__ == "__main__":
    main()
