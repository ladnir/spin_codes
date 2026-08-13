#!/usr/bin/env python3
"""Two-band outer spreading moment for low active-block counts.

For a fixed set of active EBCH data blocks, randomize the 64/64 coordinate
split and tile assignment in each band.  Excluding zero and all-ones, the
exact joint band-weight law has a small density relative to the product of
its marginals.  This decouples the coordinate supports while retaining the
tile partitions.  The one-band support moments use exact hypergeometric union
transitions.

For small active-block counts the script also enumerates both set partitions
and deletes every pair having a joint tile cell larger than two.  This models
the proposed capacity-two cross-band incidence.  Long-double arithmetic makes
the final result a design diagnostic; all combinatorics before rounding are
exact.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from punctured_ebch_outer import full_spectrum


BLOCKS = 16384
TILES = 256
BAND_COLUMNS = 64
LOCAL_CARDINALITY = 1 << 64
NONTRIVIAL_CARDINALITY = LOCAL_CARDINALITY - 2
GRAPH_CODIMENSION = 24


def band_weight_distribution() -> np.ndarray:
    """Marginal band weight for a uniform nonzero, non-all-ones word."""

    spectrum = full_spectrum()
    result = np.zeros(BAND_COLUMNS + 1, dtype=np.longdouble)
    denominator_messages = NONTRIVIAL_CARDINALITY
    for total_weight, count in enumerate(spectrum):
        if not 0 < total_weight < 128 or not count:
            continue
        denominator = math.comb(2 * BAND_COLUMNS, total_weight)
        for band_weight in range(
            max(0, total_weight - BAND_COLUMNS),
            min(BAND_COLUMNS, total_weight) + 1,
        ):
            ways = math.comb(BAND_COLUMNS, band_weight) * math.comb(
                BAND_COLUMNS, total_weight - band_weight
            )
            result[band_weight] += (
                np.longdouble(count)
                / denominator_messages
                * np.longdouble(ways)
                / denominator
            )
    if abs(float(np.sum(result)) - 1.0) > 1e-15:
        raise SystemExit("two-band spreading: band weight row does not sum to one")
    return result


def tilted_band_distribution(tilt: float) -> tuple[np.ndarray, float]:
    """Tilted marginal and max joint density relative to its product."""

    spectrum = full_spectrum()
    marginal = band_weight_distribution()
    normalizer = sum(
        float(marginal[weight]) * tilt**weight
        for weight in range(BAND_COLUMNS + 1)
    )
    tilted = np.array(
        [marginal[weight] * tilt**weight / normalizer for weight in range(BAND_COLUMNS + 1)],
        dtype=np.longdouble,
    )
    maximum = 0.0
    for left in range(BAND_COLUMNS + 1):
        for right in range(BAND_COLUMNS + 1):
            total = left + right
            if not 0 < total < 128 or not spectrum[total]:
                continue
            joint = (
                spectrum[total]
                / NONTRIVIAL_CARDINALITY
                * math.comb(BAND_COLUMNS, left)
                * math.comb(BAND_COLUMNS, right)
                / math.comb(2 * BAND_COLUMNS, total)
            )
            maximum = max(
                maximum, joint / float(tilted[left] * tilted[right])
            )
    return tilted, maximum


def cluster_moments(
    s_max: int, pole: float, weights: np.ndarray | None = None
) -> np.ndarray:
    """Moment of the union of k independent band supports in one tile."""

    weights = band_weight_distribution() if weights is None else weights
    distribution = np.zeros(BAND_COLUMNS + 1, dtype=np.longdouble)
    distribution[0] = 1.0
    moments = np.zeros(s_max + 1, dtype=np.longdouble)
    moments[0] = 1.0
    powers = np.array(
        [np.longdouble(pole) ** union for union in range(BAND_COLUMNS + 1)]
    )
    for blocks in range(1, s_max + 1):
        next_distribution = np.zeros_like(distribution)
        for union, prior in enumerate(distribution):
            if not prior:
                continue
            for weight, weight_probability in enumerate(weights):
                if not weight_probability:
                    continue
                denominator = math.comb(BAND_COLUMNS, weight)
                for intersection in range(
                    max(0, union + weight - BAND_COLUMNS),
                    min(union, weight) + 1,
                ):
                    probability = (
                        math.comb(union, intersection)
                        * math.comb(
                            BAND_COLUMNS - union, weight - intersection
                        )
                        / denominator
                    )
                    next_distribution[union + weight - intersection] += (
                        prior * weight_probability * probability
                    )
        distribution = next_distribution
        moments[blocks] = np.sum(distribution * powers)
    return moments


def tiled_moments(s_max: int, cluster: np.ndarray) -> np.ndarray:
    """Moments for s labelled blocks assigned iid to 256 tiles."""

    tile_polynomial = np.array(
        [cluster[size] / math.factorial(size) for size in range(s_max + 1)],
        dtype=np.longdouble,
    )
    coefficients = np.zeros(s_max + 1, dtype=np.longdouble)
    coefficients[0] = 1.0
    for _ in range(TILES):
        coefficients = np.convolve(coefficients, tile_polynomial)[: s_max + 1]
    result = np.zeros(s_max + 1, dtype=np.longdouble)
    for blocks in range(s_max + 1):
        result[blocks] = (
            math.factorial(blocks) * coefficients[blocks] / TILES**blocks
        )
    return result


def set_partitions(size: int):
    labels = [0] * size

    def visit(index: int, maximum: int):
        if index == size:
            yield tuple(labels)
            return
        for label in range(maximum + 2):
            labels[index] = label
            yield from visit(index + 1, max(maximum, label))

    if size == 0:
        yield ()
    else:
        labels[0] = 0
        yield from visit(1, 0)


def capped_pair_moment(
    blocks: int, cluster: np.ndarray, cell_capacity: int
) -> tuple[float, tuple[int, ...], tuple[int, ...], float]:
    """Independent partition moment with bounded joint-cell occupancy."""

    rows = []
    for labels in set_partitions(blocks):
        components = max(labels) + 1
        sizes = [0] * components
        for label in labels:
            sizes[label] += 1
        support_moment = np.longdouble(1)
        for size in sizes:
            support_moment *= cluster[size]
        tile_probability = np.longdouble(1)
        for offset in range(components):
            tile_probability *= TILES - offset
        tile_probability /= np.longdouble(TILES) ** blocks
        rows.append((labels, float(tile_probability * support_moment)))

    total = 0.0
    dominant = (-1.0, (), ())
    for left_labels, left_value in rows:
        for right_labels, right_value in rows:
            cells: dict[tuple[int, int], int] = {}
            valid = True
            for left, right in zip(left_labels, right_labels):
                key = (left, right)
                count = cells.get(key, 0) + 1
                if count > cell_capacity:
                    valid = False
                    break
                cells[key] = count
            if valid:
                value = left_value * right_value
                total += value
                if value > dominant[0]:
                    dominant = (value, left_labels, right_labels)
    return total, dominant[1], dominant[2], dominant[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s-max", type=int, default=23)
    parser.add_argument(
        "--group-pole",
        type=float,
        default=0.181,
        help="q in the target q^G group-support envelope",
    )
    parser.add_argument("--cap2-max", type=int, default=8)
    parser.add_argument("--cell-capacity", type=int, default=2)
    parser.add_argument(
        "--tilts",
        default="1,1.1,1.2,1.35,1.5,1.75,2,2.5,3,4",
        help="comma-separated exponential tilts for the band marginal",
    )
    args = parser.parse_args()
    if not 1 <= args.s_max <= 32 or not 0 < args.group_pole < 1:
        raise SystemExit("two-band spreading: invalid parameters")

    tilts = [float(value) for value in args.tilts.split(",")]
    candidates = []
    for tilt in tilts:
        weights, density = tilted_band_distribution(tilt)
        cluster = cluster_moments(args.s_max, args.group_pole, weights)
        tiled = tiled_moments(args.s_max, cluster)
        candidates.append((tilt, density, cluster, tiled))
    print("two-band striped support moment")
    print(f"group_pole={args.group_pole:.12f}")
    for tilt, density, _cluster, _tiled in candidates:
        print(
            f"tilt={tilt:.8f} split_density_log2={math.log2(density):.12f}"
        )
    rows = []
    for active_blocks in range(1, args.s_max + 1):
        # C(M,s) times the without-replacement-to-iid correction M^s/(M)_s
        # is exactly M^s/s!.  This safely adds repeated fictitious slots.
        uncapped_candidates = []
        for tilt, density, _cluster, tiled in candidates:
            value = (
                active_blocks * (
                    math.log2(BLOCKS)
                    + math.log2(NONTRIVIAL_CARDINALITY)
                    + math.log2(density)
                )
                - math.lgamma(active_blocks + 1) / math.log(2)
                - GRAPH_CODIMENSION
                + 2 * math.log2(float(tiled[active_blocks]))
            )
            uncapped_candidates.append((value, tilt))
        log2_expected, best_uncapped_tilt = min(uncapped_candidates)
        rows.append(log2_expected)
        print(
            f"active_blocks={active_blocks:2d} "
            f"outer_weighted_log2={log2_expected:.12f} "
            f"best_tilt={best_uncapped_tilt:.8f}"
        )
        if active_blocks <= args.cap2_max:
            capped_candidates = []
            for tilt, density, cluster, _tiled in candidates:
                capped, left_labels, right_labels, dominant_value = capped_pair_moment(
                    active_blocks, cluster, args.cell_capacity
                )
                value = (
                    active_blocks * (
                        math.log2(BLOCKS)
                        + math.log2(NONTRIVIAL_CARDINALITY)
                        + math.log2(density)
                    )
                    - math.lgamma(active_blocks + 1) / math.log(2)
                    - GRAPH_CODIMENSION
                    + math.log2(capped)
                )
                capped_candidates.append(
                    (
                        value,
                        tilt,
                        left_labels,
                        right_labels,
                        dominant_value,
                    )
                )
            (
                capped_log2,
                best_capped_tilt,
                best_left_labels,
                best_right_labels,
                _dominant_value,
            ) = min(capped_candidates)
            print(
                f"active_blocks={active_blocks:2d} "
                f"cap{args.cell_capacity}_outer_weighted_log2="
                f"{capped_log2:.12f} best_tilt={best_capped_tilt:.8f} "
                f"dominant_partitions={best_left_labels}/{best_right_labels}"
            )
    maximum = max(rows)
    total = maximum + math.log2(sum(2 ** (value - maximum) for value in rows))
    print(f"all_s_log2={total:.12f}")
    print("status=DIAGNOSTIC_NEEDS_CERTIFIED_GROUP_ENVELOPE_AND_OUTWARD_ARITHMETIC")


if __name__ == "__main__":
    main()
