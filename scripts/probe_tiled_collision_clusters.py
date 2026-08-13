#!/usr/bin/env python3
"""Tiled-collision diagnostic retaining all within-tile correlations.

Active data blocks are partitioned by physical tile.  Within a tile, support
unions are propagated with the exact hypergeometric law.  Marked tiles use
their actual 63 full lanes plus one punctured lane and a 127-coordinate
non-hole universe; unmarked tiles use 64 full lanes and 128 coordinates.

The injective assignment of logical blocks to physical lanes is upper-bounded
by independent lane draws with the exact factor M^s/(M)_s.  Repeated draws
only add nonnegative fictitious configurations.  Distinct tile labels are
then counted exactly with (128)_a (128)_b for a unmarked and b marked tile
components.  The local PA1 bounds are exact rational certificates, while the
final placement sums here remain floating diagnostics.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_striped_group_low_ledger import (
    DATA_BLOCKS,
    GRAPH_CODIMENSION,
    episode_parameters,
    full_spectrum,
    group_episode_log2,
    punctured_spectrum,
)


B = 64
TILES_PER_KIND = 128
TILES = 256
LOCAL_CARDINALITY = 1 << 64


def hypergeom(population: int, marked: int, draws: int, hit: int) -> float:
    if not max(0, draws - (population - marked)) <= hit <= min(marked, draws):
        return 0.0
    return (
        math.comb(marked, hit)
        * math.comb(population - marked, draws - hit)
        / math.comb(population, draws)
    )


def cluster_distributions(s_max: int, marked_tile: bool) -> list[np.ndarray]:
    """Q[k,u]: normalized mass of k active words with union weight u."""

    full = full_spectrum()
    punctured = punctured_spectrum()
    # State is (non-hole union size, hole occupied).  Unmarked tiles still use
    # the same representation, with every full coordinate choosing the
    # distinguished group according to its exact w/128 marginal.
    state = np.zeros((128, 2), dtype=np.longdouble)
    state[0, 0] = 1.0
    result = [np.zeros(129, dtype=np.longdouble) for _ in range(s_max + 1)]
    result[0][0] = 1.0

    choices: list[tuple[int, float, bool]] = []
    full_lane_probability = 63 / 64 if marked_tile else 1.0
    punctured_lane_probability = 1 / 64 if marked_tile else 0.0
    for weight, count in enumerate(full):
        if not weight or not count:
            continue
        mass = full_lane_probability * count / LOCAL_CARDINALITY
        if weight:
            choices.append((weight - 1, mass * weight / 128, True))
        if weight < 128:
            choices.append((weight, mass * (128 - weight) / 128, False))
    if punctured_lane_probability:
        for weight, count in enumerate(punctured):
            if weight and count:
                choices.append(
                    (
                        weight,
                        punctured_lane_probability * count / LOCAL_CARDINALITY,
                        False,
                    )
                )

    for blocks in range(1, s_max + 1):
        next_state = np.zeros_like(state)
        for occupied in range(128):
            for hole in range(2):
                prior = state[occupied, hole]
                if not prior:
                    continue
                for nonhole_weight, mass, uses_hole in choices:
                    for intersection in range(
                        max(0, nonhole_weight - (127 - occupied)),
                        min(occupied, nonhole_weight) + 1,
                    ):
                        probability = hypergeom(
                            127, occupied, nonhole_weight, intersection
                        )
                        new_occupied = occupied + nonhole_weight - intersection
                        next_state[
                            new_occupied, int(bool(hole) or uses_hole)
                        ] += (
                            prior * mass * probability
                        )
        state = next_state
        for occupied in range(128):
            result[blocks][occupied] += state[occupied, 0]
            result[blocks][occupied + 1] += state[occupied, 1]
    return result


def cluster_power(
    clusters: list[np.ndarray], component_max: int, s_max: int, g_max: int
) -> list[np.ndarray]:
    """P[a][s,g] = coefficient of C(x,z)^a for ordered components."""

    powers = [np.zeros((s_max + 1, g_max + 1), dtype=np.longdouble) for _ in range(component_max + 1)]
    powers[0][0, 0] = 1.0
    for components in range(1, component_max + 1):
        previous = powers[components - 1]
        current = powers[components]
        for used in range(s_max + 1):
            for group_support in range(g_max + 1):
                base = previous[used, group_support]
                if not base:
                    continue
                for size in range(1, s_max - used + 1):
                    scale = base / math.factorial(size)
                    row = clusters[size]
                    limit = min(128, g_max - group_support)
                    current[used + size, group_support : group_support + limit + 1] += (
                        scale * row[: limit + 1]
                    )
    return powers


def falling(n: int, k: int) -> int:
    value = 1
    for offset in range(k):
        value *= n - offset
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s-max", type=int, default=8)
    parser.add_argument("--g-max", type=int, default=106)
    args = parser.parse_args()
    if not 1 <= args.s_max <= 16 or not 21 <= args.g_max <= 128:
        raise SystemExit("tile cluster probe: invalid truncation")

    parameters = episode_parameters(True)
    inner = {
        group_support: group_episode_log2(group_support, parameters)
        for group_support in range(21, args.g_max + 1)
    }
    unmarked = cluster_distributions(args.s_max, False)
    marked = cluster_distributions(args.s_max, True)
    unmarked_powers = cluster_power(
        unmarked, args.s_max, args.s_max, args.g_max
    )
    marked_powers = cluster_power(marked, args.s_max, args.s_max, args.g_max)

    print("tiled collision cluster diagnostic")
    print(f"s_max={args.s_max} g_max={args.g_max}")
    all_terms: list[float] = []
    for active_blocks in range(1, args.s_max + 1):
        terms: list[float] = []
        dominant = (-math.inf, "")
        common_log2 = (
            active_blocks * (math.log2(DATA_BLOCKS) + 64 - math.log2(TILES))
            - GRAPH_CODIMENSION
        )
        for unmarked_components in range(active_blocks + 1):
            for marked_components in range(active_blocks - unmarked_components + 1):
                components = unmarked_components + marked_components
                if not components:
                    continue
                label_factor = (
                    falling(TILES_PER_KIND, unmarked_components)
                    * falling(TILES_PER_KIND, marked_components)
                    / (
                        math.factorial(unmarked_components)
                        * math.factorial(marked_components)
                    )
                )
                if not label_factor:
                    continue
                left = unmarked_powers[unmarked_components]
                right = marked_powers[marked_components]
                for left_blocks in range(active_blocks + 1):
                    right_blocks = active_blocks - left_blocks
                    left_row = left[left_blocks]
                    right_row = right[right_blocks]
                    for left_groups in np.flatnonzero(left_row):
                        limit = args.g_max - int(left_groups)
                        if limit < 0:
                            continue
                        for right_groups in np.flatnonzero(right_row[: limit + 1]):
                            total_groups = int(left_groups + right_groups)
                            if total_groups < 21:
                                continue
                            coefficient = (
                                float(left_row[left_groups])
                                * float(right_row[right_groups])
                                * label_factor
                            )
                            if coefficient:
                                value = (
                                    common_log2
                                    + math.log2(coefficient)
                                    + inner[total_groups]
                                )
                                terms.append(value)
                                if value > dominant[0]:
                                    dominant = (
                                        value,
                                        f"U={unmarked_components} M={marked_components} "
                                        f"left_blocks={left_blocks} g={total_groups}",
                                    )
        maximum = max(terms)
        total = maximum + math.log2(sum(2 ** (value - maximum) for value in terms))
        all_terms.append(total)
        print(
            f"active_blocks={active_blocks} log2_upper_partial={total:.12f} "
            f"dominant={dominant[0]:.12f}[{dominant[1]}]"
        )
    maximum = max(all_terms)
    total = maximum + math.log2(sum(2 ** (value - maximum) for value in all_terms))
    print(f"all_s_le_{args.s_max}_g_le_{args.g_max}_log2_upper={total:.12f}")
    print("status=DIAGNOSTIC_PLACEMENT_ARITHMETIC_NOT_YET_OUTWARD")


if __name__ == "__main__":
    main()
