#!/usr/bin/env python3
"""Exact/outward certificate for fixed-band motifs of sizes two through six.

The script enumerates every triple of pairwise-compatible set partitions, so
pair capacity one is enforced exactly.  Patterns with the same component
count and multiset of one-codeword vertex types are aggregated first.  Their
integer multiplicities are then combined with 160-bit outward rooted moments,
exact cell-cap greedy bounds, the exact-length puncture adjustment, and the
graph incremental-MGF factor using rational arithmetic.
"""

from __future__ import annotations

import math
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path

from certify_three_band_exact_length import (
    BAND_ZERO_COLUMNS,
    graph_factor_upper,
    puncture_selection_adjustment,
)
from outward_log2 import log2_fraction
from probe_bch_three_band_cell_cap_motifs import cell_caps, cluster_sizes
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from probe_three_band_spreading import pairwise_compatible
from probe_two_band_spreading import set_partitions
from three_band_exact_moments import (
    factor_table_upper,
    greedy_factor_bound_upper,
)


BLOCKS = 1 << 14
LANES = 1 << 6
POLE = Fraction(181, 1000)
ROOT_BITS = 160
ENVELOPE_BITS = {2: 56, 3: 62, 4: 63, 5: 60, 6: 53}
PRODUCT_GUARD = 2.0**-48
SUM_GUARD = 2.0**-40
PROFILE_HELPER = Path(__file__).with_name("enumerate_three_band_small_profiles.exe")
EXPECTED_PARTITIONS = {2: 2, 3: 5, 4: 15, 5: 52, 6: 203}


def upward_float(value: Fraction) -> float:
    result = float(value)
    if Fraction.from_float(result) < value:
        result = math.nextafter(result, math.inf)
    if Fraction.from_float(result) < value:
        raise AssertionError("small motifs: upward float conversion failed")
    return result


def pattern_profiles(size: int) -> Counter[tuple[int, tuple[tuple[int, int, int], ...]]]:
    partitions = list(set_partitions(size))
    compatibility_masks: list[int] = []
    for left in range(len(partitions)):
        mask = 0
        for right in range(len(partitions)):
            if pairwise_compatible(partitions[left], partitions[right]):
                mask |= 1 << right
        compatibility_masks.append(mask)
    sizes = [cluster_sizes(partition) for partition in partitions]
    pairs = [
        (left, right)
        for left in range(size)
        for right in range(left + 1, size)
    ]
    edge_masks = []
    for partition in partitions:
        mask = 0
        for bit, (left, right) in enumerate(pairs):
            if partition[left] == partition[right]:
                mask |= 1 << bit
        edge_masks.append(mask)

    # Only 2^15 masks occur at size six.  Precomputing the component count
    # removes Python set/DFS work from the millions-of-triples hot loop.
    component_counts = [0] * (1 << len(pairs))
    for mask in range(len(component_counts)):
        parents = list(range(size))

        def find(vertex: int) -> int:
            while parents[vertex] != vertex:
                parents[vertex] = parents[parents[vertex]]
                vertex = parents[vertex]
            return vertex

        for bit, (left, right) in enumerate(pairs):
            if mask & (1 << bit):
                left_root = find(left)
                right_root = find(right)
                if left_root != right_root:
                    parents[right_root] = left_root
        component_counts[mask] = len({find(vertex) for vertex in range(size)})

    result: Counter[tuple[int, tuple[tuple[int, int, int], ...]]] = Counter()
    for first, first_partition in enumerate(partitions):
        second_mask = compatibility_masks[first]
        while second_mask:
            second_low = second_mask & -second_mask
            second = second_low.bit_length() - 1
            second_mask ^= second_low
            third_mask = compatibility_masks[first] & compatibility_masks[second]
            while third_mask:
                third_low = third_mask & -third_mask
                third = third_low.bit_length() - 1
                third_mask ^= third_low
                components = component_counts[
                    edge_masks[first] | edge_masks[second] | edge_masks[third]
                ]
                vertex_types = tuple(
                    sorted(
                        (
                            sizes[first][vertex],
                            sizes[second][vertex],
                            sizes[third][vertex],
                        )
                        for vertex in range(size)
                    )
                )
                result[components, vertex_types] += 1
    return result


def compiled_profiles() -> dict[
    int, Counter[tuple[int, tuple[tuple[int, int, int], ...]]]
]:
    if not PROFILE_HELPER.exists():
        raise SystemExit(
            "small motifs: compile enumerate_three_band_small_profiles.cpp first"
        )
    completed = subprocess.run(
        [str(PROFILE_HELPER)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = {size: Counter() for size in range(2, 7)}
    summaries: dict[int, tuple[int, int, int]] = {}
    status = False
    for line in completed.stdout.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "summary":
            size, partitions, triples_count, profiles_count = map(int, fields[1:])
            summaries[size] = (partitions, triples_count, profiles_count)
        elif fields[0] == "profile":
            size = int(fields[1])
            components = int(fields[2])
            count = int(fields[3])
            kinds = tuple(
                tuple(map(int, encoded.split(","))) for encoded in fields[4:]
            )
            if len(kinds) != size:
                raise SystemExit("small motifs: malformed helper profile")
            result[size][components, kinds] += count
        elif fields == ["status", "EXACT_PACKED_PROFILE_ENUMERATION"]:
            status = True
        else:
            raise SystemExit(f"small motifs: unrecognized helper row {line!r}")
    if not status or set(summaries) != set(range(2, 7)):
        raise SystemExit("small motifs: incomplete helper output")
    for size, profiles in result.items():
        partitions, triples_count, profiles_count = summaries[size]
        if partitions != EXPECTED_PARTITIONS[size]:
            raise SystemExit("small motifs: Bell-count regression")
        if len(profiles) != profiles_count or sum(profiles.values()) != triples_count:
            raise SystemExit("small motifs: helper count mismatch")
    return result


def main() -> None:
    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    tables = {
        (band, degree): factor_table_upper(
            BAND_SIZES[band], degree, POLE, root_bits=ROOT_BITS
        )
        for band in range(3)
        for degree in range(1, 7)
    }
    puncture_multiplier = tuple(
        1 + (1 / POLE - 1) * Fraction(weight, BAND_ZERO_COLUMNS)
        for weight in range(BAND_ZERO_COLUMNS + 1)
    )
    factors: dict[tuple[int, int, int], Fraction] = {}

    def adjusted_factor(kind: tuple[int, int, int]) -> Fraction:
        if kind not in factors:
            full_tables = tuple(tables[band, kind[band]] for band in range(3))
            full = greedy_factor_bound_upper(
                states, caps, diagonal_mass, by_total, full_tables
            )
            punctured_tables = (
                tuple(
                    full_tables[0][weight] * puncture_multiplier[weight]
                    for weight in range(BAND_ZERO_COLUMNS + 1)
                ),
                full_tables[1],
                full_tables[2],
            )
            punctured = greedy_factor_bound_upper(
                states, caps, diagonal_mass, by_total, punctured_tables
            )
            factors[kind] = full * puncture_selection_adjustment(full, punctured)
        return factors[kind]

    graph_factor = graph_factor_upper(POLE)
    profiles_by_size = compiled_profiles()
    totals: dict[int, Fraction] = {}
    for size in range(2, 7):
        profiles = profiles_by_size[size]
        terms: list[float] = []
        for (components, vertex_types), count in profiles.items():
            placement = Fraction(
                count * BLOCKS**components * LANES ** (size - components),
                math.factorial(size),
            )
            value = upward_float(placement) * upward_float(graph_factor)
            for kind in vertex_types:
                value *= upward_float(adjusted_factor(kind))
            # At most eight positive binary64 multiplications occur here.
            # This relative guard is over four times their naive gamma_n.
            value = math.nextafter(value * (1 + PRODUCT_GUARD), math.inf)
            terms.append(value)
        # CPython's fsum uses compensated partial sums.  Retain a much larger
        # explicit guard than one final ulp before converting back to an exact
        # Fraction for the threshold comparison.
        total_float = math.nextafter(
            math.fsum(terms) * (1 + SUM_GUARD), math.inf
        )
        total = Fraction.from_float(total_float)
        totals[size] = total
        print(
            f"size={size} profiles={len(profiles)} "
            f"log2_upper={log2_fraction(total).hi} "
            f"target=2^-{ENVELOPE_BITS[size]}"
        )
        if total > Fraction(1, 1 << ENVELOPE_BITS[size]):
            raise SystemExit(
                f"small motifs: size {size} exceeds 2^-{ENVELOPE_BITS[size]}"
            )
        print(f"size_{size}_le_2^-{ENVELOPE_BITS[size]}=PASS")
    combined_float = math.nextafter(
        math.fsum(upward_float(value) for value in totals.values())
        * (1 + SUM_GUARD),
        math.inf,
    )
    combined = Fraction.from_float(combined_float)
    print(f"sizes_2_6_sum_log2_upper={log2_fraction(combined).hi}")
    print("status=EXACT_COMBINATORICS_DIRECTED_OUTWARD_PROFILE_SUM")


if __name__ == "__main__":
    main()
