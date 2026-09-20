#!/usr/bin/env python3
"""Exact one-block late-suffix lower bound for ideal Riffle g=8.

This diagnostic uses the fast sequential outer construction:

* 16,384 user-message blocks and one graph block;
* one componentwise-XOR parity block;
* one puncture in each of the first 256 outer blocks; and
* a uniform permutation of the N/8 eight-bit packets.

Only the 16,128 unpunctured user-message blocks enter the bound.  For one
such block, fix a nonzero local message whose graph syndrome is zero.  Its
data and parity BCH codewords are identical.  If all their active packets
land in the last floor(d/64) inner nodes, the emitted word has weight at most
d, regardless of the deterministic recursive inner map.

The script computes the exact packet-support enumerator of the binary
extended BCH [128,64,22] code.  For a subset S of the 16 byte packets, the
number of codewords supported in S is 2^dim(C_S).  Subset Moebius inversion
then gives the number with exact packet support S.

All ranks, counts, probabilities, and sums are exact integers or Fractions.
Only the printed base-two logarithms are approximate.  The result is a lower
bound on E[Z_d], not an upper-bound certificate.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction


N = 1 << 21
DISTANCE_THRESHOLD = 188_743
PACKET_BITS = 8
LOCAL_LENGTH = 128
LOCAL_DIMENSION = 64
GRAPH_CODIMENSION = 24
USER_BLOCKS = 1 << 14
PUNCTURED_USER_BLOCKS = 256
UNPUNCTURED_OUTER_BLOCKS = 16_386
GENERATOR = 0xF4845518B9582A1F


def binary_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def generator_rows() -> tuple[int, ...]:
    return tuple(
        (GENERATOR << row) | (1 << (LOCAL_LENGTH - 1))
        for row in range(LOCAL_DIMENSION)
    )


def packet_support_enumerator() -> tuple[int, ...]:
    packet_count = LOCAL_LENGTH // PACKET_BITS
    subset_count = 1 << packet_count
    full_mask = (1 << LOCAL_LENGTH) - 1
    packet_masks = tuple(
        ((1 << PACKET_BITS) - 1) << (PACKET_BITS * packet)
        for packet in range(packet_count)
    )
    rows = generator_rows()

    supported_counts = [0] * subset_count
    inside_masks = [0] * subset_count
    for subset in range(1, subset_count):
        low_bit = subset & -subset
        packet = low_bit.bit_length() - 1
        inside_masks[subset] = inside_masks[subset ^ low_bit] | packet_masks[packet]

    for subset, inside in enumerate(inside_masks):
        outside = full_mask ^ inside
        shortened_dimension = LOCAL_DIMENSION - binary_rank(
            [row & outside for row in rows]
        )
        supported_counts[subset] = 1 << shortened_dimension

    exact_counts = supported_counts[:]
    for packet in range(packet_count):
        bit = 1 << packet
        for subset in range(subset_count):
            if subset & bit:
                exact_counts[subset] -= exact_counts[subset ^ bit]

    enumerator = [0] * (packet_count + 1)
    for subset, count in enumerate(exact_counts):
        if count < 0:
            raise SystemExit("Riffle g=8 packet-support inversion became negative")
        enumerator[subset.bit_count()] += count
    if sum(enumerator) != 1 << LOCAL_DIMENSION:
        raise SystemExit("Riffle g=8 packet-support enumerator has wrong mass")
    return tuple(enumerator)


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        return -math.inf
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--unpunctured",
        action="store_true",
        help=(
            "use all 16,386 full outer BCH blocks, N=2^21+256, and "
            "d=floor(0.09N)"
        ),
    )
    args = parser.parse_args()

    code_length = UNPUNCTURED_OUTER_BLOCKS * LOCAL_LENGTH if args.unpunctured else N
    distance_threshold = (
        9 * code_length // 100 if args.unpunctured else DISTANCE_THRESHOLD
    )
    punctured_user_blocks = 0 if args.unpunctured else PUNCTURED_USER_BLOCKS
    packet_count = code_length // PACKET_BITS
    inner_nodes = code_length // 64
    late_nodes = distance_threshold // 64
    late_packets = late_nodes * (64 // PACKET_BITS)
    if 64 * late_nodes > distance_threshold:
        raise SystemExit("Riffle g=8 late suffix exceeds the distance threshold")

    enumerator = packet_support_enumerator()
    local_sum = Fraction(0)
    contributions: list[tuple[int, int, Fraction]] = []
    for support, count in enumerate(enumerator):
        if support == 0 or count == 0:
            continue
        late_probability = Fraction(
            math.comb(late_packets, 2 * support),
            math.comb(packet_count, 2 * support),
        )
        contribution = count * late_probability
        local_sum += contribution
        contributions.append((support, count, contribution))

    unpunctured_user_blocks = USER_BLOCKS - punctured_user_blocks
    expectation_lower = (
        Fraction(unpunctured_user_blocks, 1 << GRAPH_CODIMENSION) * local_sum
    )

    print("Riffle g=8 exact one-block late-suffix diagnostic")
    print(
        f"variant={'unpunctured' if args.unpunctured else 'power2-punctured'} "
        f"N={code_length} d={distance_threshold} inner_nodes={inner_nodes}"
    )
    print(f"packets={packet_count} late_nodes={late_nodes} late_packets={late_packets}")
    print(f"unpunctured_user_blocks={unpunctured_user_blocks}")
    print("packet_support,count,global_log2_contribution")
    for support, count, contribution in contributions:
        global_contribution = Fraction(
            unpunctured_user_blocks, 1 << GRAPH_CODIMENSION
        ) * contribution
        print(f"{support},{count},{log2_fraction(global_contribution):.12f}")
    print(f"packet_enumerator_mass={sum(enumerator)}")
    print(f"local_late_sum_log2={log2_fraction(local_sum):.12f}")
    print(f"expected_bad_messages_lower_log2={log2_fraction(expectation_lower):.12f}")
    print(
        "without_graph_charge_lower_log2="
        f"{log2_fraction(expectation_lower) + GRAPH_CODIMENSION:.12f}"
    )
    print("status=EXACT_RIFFLE_G8_ONE_BLOCK_LATE_SUFFIX_LOWER_BOUND")


if __name__ == "__main__":
    main()
