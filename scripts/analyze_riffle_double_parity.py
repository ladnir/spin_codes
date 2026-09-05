#!/usr/bin/env python3
"""Exact low-support diagnostics for Riffle Double-Parity g=8.

The script analyzes the unpunctured length N=(2^14+2)*128.  It compares
three proof-relevant families:

* one-block messages in Graph-24 and Double-Parity;
* equal two-block messages in Graph-24; and
* parity-zero messages on three consecutive Double-Parity blocks.

The consecutive triple has local messages (x*t, (1+x)*t, t), up to a common
nonzero field multiplier.  The script exactly enumerates every triple whose
total packet support is at most 26.  Such a triple has at least one BCH word
with packet support at most eight, so shortened-code enumeration is complete
for that range.

All combinatorial values are exact integers or Fractions.  Printed logarithms
are diagnostic base-two approximations.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter
from fractions import Fraction

import analyze_riffle_g8_late_suffix as base


FIELD_REDUCTION = 0x1B
FIELD_ORDER = (1 << 64) - 1
USER_BLOCKS = 1 << 14
MAX_EXACT_LOCAL_PACKET_SUPPORT = 8


def field_multiply_x(value: int) -> int:
    return ((value << 1) & FIELD_ORDER) ^ (FIELD_REDUCTION if value >> 63 else 0)


def field_multiply(left: int, right: int) -> int:
    result = 0
    for _ in range(64):
        if right & 1:
            result ^= left
        right >>= 1
        left = field_multiply_x(left)
    return result


def field_power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = field_multiply(result, value)
        value = field_multiply(value, value)
        exponent >>= 1
    return result


def field_inverse(value: int) -> int:
    if value == 0:
        raise ValueError("zero has no field inverse")
    return field_power(value, FIELD_ORDER - 1)


GENERATOR_ROWS = base.generator_rows()


def encode_message(message: int) -> int:
    codeword = 0
    while message:
        bit = (message & -message).bit_length() - 1
        codeword ^= GENERATOR_ROWS[bit]
        message &= message - 1
    return codeword


def packet_mask(codeword: int) -> int:
    result = 0
    for packet in range(16):
        if (codeword >> (8 * packet)) & 0xFF:
            result |= 1 << packet
    return result


def output_constraint_rows() -> tuple[int, ...]:
    rows = []
    for output_bit in range(128):
        row = 0
        for message_bit, generator_row in enumerate(GENERATOR_ROWS):
            if (generator_row >> output_bit) & 1:
                row |= 1 << message_bit
        rows.append(row)
    return tuple(rows)


CONSTRAINT_ROWS = output_constraint_rows()


def nullspace_basis(rows: list[int]) -> tuple[int, ...]:
    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
                continue
            pivots[pivot] = row
            break

    # Each pivot is the highest bit in its row.  Remove lower pivot columns
    # before constructing the nullspace vectors.
    for pivot in sorted(pivots):
        row = pivots[pivot]
        for lower in sorted(column for column in pivots if column < pivot):
            if (row >> lower) & 1:
                row ^= pivots[lower]
        pivots[pivot] = row

    free = [column for column in range(64) if column not in pivots]
    basis = []
    for column in free:
        vector = 1 << column
        for pivot, row in pivots.items():
            if (row >> column) & 1:
                vector |= 1 << pivot
        basis.append(vector)

    for vector in basis:
        if any((row & vector).bit_count() & 1 for row in rows):
            raise RuntimeError("shortened BCH nullspace verification failed")
    return tuple(basis)


def exact_low_support_messages(maximum_support: int) -> set[int]:
    messages: set[int] = set()
    for support in range(1, maximum_support + 1):
        for packets in itertools.combinations(range(16), support):
            allowed = sum(1 << packet for packet in packets)
            constraints = [
                CONSTRAINT_ROWS[bit]
                for bit in range(128)
                if not (allowed >> (bit // 8)) & 1
            ]
            basis = nullspace_basis(constraints)
            for selector in range(1, 1 << len(basis)):
                message = 0
                for index, vector in enumerate(basis):
                    if (selector >> index) & 1:
                        message ^= vector
                if packet_mask(encode_message(message)) == allowed:
                    messages.add(message)
    return messages


def late_probability(packet_count: int) -> Fraction:
    code_length = base.UNPUNCTURED_OUTER_BLOCKS * base.LOCAL_LENGTH
    distance = 9 * code_length // 100
    all_packets = code_length // base.PACKET_BITS
    late_packets = (distance // 64) * (64 // base.PACKET_BITS)
    return Fraction(
        math.comb(late_packets, packet_count),
        math.comb(all_packets, packet_count),
    )


def one_block_bounds(enumerator: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    lower = Fraction()
    upper = Fraction()
    for support, count in enumerate(enumerator):
        if support == 0 or count == 0:
            continue
        lower += count * late_probability(2 * support + 16)
        upper += count * late_probability(2 * support + 7)
    return USER_BLOCKS * lower, USER_BLOCKS * upper


def mds_block_weight_count(block_weight: int) -> int:
    field_size = 1 << 64
    block_length = USER_BLOCKS + 2
    distance = 3
    if block_weight < distance or block_weight > block_length:
        return 0
    full_support = (field_size - 1) * sum(
        (-1) ** index
        * math.comb(block_weight - 1, index)
        * field_size ** (block_weight - distance - index)
        for index in range(block_weight - distance + 1)
    )
    return math.comb(block_length, block_weight) * full_support


def mds_terminal_lower_bound(block_weight: int) -> Fraction:
    # Every active BCH block occupies at most 16 eight-bit packets.  Requiring
    # all 16*block_weight possible packets to land in the terminal set is a
    # subset of the bad event for every codeword of this block weight.
    return mds_block_weight_count(block_weight) * late_probability(16 * block_weight)


def randomized_one_block_exact(enumerator: tuple[int, ...]) -> Fraction:
    result = Fraction()
    for support0, count0 in enumerate(enumerator):
        if support0 == 0 or count0 == 0:
            continue
        for support1, count1 in enumerate(enumerator):
            if support1 == 0 or count1 == 0:
                continue
            result += (
                Fraction(count0 * count1, FIELD_ORDER)
                * late_probability(2 * support0 + support1)
            )
    return USER_BLOCKS * result


def graph_families(enumerator: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    local = Fraction()
    for support, count in enumerate(enumerator):
        if support and count:
            local += count * late_probability(2 * support)
    one = Fraction(USER_BLOCKS, 1 << base.GRAPH_CODIMENSION) * local
    pair = Fraction(math.comb(USER_BLOCKS, 2), 1 << base.GRAPH_CODIMENSION) * local
    return one, pair


def graph_power2_unpunctured_pair(enumerator: tuple[int, ...]) -> Fraction:
    all_packets = base.N // base.PACKET_BITS
    late_packets = (
        base.DISTANCE_THRESHOLD // 64
    ) * (64 // base.PACKET_BITS)
    local = Fraction()
    for support, count in enumerate(enumerator):
        if support and count:
            local += count * Fraction(
                math.comb(late_packets, 2 * support),
                math.comb(all_packets, 2 * support),
            )
    blocks = USER_BLOCKS - base.PUNCTURED_USER_BLOCKS
    return Fraction(math.comb(blocks, 2), 1 << base.GRAPH_CODIMENSION) * local


def consecutive_triple_low_tail(low_messages: set[int]) -> tuple[Counter[int], Fraction]:
    inverse_x = field_inverse(2)
    inverse_one_plus_x = field_inverse(3)
    parameters: set[int] = set()
    for message in low_messages:
        parameters.add(message)
        parameters.add(field_multiply(inverse_x, message))
        parameters.add(field_multiply(inverse_one_plus_x, message))

    enumerator: Counter[int] = Counter()
    contribution = Fraction()
    for parameter in parameters:
        shifted = field_multiply_x(parameter)
        messages = (shifted, shifted ^ parameter, parameter)
        total_support = sum(
            packet_mask(encode_message(message)).bit_count()
            for message in messages
        )
        if total_support <= 26:
            enumerator[total_support] += 1
            contribution += late_probability(total_support)

    # There are USER_BLOCKS-2 consecutive index triples.  Their message sets
    # are disjoint because their outer-block supports differ.
    return enumerator, (USER_BLOCKS - 2) * contribution


def fixed_gap_triple_low_tail(
    low_messages: set[int], gap1: int, gap2: int
) -> tuple[Counter[int], Fraction]:
    if not 0 < gap1 < gap2 < USER_BLOCKS:
        raise ValueError("triple gaps must satisfy 0 < gap1 < gap2 < B")
    powers = [1]
    for _ in range(1, gap2 + 1):
        powers.append(field_multiply_x(powers[-1]))
    coefficients = (
        powers[gap1] ^ powers[gap2],
        powers[gap2] ^ 1,
        1 ^ powers[gap1],
    )
    inverses = tuple(field_inverse(value) for value in coefficients)
    parameters = {
        field_multiply(inverse, message)
        for inverse in inverses
        for message in low_messages
    }
    enumerator: Counter[int] = Counter()
    contribution = Fraction()
    for parameter in parameters:
        total_support = sum(
            packet_mask(encode_message(field_multiply(coefficient, parameter))).bit_count()
            for coefficient in coefficients
        )
        if total_support <= 26:
            enumerator[total_support] += 1
            contribution += late_probability(total_support)
    return enumerator, (USER_BLOCKS - gap2) * contribution


def all_triple_support_21(low_messages: set[int]) -> tuple[int, Fraction]:
    support_seven = {
        message
        for message in low_messages
        if packet_mask(encode_message(message)).bit_count() == 7
    }
    additive_triples = {
        tuple(sorted((left, right, left ^ right)))
        for left in support_seven
        for right in support_seven
        if left != right and (left ^ right) in support_seven
    }

    powers = [1]
    for _ in range(1, USER_BLOCKS):
        powers.append(field_multiply_x(powers[-1]))
    discrete_logs = {value: exponent for exponent, value in enumerate(powers)}

    codeword_count = 0
    for triple in additive_triples:
        for value0, value1, value2 in itertools.permutations(triple):
            inverse2 = field_inverse(value2)
            for gap1, alpha1 in enumerate(powers[1:-1], 1):
                alpha2 = field_multiply(
                    value0 ^ field_multiply(alpha1, value1), inverse2
                )
                gap2 = discrete_logs.get(alpha2)
                if gap2 is not None and gap2 > gap1:
                    codeword_count += USER_BLOCKS - gap2

    return codeword_count, codeword_count * late_probability(21)


def sample_consecutive_triple(
    sample_count: int, batch_size: int = 100_000
) -> tuple[Counter[int], float]:
    import numpy as np

    generator_low = np.asarray(
        [row & FIELD_ORDER for row in GENERATOR_ROWS], dtype=np.uint64
    )
    generator_high = np.asarray(
        [row >> 64 for row in GENERATOR_ROWS], dtype=np.uint64
    )
    maximum = np.uint64(FIELD_ORDER)

    def support(messages: np.ndarray) -> np.ndarray:
        low = np.zeros(messages.size, dtype=np.uint64)
        high = np.zeros(messages.size, dtype=np.uint64)
        for bit in range(64):
            active = (messages >> np.uint64(bit)) & np.uint64(1)
            mask = np.where(active != 0, maximum, np.uint64(0))
            low ^= mask & generator_low[bit]
            high ^= mask & generator_high[bit]
        return (
            np.count_nonzero(low.view(np.uint8).reshape(-1, 8), axis=1)
            + np.count_nonzero(high.view(np.uint8).reshape(-1, 8), axis=1)
        )

    probabilities = {
        packets: float(late_probability(packets)) for packets in range(21, 49)
    }
    rng = np.random.default_rng(0xD0B1E)
    histogram: Counter[int] = Counter()
    probability_sum = 0.0
    completed = 0
    while completed < sample_count:
        count = min(batch_size, sample_count - completed)
        parameters = rng.integers(0, 1 << 64, size=count, dtype=np.uint64)
        shifted = (parameters << np.uint64(1)) ^ np.where(
            parameters >> np.uint64(63), np.uint64(FIELD_REDUCTION), np.uint64(0)
        )
        totals = support(parameters) + support(shifted) + support(parameters ^ shifted)
        values, counts = np.unique(totals, return_counts=True)
        for packets, frequency in zip(values.tolist(), counts.tolist()):
            histogram[int(packets)] += int(frequency)
            if packets >= 27:
                probability_sum += frequency * probabilities[int(packets)]
        completed += count
    return histogram, probability_sum / sample_count


def sample_support_three_category(
    category: str, sample_count: int, batch_size: int = 100_000
) -> tuple[Counter[int], float]:
    import numpy as np

    generator_low = np.asarray(
        [row & FIELD_ORDER for row in GENERATOR_ROWS], dtype=np.uint64
    )
    generator_high = np.asarray(
        [row >> 64 for row in GENERATOR_ROWS], dtype=np.uint64
    )
    powers = [1]
    for _ in range(1, USER_BLOCKS):
        powers.append(field_multiply_x(powers[-1]))
    coefficients = np.asarray(powers, dtype=np.uint64)
    maximum = np.uint64(FIELD_ORDER)

    def multiply(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        result = np.zeros(left.size, dtype=np.uint64)
        left = left.copy()
        right = right.copy()
        for _ in range(64):
            result ^= np.where(right & np.uint64(1), left, np.uint64(0))
            left = (left << np.uint64(1)) ^ np.where(
                left >> np.uint64(63),
                np.uint64(FIELD_REDUCTION),
                np.uint64(0),
            )
            right >>= np.uint64(1)
        return result

    def support(messages: np.ndarray) -> np.ndarray:
        low = np.zeros(messages.size, dtype=np.uint64)
        high = np.zeros(messages.size, dtype=np.uint64)
        for bit in range(64):
            active = (messages >> np.uint64(bit)) & np.uint64(1)
            mask = np.where(active != 0, maximum, np.uint64(0))
            low ^= mask & generator_low[bit]
            high ^= mask & generator_high[bit]
        return (
            np.count_nonzero(low.view(np.uint8).reshape(-1, 8), axis=1)
            + np.count_nonzero(high.view(np.uint8).reshape(-1, 8), axis=1)
        )

    rng = np.random.default_rng(
        {"one": 0x51, "pair_p0": 0x52, "pair_p1": 0x53, "triple": 0x54}[category]
    )
    probabilities = {
        packets: float(late_probability(packets)) for packets in range(21, 49)
    }
    histogram: Counter[int] = Counter()
    probability_sum = 0.0
    completed = 0
    while completed < sample_count:
        count = min(batch_size, sample_count - completed)
        parameter = rng.integers(1, 1 << 64, size=count, dtype=np.uint64)
        index_count = {"one": 1, "pair_p0": 2, "pair_p1": 2, "triple": 3}[category]
        indices = rng.integers(0, USER_BLOCKS, size=(count, index_count), dtype=np.int64)
        if index_count > 1:
            collision = np.any(np.diff(np.sort(indices, axis=1), axis=1) == 0, axis=1)
            while np.any(collision):
                indices[collision] = rng.integers(
                    0, USER_BLOCKS, size=(int(np.sum(collision)), index_count)
                )
                collision = np.any(
                    np.diff(np.sort(indices, axis=1), axis=1) == 0, axis=1
                )
            indices.sort(axis=1)
        alpha = [coefficients[indices[:, position]] for position in range(index_count)]

        if category == "one":
            messages = (parameter, parameter, multiply(alpha[0], parameter))
        elif category == "pair_p0":
            messages = (
                multiply(alpha[1], parameter),
                multiply(alpha[0], parameter),
                multiply(alpha[0] ^ alpha[1], parameter),
            )
        elif category == "pair_p1":
            messages = (
                parameter,
                parameter,
                multiply(alpha[0] ^ alpha[1], parameter),
            )
        elif category == "triple":
            messages = (
                multiply(alpha[1] ^ alpha[2], parameter),
                multiply(alpha[2] ^ alpha[0], parameter),
                multiply(alpha[0] ^ alpha[1], parameter),
            )
        else:
            raise ValueError(f"unknown support-three category {category}")

        totals = support(messages[0]) + support(messages[1]) + support(messages[2])
        values, counts = np.unique(totals, return_counts=True)
        for packets, frequency in zip(values.tolist(), counts.tolist()):
            histogram[int(packets)] += int(frequency)
            probability_sum += frequency * probabilities[int(packets)]
        completed += count
    return histogram, probability_sum / sample_count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        type=int,
        default=0,
        help="sample this many consecutive-triple parameters after the exact tail",
    )
    parser.add_argument(
        "--gap",
        action="append",
        default=[],
        help="exactly scan one three-data gap pair a,b through 26 packets",
    )
    parser.add_argument(
        "--support3-samples",
        type=int,
        default=0,
        help="sample this many codewords in each outer support-three category",
    )
    args = parser.parse_args()
    if args.samples < 0 or args.support3_samples < 0:
        raise SystemExit("sample counts must be nonnegative")

    packet_enumerator = base.packet_support_enumerator()
    graph_one, graph_pair = graph_families(packet_enumerator)
    graph_power2_pair = graph_power2_unpunctured_pair(packet_enumerator)
    double_lower, double_upper = one_block_bounds(packet_enumerator)
    randomized_exact = randomized_one_block_exact(packet_enumerator)

    low_messages = exact_low_support_messages(MAX_EXACT_LOCAL_PACKET_SUPPORT)
    expected_low_count = sum(
        packet_enumerator[1 : MAX_EXACT_LOCAL_PACKET_SUPPORT + 1]
    )
    if len(low_messages) != expected_low_count:
        raise SystemExit(
            "low-support BCH enumeration has wrong mass: "
            f"{len(low_messages)} != {expected_low_count}"
        )
    triple_enumerator, triple_contribution = consecutive_triple_low_tail(low_messages)
    support_21_count, support_21_contribution = all_triple_support_21(low_messages)

    print("Riffle Double-Parity g=8 late-suffix diagnostic")
    print(f"graph_one_block_log2={base.log2_fraction(graph_one):.12f}")
    print(f"graph_equal_pair_log2={base.log2_fraction(graph_pair):.12f}")
    print(
        "graph_power2_unpunctured_equal_pair_log2="
        f"{base.log2_fraction(graph_power2_pair):.12f}"
    )
    print(f"double_one_block_lower_log2={base.log2_fraction(double_lower):.12f}")
    print(f"double_one_block_upper_log2={base.log2_fraction(double_upper):.12f}")
    print(
        "random_global_multiplier_one_block_exact_log2="
        f"{base.log2_fraction(randomized_exact):.12f}"
    )
    print("double_mds_block_weight,max_packet_terminal_lower_log2")
    for block_weight in range(3, 9):
        print(
            f"{block_weight},"
            f"{base.log2_fraction(mds_terminal_lower_bound(block_weight)):.12f}"
        )
    print(f"exact_bch_messages_support_le_8={len(low_messages)}")
    print("consecutive_triple_total_packets,count")
    for packets, count in sorted(triple_enumerator.items()):
        print(f"{packets},{count}")
    print(
        "consecutive_triple_support_le_26_lower_log2="
        f"{base.log2_fraction(triple_contribution):.12f}"
    )
    print(f"all_triples_total_packet_support_21_count={support_21_count}")
    print(
        "all_triples_total_packet_support_21_log2="
        f"{base.log2_fraction(support_21_contribution):.12f}"
    )
    if args.samples:
        histogram, sampled_mean = sample_consecutive_triple(args.samples)
        exact_low = Fraction()
        for packets, count in triple_enumerator.items():
            exact_low += count * late_probability(packets)
        estimated = (USER_BLOCKS - 2) * (
            float(exact_low) + FIELD_ORDER * sampled_mean
        )
        print(f"consecutive_triple_samples={args.samples}")
        print("sampled_total_packets,count")
        for packets, count in sorted(histogram.items()):
            print(f"{packets},{count}")
        print(
            "consecutive_triple_exact_le_26_plus_sampled_ge_27_log2="
            f"{math.log2(estimated):.12f}"
        )
    for gap_text in args.gap:
        fields = gap_text.split(",")
        if len(fields) != 2:
            raise SystemExit("each --gap must have the form a,b")
        gap1, gap2 = (int(field) for field in fields)
        enumerator, contribution = fixed_gap_triple_low_tail(
            low_messages, gap1, gap2
        )
        print(f"fixed_gap={gap1},{gap2}")
        print(
            "fixed_gap_support_le_26_counts="
            + ",".join(f"{packets}:{count}" for packets, count in sorted(enumerator.items()))
        )
        print(
            "fixed_gap_support_le_26_log2="
            f"{base.log2_fraction(contribution):.12f}"
        )
    if args.support3_samples:
        category_counts = {
            "one": USER_BLOCKS,
            "pair_p0": math.comb(USER_BLOCKS, 2),
            "pair_p1": math.comb(USER_BLOCKS, 2),
            "triple": math.comb(USER_BLOCKS, 3),
        }
        estimated_total = 0.0
        print(f"support3_samples_per_category={args.support3_samples}")
        for category, support_count in category_counts.items():
            histogram, sampled_mean = sample_support_three_category(
                category, args.support3_samples
            )
            estimate = support_count * FIELD_ORDER * sampled_mean
            estimated_total += estimate
            print(f"sampled_{category}_log2={math.log2(estimate):.12f}")
            print(
                f"sampled_{category}_packet_range="
                f"{min(histogram)},{max(histogram)}"
            )
        print(f"sampled_all_support3_log2={math.log2(estimated_total):.12f}")
    print("status=EXACT_FIRST_MOMENT_OBSTRUCTION_DOUBLE_PARITY_G8")


if __name__ == "__main__":
    main()
