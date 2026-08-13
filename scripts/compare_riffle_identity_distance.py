#!/usr/bin/env python3
"""Exact small-dimension distance comparison: Singer versus identity Riffle.

The outer is the same explicit one-parity-block construction used by the
implementation.  For K<=64 it consists of one BCH-encoded data block and one
identical parity block.  Each trial shares the outer permutation and all local
split permutations between the two inner variants; only the Singer maps are
removed.  After constructing the K generator rows, Gray-code enumeration
examines all 2^K codewords exactly.

This is evidence about true distance, not a replacement for the large-K
first-moment certificate.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass


B = 64
LOCAL_N = 128
GENERATOR = 0xF4845518B9582A1F
BCH_ROWS = tuple((GENERATOR << row) | (1 << 127) for row in range(B))


def bch_encode(message: int) -> int:
    result = 0
    while message:
        bit = (message & -message).bit_length() - 1
        result ^= BCH_ROWS[bit]
        message &= message - 1
    return result


def gf_multiply(left: int, right: int) -> int:
    result = 0
    for _ in range(64):
        if right & 1:
            result ^= left
        carry = left >> 63
        left = (left << 1) & ((1 << 64) - 1)
        if carry:
            left ^= 0x1B
        right >>= 1
    return result


def singer_apply(multiplier: int, value: int) -> int:
    result = 0
    for row in range(64):
        result |= ((value & gf_multiply(multiplier, 1 << row)).bit_count() & 1) << row
    return result


def select_bits(value: int, permutation: tuple[int, ...], begin: int, length: int) -> int:
    result = 0
    for output in range(length):
        result |= ((value >> permutation[begin + output]) & 1) << output
    return result


@dataclass(frozen=True)
class Instance:
    outer_permutation: tuple[int, ...]
    multipliers: tuple[int, ...]
    splits: tuple[tuple[int, ...], ...]


def sample_instance(rng: random.Random, n: int) -> Instance:
    outer = list(range(n))
    rng.shuffle(outer)
    multipliers: list[int] = []
    splits: list[tuple[int, ...]] = []
    for _ in range(n // B):
        multipliers.append(rng.randrange(1, 1 << 64))
        split = list(range(LOCAL_N))
        rng.shuffle(split)
        splits.append(tuple(split))
    return Instance(tuple(outer), tuple(multipliers), tuple(splits))


def encode(message: int, instance: Instance, *, singer: bool) -> int:
    # K<=64: the parity block repeats the only data block.
    local = bch_encode(message)
    n = len(instance.outer_permutation)
    inner_input = 0
    for outer_block in range(2):
        base = outer_block * LOCAL_N
        for coordinate in range(LOCAL_N):
            if (local >> coordinate) & 1:
                inner_input |= 1 << instance.outer_permutation[base + coordinate]

    output = 0
    state = 0
    for block_index, split in enumerate(instance.splits):
        drive = ((inner_input >> (B * block_index)) & ((1 << B) - 1)) ^ state
        if singer and drive:
            drive = singer_apply(instance.multipliers[block_index], drive)
        codeword = bch_encode(drive)
        emitted = select_bits(codeword, split, 0, B)
        state = select_bits(codeword, split, B, B)
        output |= emitted << (B * block_index)
    return output


def generator_rows(k: int, instance: Instance, *, singer: bool) -> tuple[int, ...]:
    return tuple(encode(1 << row, instance, singer=singer) for row in range(k))


def exact_distance(rows: tuple[int, ...]) -> tuple[int, int]:
    minimum = sum(1 for _ in range(1)) * 10**9
    multiplicity = 0
    word = 0
    previous_gray = 0
    for index in range(1, 1 << len(rows)):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        row = (changed & -changed).bit_length() - 1
        word ^= rows[row]
        weight = word.bit_count()
        if weight < minimum:
            minimum = weight
            multiplicity = 1
        elif weight == minimum:
            multiplicity += 1
        previous_gray = gray
    return minimum, multiplicity


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=20)
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0x524946464C45)
    args = parser.parse_args()
    if not 1 <= args.k <= 24:
        raise SystemExit("small-distance comparison expects 1 <= K <= 24")

    n = 2 * LOCAL_N
    rng = random.Random(args.seed)
    singer_distances: list[int] = []
    identity_distances: list[int] = []
    print(f"K={args.k} N={n} trials={args.trials}")
    for trial in range(args.trials):
        instance = sample_instance(rng, n)
        singer_distance, singer_count = exact_distance(
            generator_rows(args.k, instance, singer=True)
        )
        identity_distance, identity_count = exact_distance(
            generator_rows(args.k, instance, singer=False)
        )
        singer_distances.append(singer_distance)
        identity_distances.append(identity_distance)
        print(
            f"trial={trial} singer_d={singer_distance} singer_A_d={singer_count} "
            f"identity_d={identity_distance} identity_A_d={identity_count} "
            f"delta={identity_distance - singer_distance}"
        )

    def summary(values: list[int]) -> str:
        ordered = sorted(values)
        return (
            f"min={ordered[0]} median={ordered[len(ordered)//2]} max={ordered[-1]} "
            f"mean={sum(ordered)/len(ordered):.6f}"
        )

    print(f"singer_summary {summary(singer_distances)}")
    print(f"identity_summary {summary(identity_distances)}")
    wins = sum(left > right for left, right in zip(identity_distances, singer_distances))
    ties = sum(left == right for left, right in zip(identity_distances, singer_distances))
    losses = args.trials - wins - ties
    print(f"identity_vs_singer wins={wins} ties={ties} losses={losses}")


if __name__ == "__main__":
    main()
