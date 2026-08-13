#!/usr/bin/env python3
"""Exact small-dimension distance probe for the fixed-split group-chain Riffle.

This is adversarial evidence, not a distance proof.  For K+r<=64 the outer is
one graph-extended EBCH data block followed by its identical parity block.
There are four physical 64-coordinate groups, so every one of their 24 chain
orders can be checked.  Gray-code enumeration then computes the exact minimum
distance for each order.
"""

from __future__ import annotations

import argparse
import itertools
import random
from dataclasses import dataclass


B = 64
MASK = (1 << B) - 1
GENERATOR = 0xF4845518B9582A1F
BCH_ROWS = tuple((GENERATOR << row) | (1 << 127) for row in range(B))


def bch_encode(message: int) -> int:
    word = 0
    while message:
        bit = message & -message
        word ^= BCH_ROWS[bit.bit_length() - 1]
        message ^= bit
    return word


@dataclass(frozen=True)
class Graph:
    rows: tuple[int, ...]

    def encode(self, message: int) -> int:
        result = 0
        for index, row in enumerate(self.rows):
            result |= ((message & row).bit_count() & 1) << index
        return result


def sample_graph(k: int, r: int, rng: random.Random) -> Graph:
    return Graph(tuple(rng.getrandbits(k) for _ in range(r)))


def outer_groups(message: int, graph: Graph, k: int) -> tuple[int, ...]:
    local_message = message | (graph.encode(message) << k)
    word = bch_encode(local_message)
    left = word & MASK
    right = word >> B
    return left, right, left, right


def encode(message: int, graph: Graph, k: int, order: tuple[int, ...]) -> int:
    groups = outer_groups(message, graph, k)
    state = 0
    output = 0
    for position, physical in enumerate(order):
        word = bch_encode(groups[physical] ^ state)
        output |= (word & MASK) << (B * position)
        state = word >> B
    return output


def generator_rows(graph: Graph, k: int, order: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(encode(1 << row, graph, k, order) for row in range(k))


def exact_distance(rows: tuple[int, ...]) -> tuple[int, int]:
    minimum = 1 << 60
    multiplicity = 0
    word = 0
    previous_gray = 0
    for index in range(1, 1 << len(rows)):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        word ^= rows[(changed & -changed).bit_length() - 1]
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
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--codim", type=int, default=24)
    parser.add_argument("--graphs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0x47524F555043484E)
    args = parser.parse_args()
    if not 1 <= args.k <= 24:
        raise SystemExit("exact Gray-code probe expects 1 <= K <= 24")
    if not 0 <= args.codim <= B - args.k:
        raise SystemExit("small outer requires K+codim <= 64")

    orders = tuple(itertools.permutations(range(4)))
    rng = random.Random(args.seed)
    all_distances: list[int] = []
    print(
        f"fixed-split group-chain exact probe K={args.k} r={args.codim} "
        f"graphs={args.graphs} orders={len(orders)} N=256"
    )
    for trial in range(args.graphs):
        graph = sample_graph(args.k, args.codim, rng)
        distances: list[int] = []
        worst_orders: list[tuple[tuple[int, ...], int]] = []
        trial_min = 1 << 60
        for order in orders:
            distance, multiplicity = exact_distance(generator_rows(graph, args.k, order))
            distances.append(distance)
            if distance < trial_min:
                trial_min = distance
                worst_orders = [(order, multiplicity)]
            elif distance == trial_min:
                worst_orders.append((order, multiplicity))
        all_distances.extend(distances)
        ordered = sorted(distances)
        print(
            f"graph={trial} min={ordered[0]} median={ordered[len(ordered)//2]} "
            f"max={ordered[-1]} mean={sum(ordered)/len(ordered):.6f} "
            f"worst_orders={worst_orders}"
        )

    ordered = sorted(all_distances)
    print(
        f"all min={ordered[0]} median={ordered[len(ordered)//2]} max={ordered[-1]} "
        f"mean={sum(ordered)/len(ordered):.6f}"
    )


if __name__ == "__main__":
    main()
