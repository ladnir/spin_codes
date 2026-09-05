#!/usr/bin/env python3
"""Decompose a uniform permutation into a three-stage local route.

For N=R*C positions, group inputs and outputs into R groups of C positions.
A permutation induces a C-regular bipartite multigraph between input and
output groups.  Repeated perfect matchings color its edges with C colors.
The colors define:

1. R input permutations of size C;
2. C middle permutations of size R;
3. R output permutations of size C.

The script constructs and verifies the target R=256, C=32 route.  It validates
the factorization, not its runtime performance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import struct
from collections import deque
from pathlib import Path


DEFAULT_OUTPUT = Path(
    "constructions/riffle_exactperm_clos256x32_fieldcheckpoint/"
    "receipts/route_validation.json"
)


def perfect_matching(buckets: list[list[list[int]]]) -> tuple[list[int], list[int]]:
    """Return one perfect matching in the support of a regular multigraph."""
    size = len(buckets)
    adjacency = [
        [right for right in range(size) if buckets[left][right]]
        for left in range(size)
    ]
    match_left = [-1] * size
    match_right = [-1] * size
    distance = [0] * size

    def bfs() -> bool:
        queue: deque[int] = deque()
        found = False
        for left in range(size):
            if match_left[left] < 0:
                distance[left] = 0
                queue.append(left)
            else:
                distance[left] = -1
        while queue:
            left = queue.popleft()
            for right in adjacency[left]:
                mate = match_right[right]
                if mate < 0:
                    found = True
                elif distance[mate] < 0:
                    distance[mate] = distance[left] + 1
                    queue.append(mate)
        return found

    def dfs(left: int) -> bool:
        for right in adjacency[left]:
            mate = match_right[right]
            if mate < 0 or (
                distance[mate] == distance[left] + 1 and dfs(mate)
            ):
                match_left[left] = right
                match_right[right] = left
                return True
        distance[left] = -1
        return False

    while bfs():
        for left in range(size):
            if match_left[left] < 0:
                dfs(left)
    if any(right < 0 for right in match_left):
        raise AssertionError("regular bipartite graph lacks a perfect matching")
    return match_left, match_right


def edge_coloring(permutation: list[int], rows: int, columns: int) -> list[int]:
    buckets = [[[] for _ in range(rows)] for _ in range(rows)]
    for source, destination in enumerate(permutation):
        buckets[source // columns][destination // columns].append(source)

    colors = [-1] * len(permutation)
    for color in range(columns):
        match_left, _ = perfect_matching(buckets)
        for left, right in enumerate(match_left):
            edge = buckets[left][right].pop()
            colors[edge] = color
    if any(bucket for row in buckets for bucket in row):
        raise AssertionError("edge coloring did not consume every edge")
    if any(color < 0 for color in colors):
        raise AssertionError("edge coloring left an uncolored edge")
    return colors


def route_schedules(
    permutation: list[int], colors: list[int], rows: int, columns: int
) -> tuple[list[int], list[int], list[int]]:
    size = rows * columns
    first = [-1] * size
    middle = [-1] * size
    last = [-1] * size
    for source, destination in enumerate(permutation):
        input_group = source // columns
        output_group, output_offset = divmod(destination, columns)
        color = colors[source]
        first[input_group * columns + color] = source
        middle[color * rows + output_group] = input_group * columns + color
        last[output_group * columns + output_offset] = color * rows + output_group
    if any(index < 0 for schedule in (first, middle, last) for index in schedule):
        raise AssertionError("route contains an unassigned slot")
    return first, middle, last


def apply_route(
    first: list[int], middle: list[int], last: list[int], rows: int, columns: int
) -> list[int]:
    source = list(range(rows * columns))
    first_values = [source[index] for index in first]
    middle_values = [first_values[index] for index in middle]
    return [middle_values[index] for index in last]


def schedule_hash(schedule: list[int]) -> str:
    digest = hashlib.sha256()
    for value in schedule:
        digest.update(struct.pack("<I", value))
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=256)
    parser.add_argument("--columns", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0xC10525632)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    size = args.rows * args.columns
    permutation = list(range(size))
    random.Random(args.seed).shuffle(permutation)
    colors = edge_coloring(permutation, args.rows, args.columns)
    first, middle, last = route_schedules(
        permutation, colors, args.rows, args.columns
    )
    routed = apply_route(first, middle, last, args.rows, args.columns)
    expected = [-1] * size
    for source, destination in enumerate(permutation):
        expected[destination] = source
    if routed != expected:
        raise AssertionError("three-stage route does not equal the permutation")

    input_color_checks = all(
        len({colors[group * args.columns + offset] for offset in range(args.columns)})
        == args.columns
        for group in range(args.rows)
    )
    output_color_checks = all(
        len(
            {
                colors[source]
                for source, destination in enumerate(permutation)
                if destination // args.columns == group
            }
        )
        == args.columns
        for group in range(args.rows)
    )
    if not input_color_checks or not output_color_checks:
        raise AssertionError("edge colors do not define local permutations")

    payload = {
        "schema": "riffle-exactperm-clos-route-validation-v1",
        "rows": args.rows,
        "columns": args.columns,
        "positions": size,
        "seed": args.seed,
        "correct": True,
        "distribution_statement": (
            "Sampling the input permutation uniformly and decomposing it "
            "preserves the exact uniform-permutation distribution."
        ),
        "stages": [
            {"permutation_count": args.rows, "permutation_size": args.columns},
            {"permutation_count": args.columns, "permutation_size": args.rows},
            {"permutation_count": args.rows, "permutation_size": args.columns},
        ],
        "schedule_sha256": {
            "first": schedule_hash(first),
            "middle": schedule_hash(middle),
            "last": schedule_hash(last),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"correct,{payload['correct']}")
    print(f"positions,{size}")
    print(f"stages,{args.rows}x{args.columns},{args.columns}x{args.rows},{args.rows}x{args.columns}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
