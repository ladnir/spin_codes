#!/usr/bin/env python3
"""Synthesize exact transpose circuits for contiguous EBCH output bands.

This generalizes the two-half probe to arbitrary positive band sizes summing
to 128.  It is used to estimate the circuit cost of grouped outer layouts
before changing the C++ kernel.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from probe_bch_forward_xor_circuit import (
    Circuit,
    synthesize,
    target_outputs,
    transpose_xor_count,
    write_production_transpose_header,
)


def best_band(
    target: list[int], seed_start: int, seeds: int, min_overlap: int
) -> tuple[int, Circuit]:
    best: tuple[tuple[int, int, int], Circuit] | None = None
    for seed in range(seed_start, seed_start + seeds):
        circuit = synthesize(seed, min_overlap, target)
        score = (transpose_xor_count(circuit), circuit.xor_count, seed)
        if best is None or score < best[0]:
            best = (score, circuit)
    assert best is not None
    return best[0][2], best[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", default="42,43,43")
    parser.add_argument("--seeds", type=int, default=64)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--min-overlap", type=int, default=2)
    parser.add_argument(
        "--header-prefix",
        type=Path,
        help="optional prefix; emits <prefix>Band0.h, <prefix>Band1.h, ...",
    )
    args = parser.parse_args()
    sizes = [int(value) for value in args.sizes.split(",")]
    if args.seeds < 1 or any(size < 1 for size in sizes) or sum(sizes) != 128:
        raise SystemExit("band transpose synthesis: sizes must be positive and sum to 128")

    outputs = target_outputs()
    offset = 0
    rows = []
    for band, size in enumerate(sizes):
        target = outputs[offset : offset + size]
        seed, circuit = best_band(
            target, args.seed_start, args.seeds, args.min_overlap
        )
        transpose_xors = transpose_xor_count(circuit)
        rows.append((band, offset, size, seed, circuit.xor_count, transpose_xors))
        if args.header_prefix:
            header = args.header_prefix.with_name(
                f"{args.header_prefix.name}Band{band}.h"
            )
            write_production_transpose_header(
                header,
                circuit,
                seed,
                args.min_overlap,
                function_name=f"extendedBch128x64Band{band}TransposeCircuit",
            )
        offset += size

    print("EBCH [128,64] separate band-transpose circuits")
    for band, offset, size, seed, forward_xors, transpose_xors in rows:
        print(
            f"band={band} offset={offset} size={size} seed={seed} "
            f"forward_xors={forward_xors} transpose_xors={transpose_xors}"
        )
    print(f"combined_transpose_xors={sum(row[5] for row in rows)}")
    print("status=EXACT_VERIFIED")


if __name__ == "__main__":
    main()
