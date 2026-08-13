#!/usr/bin/env python3
"""Synthesize exact separate L^T and R^T circuits for EBCH [128,64].

Writing the two 64-coordinate halves separately permits a two-band outer
layout to run each transposed half immediately after its cache-local 64x64
tile transpose.  Candidates are selected by transposed XOR count and verified
exactly by the shared circuit synthesizer.
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


def best_half(
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
    parser.add_argument("--seeds", type=int, default=64)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--min-overlap", type=int, default=2)
    parser.add_argument("--left-header", type=Path)
    parser.add_argument("--right-header", type=Path)
    args = parser.parse_args()
    if args.seeds < 1:
        raise SystemExit("half transpose synthesis: seeds must be positive")
    outputs = target_outputs()
    rows = []
    for name, target, header in (
        ("left", outputs[:64], args.left_header),
        ("right", outputs[64:], args.right_header),
    ):
        seed, circuit = best_half(
            target, args.seed_start, args.seeds, args.min_overlap
        )
        transpose_xors = transpose_xor_count(circuit)
        rows.append((name, seed, circuit.xor_count, transpose_xors))
        if header:
            write_production_transpose_header(
                header,
                circuit,
                seed,
                args.min_overlap,
                function_name=f"extendedBch128x64{name.title()}TransposeCircuit",
            )
    print("EBCH [128,64] separate half-transpose circuits")
    for name, seed, forward_xors, transpose_xors in rows:
        print(
            f"half={name} seed={seed} forward_xors={forward_xors} "
            f"transpose_xors={transpose_xors}"
        )
    print(f"combined_transpose_xors={sum(row[3] for row in rows)}")
    print("status=EXACT_VERIFIED")


if __name__ == "__main__":
    main()
