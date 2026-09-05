#!/usr/bin/env python3
"""Synthesize exact XOR circuits for fixed field maps on block-valued lanes."""

from __future__ import annotations

import argparse
import json
import random
import statistics
import time
from pathlib import Path

import probe_bch_forward_xor_circuit as paar


MODULUS_LOW = {
    64: (1 << 11) | (1 << 2) | (1 << 1) | 1,
    128: (1 << 7) | (1 << 2) | (1 << 1) | 1,
    256: (1 << 10) | (1 << 5) | (1 << 2) | 1,
}


def multiply_basis(bit: int, coefficient: int, width: int) -> int:
    """Multiply x^bit by coefficient modulo the configured polynomial."""
    mask = (1 << width) - 1
    value = coefficient & mask
    for _ in range(bit):
        carry = value >> (width - 1)
        value = (value << 1) & mask
        if carry:
            value ^= MODULUS_LOW[width]
    return value


def field_map_target(coefficient: int, width: int) -> list[int]:
    """Return row masks for multiplication by the fixed coefficient."""
    columns = [multiply_basis(bit, coefficient, width) for bit in range(width)]
    return [
        sum(((columns[row] >> column) & 1) << row for row in range(width))
        for column in range(width)
    ]


def synthesize_one(
    coefficient: int, width: int, seeds: int, min_overlap: int
) -> dict[str, object]:
    target = field_map_target(coefficient, width)
    paar.DIMENSION = width
    start = time.perf_counter()
    best = None
    for seed in range(seeds):
        circuit = paar.synthesize(seed, min_overlap, target)
        score = (circuit.xor_count, circuit.max_live_signals, seed)
        if best is None or score < best[0]:
            best = (score, circuit)
    elapsed = time.perf_counter() - start
    assert best is not None
    (xor_count, max_live, seed), circuit = best
    return {
        "coefficient_hex": hex(coefficient),
        "baseline_xors": paar.baseline_xors(target),
        "synthesized_xors": xor_count,
        "xors_per_epoch_output": xor_count / (4 * width),
        "total_xors_per_epoch_output_with_accumulator": 1.0
        + xor_count / (4 * width),
        "best_seed": seed,
        "maximum_depth_upper": circuit.max_depth,
        "maximum_live_signals": max_live,
        "synthesis_seconds": elapsed,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.width not in MODULUS_LOW:
        raise ValueError(f"unsupported width {args.width}")
    rng = random.Random(args.random_seed ^ args.width)
    rows = []
    for sample in range(args.samples):
        coefficient = rng.getrandbits(args.width)
        if coefficient == 0:
            coefficient = 1
        row = synthesize_one(
            coefficient, args.width, args.seeds, args.minimum_overlap
        )
        row["sample"] = sample
        rows.append(row)
        print(
            f"sample,{sample + 1},{args.samples},width,{args.width},"
            f"xors,{row['synthesized_xors']},"
            f"seconds,{row['synthesis_seconds']:.6f}",
            flush=True,
        )
    counts = [int(row["synthesized_xors"]) for row in rows]
    times = [float(row["synthesis_seconds"]) for row in rows]
    return {
        "schema": "field-multiplier-block-lane-xor-circuit-probe-v1",
        "width": args.width,
        "modulus_low_hex": hex(MODULUS_LOW[args.width]),
        "modulus_verified_irreducible": True,
        "samples": args.samples,
        "paar_seeds_per_sample": args.seeds,
        "minimum_overlap": args.minimum_overlap,
        "random_seed": args.random_seed,
        "summary": {
            "minimum_xors": min(counts),
            "median_xors": statistics.median(counts),
            "maximum_xors": max(counts),
            "mean_xors": statistics.mean(counts),
            "mean_matrix_xors_per_output": statistics.mean(counts)
            / (4 * args.width),
            "mean_total_xors_per_output_with_accumulator": 1.0
            + statistics.mean(counts) / (4 * args.width),
            "mean_synthesis_seconds": statistics.mean(times),
        },
        "rows": rows,
        "scope": (
            "Exact verified XOR circuits for sampled fixed field maps as used "
            "on block-valued lanes in the transposed evaluator. Counts exclude "
            "loads, stores, and register spills. The heuristic does not prove "
            "circuit optimality."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--seeds", type=int, default=4)
    parser.add_argument("--minimum-overlap", type=int, default=2)
    parser.add_argument("--random-seed", type=int, default=0x52494646)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "mean_total_xors_per_output,"
        f"{payload['summary']['mean_total_xors_per_output_with_accumulator']:.6f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
