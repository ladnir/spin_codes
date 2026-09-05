#!/usr/bin/env python3
"""Sample a fixed-row-weight binary map and optimize its exact XOR circuit.

The sampler is deterministic from the recorded seed.  Rank and circuit
identities are checked with exact integer arithmetic.  The optimizer is the
repository's Paar common-subexpression heuristic; its result is evidence
about implementation cost, not part of a spectrum certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import probe_bch_forward_xor_circuit as paar  # noqa: E402


def gf2_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            prior = pivots.get(pivot)
            if prior is None:
                pivots[pivot] = row
                break
            row ^= prior
    return len(pivots)


def sample_rows(
    input_bits: int, output_bits: int, degree: int, rng: random.Random
) -> list[int]:
    rows = []
    for _ in range(output_bits):
        row = 0
        for coordinate in rng.sample(range(input_bits), degree):
            row |= 1 << coordinate
        rows.append(row)
    return rows


def rows_hash(rows: list[int], input_bits: int) -> str:
    width = (input_bits + 7) // 8
    return hashlib.sha256(
        b"".join(row.to_bytes(width, "little") for row in rows)
    ).hexdigest()


def circuit_hash(circuit: paar.Circuit) -> str:
    encoded = json.dumps(
        {"gates": circuit.gates, "outputs": circuit.output_signals},
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bits", type=int, required=True)
    parser.add_argument("--output-bits", type=int, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--matrix-seed", type=int, default=1)
    parser.add_argument("--require-full-rank", action="store_true")
    parser.add_argument("--maximum-matrix-attempts", type=int, default=10000)
    parser.add_argument("--optimizer-seed", type=int, default=1)
    parser.add_argument("--optimizer-trials", type=int, default=1)
    parser.add_argument("--minimum-overlap", type=int, default=2)
    parser.add_argument("--rank-samples", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not 0 < args.degree <= args.input_bits:
        raise ValueError("degree must lie in [1,input_bits]")
    maximum_rank = min(args.input_bits, args.output_bits)
    rng = random.Random(args.matrix_seed)

    rank_histogram: dict[int, int] = {}
    for _ in range(args.rank_samples):
        rank = gf2_rank(
            sample_rows(args.input_bits, args.output_bits, args.degree, rng)
        )
        rank_histogram[rank] = rank_histogram.get(rank, 0) + 1

    selected_rows: list[int] | None = None
    selected_rank = -1
    selected_attempt = -1
    for attempt in range(1, args.maximum_matrix_attempts + 1):
        rows = sample_rows(args.input_bits, args.output_bits, args.degree, rng)
        rank = gf2_rank(rows)
        if not args.require_full_rank or rank == maximum_rank:
            selected_rows = rows
            selected_rank = rank
            selected_attempt = attempt
            break
    if selected_rows is None:
        raise RuntimeError("no acceptable matrix was sampled")

    paar.DIMENSION = args.input_bits
    paar.LENGTH = args.output_bits
    baseline = paar.baseline_xors(selected_rows)
    best: tuple[tuple[int, int], paar.Circuit] | None = None
    for seed in range(
        args.optimizer_seed, args.optimizer_seed + args.optimizer_trials
    ):
        circuit = paar.synthesize(seed, args.minimum_overlap, selected_rows)
        score = (circuit.xor_count, seed)
        if best is None or score < best[0]:
            best = (score, circuit)
    assert best is not None
    (optimized, best_seed), circuit = best

    payload = {
        "schema": "sampled-fixed-row-weight-xor-circuit-v1",
        "status": "EXACT_MATRIX_RANK_AND_CIRCUIT_IDENTITY",
        "parameters": {
            "input_bits": args.input_bits,
            "output_bits": args.output_bits,
            "row_degree": args.degree,
            "matrix_seed": args.matrix_seed,
            "require_full_rank": args.require_full_rank,
            "rank_samples_before_selection": args.rank_samples,
            "optimizer_seed_start": args.optimizer_seed,
            "optimizer_trials": args.optimizer_trials,
            "minimum_overlap": args.minimum_overlap,
        },
        "rank_sample_histogram": {
            str(rank): count for rank, count in sorted(rank_histogram.items())
        },
        "selected_matrix": {
            "attempt_after_rank_samples": selected_attempt,
            "rank": selected_rank,
            "rows_sha256": rows_hash(selected_rows, args.input_bits),
            "rows_hex": [hex(row) for row in selected_rows],
        },
        "circuit": {
            "baseline_forward_xors": baseline,
            "optimized_forward_xors": optimized,
            "forward_xors_saved": baseline - optimized,
            "forward_fraction_saved": (baseline - optimized) / baseline,
            "transpose_xors": paar.transpose_xor_count(circuit),
            "shared_gates": len(circuit.gates),
            "max_depth_upper": circuit.max_depth,
            "max_live_signals": circuit.max_live_signals,
            "selected_optimizer_seed": best_seed,
            "circuit_sha256": circuit_hash(circuit),
        },
        "scope": [
            "Every sampled row has exactly the recorded fixed degree.",
            "The selected rank is computed by exact GF(2) elimination.",
            "The forward and transposed straight-line circuits are checked exactly against the sampled matrix.",
            "The rank histogram and heuristic XOR count are measurements, not probabilistic certificates.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "rank_histogram": payload["rank_sample_histogram"],
        "selected_attempt": selected_attempt,
        "selected_rank": selected_rank,
        "baseline_forward_xors": baseline,
        "optimized_forward_xors": optimized,
        "forward_fraction_saved": payload["circuit"]["forward_fraction_saved"],
        "transpose_xors": payload["circuit"]["transpose_xors"],
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
