#!/usr/bin/env python3
"""Measure exact low-tail growth in the GF(16)/eBCH8 Riffle proof gym."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from analyze_riffle_small_exact_vs_type_bound import (
    cumulative,
    exact_output_spectrum,
    first_at_least_one,
    fraction_payload,
    outer_histogram_distribution,
    target_accumulator_enumerators,
)


PARITY_SYMBOLS = 2
DEFAULT_H_MAX = 16
DEFAULT_OUTPUT = Path(
    "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/"
    "receipts/goal27_small_exact_growth_b1_through_b15.json"
)


def source_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def one_size(data_blocks: int, h_max: int) -> dict[str, object]:
    outer, signatures = outer_histogram_distribution(data_blocks, PARITY_SYMBOLS)
    packet_positions = 2 * (data_blocks + PARITY_SYMBOLS)
    accumulator = target_accumulator_enumerators(
        packet_positions,
        set(outer),
        h_max,
    )
    spectrum = exact_output_spectrum(outer, accumulator, h_max)
    tails = cumulative(spectrum)
    crossing = first_at_least_one(tails)
    first_nonzero = next(
        (weight for weight, value in enumerate(tails) if value),
        None,
    )
    expected_words = (1 << (4 * data_blocks)) - 1
    if sum(outer.values(), Fraction(0)) != expected_words:
        raise RuntimeError("outer histogram mass mismatch")
    return {
        "data_blocks": data_blocks,
        "message_bits": 4 * data_blocks,
        "physical_blocks": data_blocks + PARITY_SYMBOLS,
        "packet_positions": packet_positions,
        "binary_output_length": 4 * packet_positions,
        "outer_packet_histograms": len(outer),
        "message_signatures": len(signatures),
        "exact_minimum_nonzero_output_weight": first_nonzero,
        "exact_expected_count_crossing": crossing,
        "crossing_relative_to_output_length": (
            crossing / (4 * packet_positions) if crossing is not None else None
        ),
        "exact_cumulative_through_h_max": [
            fraction_payload(value) for value in tails
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-blocks",
        type=int,
        nargs="+",
        default=list(range(1, 16)),
    )
    parser.add_argument("--h-max", type=int, default=DEFAULT_H_MAX)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if any(not 1 <= blocks <= 15 for blocks in args.data_blocks):
        raise ValueError("the GF(16) schedule requires 1 <= B <= 15")
    if args.h_max < 0:
        raise ValueError("h-max must be nonnegative")

    rows = [one_size(blocks, args.h_max) for blocks in args.data_blocks]
    payload = {
        "schema": "riffle-small-exact-growth-v1",
        "evidence_label": "EXACT_RATIONAL_LOW_TAIL",
        "construction": "Riffle ShiftAlpha4-eBCH8BlockPerm-ParallelAcc g=4",
        "parity_symbols": PARITY_SYMBOLS,
        "h_max": args.h_max,
        "rows": rows,
        "source_sha256": source_sha256(),
        "scope": (
            "Every value averages the independent local bit permutations and "
            "the global packet permutation over all nonzero GF(16)^B messages. "
            "Clipping is exact through h-max because accumulator output weight "
            "never decreases. This finite GF(16) family ends at B=15 and is not "
            "an asymptotic model of the target GF(2^64)/eBCH128 construction."
        ),
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "rows": [
                    {
                        "B": row["data_blocks"],
                        "n": row["binary_output_length"],
                        "crossing": row["exact_expected_count_crossing"],
                        "relative_crossing": row[
                            "crossing_relative_to_output_length"
                        ],
                    }
                    for row in rows
                ],
                "status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
