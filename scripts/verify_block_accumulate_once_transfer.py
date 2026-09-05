#!/usr/bin/env python3
"""Exhaustively verify the BlockAccumulateOnce input-output formula."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter


def choose(n: int, k: int) -> int:
    return math.comb(n, k) if 0 <= k <= n else 0


def formula(length: int, input_weight: int, output_weight: int) -> int:
    return choose(output_weight - 1, (input_weight + 1) // 2 - 1) * choose(
        length - output_weight, input_weight // 2
    )


def accumulate(word: int, length: int) -> int:
    state = 0
    output = 0
    for position in range(length):
        state ^= (word >> position) & 1
        output |= state << position
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=12)
    args = parser.parse_args()
    if not 1 <= args.maximum_length <= 24:
        raise ValueError("maximum length must lie between 1 and 24")

    checked_pairs = 0
    checked_words = 0
    rows = []
    for length in range(1, args.maximum_length + 1):
        observed: Counter[tuple[int, int]] = Counter()
        for word in range(1, 1 << length):
            observed[(word.bit_count(), accumulate(word, length).bit_count())] += 1
            checked_words += 1
        for input_weight in range(1, length + 1):
            total = 0
            for output_weight in range(1, length + 1):
                expected = formula(length, input_weight, output_weight)
                if observed[(input_weight, output_weight)] != expected:
                    raise RuntimeError(
                        "formula mismatch at "
                        f"N={length}, h={input_weight}, w={output_weight}"
                    )
                total += expected
                checked_pairs += 1
            if total != math.comb(length, input_weight):
                raise RuntimeError("input-weight row does not normalize")
        rows.append(
            {
                "length": length,
                "nonzero_inputs": (1 << length) - 1,
                "input_output_cells": length * length,
            }
        )

    payload = {
        "schema": "block-accumulate-once-transfer-audit-v1",
        "maximum_length": args.maximum_length,
        "checked_nonzero_inputs": checked_words,
        "checked_input_output_cells": checked_pairs,
        "length_rows": rows,
        "validation": {
            "exhaustive_input_output_counts_match_formula": True,
            "every_input_weight_row_sums_to_binomial": True,
        },
        "scope": (
            "Exhaustive finite-state audit of the accumulator count T_N(h,w). "
            "The union-bound theorem is an analytic consequence of this count."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
