#!/usr/bin/env python3
"""Verify the exact pair-type kernel of a permuted prefix accumulator.

For an ordered pair of binary words, record the counts of coordinate symbols
00, 01, 10, and 11.  A common uniform coordinate permutation makes every
symbol sequence with those counts equiprobable.  Prefix accumulation turns
the running XOR state into the output-pair symbol.  A finite dynamic program
counts the resulting output types exactly.

The verifier compares that dynamic program with all labelled coordinate
permutations at a small length.  The identity is length independent; the
small instance only makes exhaustive verification practical.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import itertools
import json
import math
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "accumulator_pair_type_kernel_small_exact.json"


Type4 = tuple[int, int, int, int]


def word_from_symbols(symbols: tuple[int, ...], component: int) -> int:
    result = 0
    for index, symbol in enumerate(symbols):
        shift = 1 - component
        result |= ((symbol >> shift) & 1) << index
    return result


def permute_word(word: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for output, source in enumerate(permutation):
        result |= ((word >> source) & 1) << output
    return result


def accumulate_word(word: int, length: int) -> int:
    result = word
    shift = 1
    while shift < length:
        result ^= result << shift
        shift <<= 1
    return result & ((1 << length) - 1)


def pair_type(first: int, second: int, length: int) -> Type4:
    counts = [0, 0, 0, 0]
    for index in range(length):
        symbol = (((first >> index) & 1) << 1) | ((second >> index) & 1)
        counts[symbol] += 1
    return tuple(counts)  # type: ignore[return-value]


def kernel_counts(input_type: Type4) -> Counter[Type4]:
    """Count output types over distinct input-symbol sequences exactly."""

    # State keys contain used input counts, current XOR state, and output
    # counts.  Choosing a symbol once counts each distinct symbol sequence
    # once.  Uniform labelled permutations induce the uniform law on these
    # sequences because every sequence has prod_a n_a! labelled preimages.
    zero: Type4 = (0, 0, 0, 0)
    table: dict[tuple[Type4, int, Type4], int] = {(zero, 0, zero): 1}
    length = sum(input_type)
    for _ in range(length):
        following: defaultdict[tuple[Type4, int, Type4], int] = defaultdict(int)
        for (used, state, outputs), multiplicity in table.items():
            for symbol in range(4):
                if used[symbol] == input_type[symbol]:
                    continue
                next_used_list = list(used)
                next_used_list[symbol] += 1
                next_used = tuple(next_used_list)
                next_state = state ^ symbol
                next_outputs_list = list(outputs)
                next_outputs_list[next_state] += 1
                next_outputs = tuple(next_outputs_list)
                following[(next_used, next_state, next_outputs)] += multiplicity
        table = dict(following)

    result: Counter[Type4] = Counter()
    for (used, _state, outputs), multiplicity in table.items():
        if used != input_type:
            raise AssertionError("dynamic program ended at the wrong input type")
        result[outputs] += multiplicity
    return result


def brute_labelled_counts(input_type: Type4) -> Counter[Type4]:
    symbols = tuple(
        symbol for symbol, multiplicity in enumerate(input_type) for _ in range(multiplicity)
    )
    length = len(symbols)
    first = word_from_symbols(symbols, 0)
    second = word_from_symbols(symbols, 1)
    result: Counter[Type4] = Counter()
    for permutation in itertools.permutations(range(length)):
        first_out = accumulate_word(permute_word(first, permutation), length)
        second_out = accumulate_word(permute_word(second, permutation), length)
        result[pair_type(first_out, second_out, length)] += 1
    duplicate_factor = math.prod(math.factorial(value) for value in input_type)
    if any(value % duplicate_factor for value in result.values()):
        raise AssertionError("labelled permutation counts have a nonintegral quotient")
    return Counter({key: value // duplicate_factor for key, value in result.items()})


def multinomial(counts: Type4) -> int:
    result = math.factorial(sum(counts))
    for value in counts:
        result //= math.factorial(value)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.length != 8:
        parser.error("the checked-in representative types are defined for length 8")

    representatives: tuple[Type4, ...] = (
        (2, 2, 2, 2),
        (4, 2, 1, 1),
        (3, 3, 2, 0),
        (1, 1, 3, 3),
    )
    rows = []
    for input_type in representatives:
        dynamic = kernel_counts(input_type)
        brute = brute_labelled_counts(input_type)
        if dynamic != brute:
            raise AssertionError(f"pair kernel mismatch for input type {input_type}")
        mass = sum(dynamic.values())
        expected_mass = multinomial(input_type)
        if mass != expected_mass:
            raise AssertionError("kernel row has the wrong mass")
        rows.append(
            {
                "input_type_00_01_10_11": list(input_type),
                "distinct_input_sequences": mass,
                "nonzero_output_types": len(dynamic),
                "output_counts": [
                    {
                        "output_type_00_01_10_11": list(output_type),
                        "count": count,
                    }
                    for output_type, count in sorted(dynamic.items())
                ],
            }
        )
        print(
            f"input_type,{input_type},mass,{mass},output_types,{len(dynamic)}",
            flush=True,
        )

    result = {
        "schema": "accumulator-pair-type-kernel-small-exact-v1",
        "status": "EXACT_INTEGER_EXHAUSTIVE_VERIFICATION",
        "length": args.length,
        "symbol_order": ["00", "01", "10", "11"],
        "kernel": {
            "input": "one common uniform coordinate permutation of an ordered word pair",
            "state_update": "s_i = s_{i-1} xor a_i with s_{-1}=00",
            "output_symbol": "s_i",
            "row_denominator": "multinomial(length; input_type)",
        },
        "rows": rows,
        "claim": "dynamic path counts equal exhaustive labelled-permutation counts after quotienting by duplicate labels",
        "limitations": [
            "The exhaustive check is at length 8.",
            "The counting identity is length independent, but this receipt does not bound the length-512 operator norm.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
