#!/usr/bin/env python3
"""Select a fixed sparse nested A compatible with the recorded B32."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path


DEFAULT_B = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_selection.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_a32_selection.json"
)


def sample_covering_columns(
    rng: random.Random,
    column_weight: int = 7,
    require_coverage: bool = True,
) -> list[int]:
    while True:
        columns = []
        used = set()
        while len(columns) < 32:
            positions = rng.sample(range(64), column_weight)
            value = sum(1 << position for position in positions)
            if value not in used:
                used.add(value)
                columns.append(value)
        coverage = 0
        for column in columns:
            coverage |= column
        if not require_coverage or coverage == (1 << 64) - 1:
            return columns


def prefix_parities(value: int, bits: int = 64) -> int:
    mask = (1 << bits) - 1
    result = value
    shift = 1
    while shift < bits:
        result ^= result << shift
        result &= mask
        shift <<= 1
    return result


def permute_bits(value: int, permutation: list[int]) -> int:
    result = 0
    while value:
        bit = value & -value
        result |= 1 << permutation[bit.bit_length() - 1]
        value ^= bit
    return result


def accumulated_column(
    value: int, permutation: list[int], accumulator_depth: int = 2
) -> int:
    result = value
    for layer in range(accumulator_depth):
        result = prefix_parities(result)
        if layer + 1 < accumulator_depth:
            result = permute_bits(result, permutation)
    return result


def build_a_columns(
    p_columns: list[int],
    permutation: list[int],
    b_mixed: list[int],
    accumulator_depth: int = 2,
) -> list[int]:
    result = []
    for information in range(32):
        auxiliary = accumulated_column(
            p_columns[information], permutation, accumulator_depth
        )
        mixed_word = (1 << information) | (auxiliary << 32)
        syndrome = 0
        selector = mixed_word
        while selector:
            bit = selector & -selector
            syndrome ^= b_mixed[bit.bit_length() - 1]
            selector ^= bit
        result.append(mixed_word | (syndrome << 96))
    return result


def coordinate_forms(a_columns: list[int]) -> list[int]:
    forms = []
    for coordinate in range(128):
        form = 0
        for information, codeword in enumerate(a_columns):
            form |= ((codeword >> coordinate) & 1) << information
        forms.append(form)
    return forms


def coordinate_audit(forms: list[int]) -> tuple[bool, dict[str, int]]:
    zero = sum(form == 0 for form in forms)
    distinct = len(set(forms))
    if zero or distinct != len(forms):
        return (
            False,
            {
                "zero_coordinate_forms": zero,
                "distinct_coordinate_forms": distinct,
                "triple_dependencies": -1,
            },
        )
    form_set = set(forms)
    triple_dependencies = 0
    for left in range(len(forms)):
        for right in range(left + 1, len(forms)):
            target = forms[left] ^ forms[right]
            if target in form_set:
                triple_dependencies += 1
    # Every dependent unordered triple is discovered through its three pairs.
    if triple_dependencies % 3:
        raise ArithmeticError("triple dependency count lost divisibility")
    triple_dependencies //= 3
    return (
        zero == 0 and distinct == len(forms) and triple_dependencies == 0,
        {
            "zero_coordinate_forms": zero,
            "distinct_coordinate_forms": distinct,
            "triple_dependencies": triple_dependencies,
        },
    )


def encode(selector: int, a_columns: list[int]) -> int:
    result = 0
    while selector:
        bit = selector & -selector
        result ^= a_columns[bit.bit_length() - 1]
        selector ^= bit
    return result


def low_input_audit(a_columns: list[int], maximum: int = 4) -> dict[str, object]:
    by_weight = []
    overall = 129
    witness = 0
    for information_weight in range(1, maximum + 1):
        minimum = 129
        count = 0
        for support in itertools.combinations(range(32), information_weight):
            selector = sum(1 << position for position in support)
            output_weight = encode(selector, a_columns).bit_count()
            if output_weight < minimum:
                minimum = output_weight
            if output_weight < overall:
                overall = output_weight
                witness = selector
            count += 1
        by_weight.append(
            {
                "information_weight": information_weight,
                "minimum_output_weight": minimum,
                "support_count": count,
            }
        )
    return {
        "maximum_information_weight": maximum,
        "minimum_output_weight": overall,
        "minimum_witness_hex": f"{witness:08x}",
        "by_information_weight": by_weight,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    b_payload = json.loads(args.b.read_text(encoding="utf-8"))
    b_mixed = [
        int(value, 16)
        for value in b_payload["selected"]["mixed_columns_hex"]
    ]
    best = None
    accepted = 0
    for trial in range(args.trials):
        seed = args.seed + trial
        rng = random.Random(seed)
        p_columns = sample_covering_columns(
            rng, args.p_column_weight, not args.no_require_p_coverage
        )
        permutation = list(range(64))
        rng.shuffle(permutation)
        a_columns = build_a_columns(
            p_columns, permutation, b_mixed, args.accumulator_depth
        )
        forms = coordinate_forms(a_columns)
        valid, form_audit = coordinate_audit(forms)
        if not valid:
            continue
        accepted += 1
        low = low_input_audit(a_columns)
        score = (
            -int(low["minimum_output_weight"]),
            max(codeword.bit_count() for codeword in a_columns),
            seed,
        )
        if best is None or score < best[0]:
            best = (
                score,
                seed,
                p_columns,
                permutation,
                a_columns,
                forms,
                form_audit,
                low,
            )
    if best is None:
        raise RuntimeError("no candidate passed the coordinate audit")
    _, seed, p_columns, permutation, a_columns, forms, form_audit, low = best
    # BA=0 is checked directly using all 128 B columns.
    b_columns = [
        int(value, 16) for value in b_payload["selected"]["columns_hex"]
    ]
    ba_syndromes = []
    for codeword in a_columns:
        syndrome = 0
        selector = codeword
        while selector:
            bit = selector & -selector
            syndrome ^= b_columns[bit.bit_length() - 1]
            selector ^= bit
        ba_syndromes.append(syndrome)
    ba_zero = all(syndrome == 0 for syndrome in ba_syndromes)
    if not ba_zero:
        raise ArithmeticError("constructed A does not satisfy BA=0")
    return {
        "schema": "riffle-splitstate-fixed-a32-selection-v1",
        "candidate": "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
        "parameters": {
            "seed_start": args.seed,
            "trials": args.trials,
            "accepted_coordinate_profiles": accepted,
            "p_column_weight": args.p_column_weight,
            "accumulator_depth": args.accumulator_depth,
        },
        "selected": {
            "seed": seed,
            "p_columns_hex": [f"{value:016x}" for value in p_columns],
            "auxiliary_permutation": permutation,
            "a_columns_hex": [f"{value:032x}" for value in a_columns],
            "coordinate_forms_hex": [f"{value:08x}" for value in forms],
            "coordinate_audit": form_audit,
            "low_input_audit": low,
            "ba_zero": ba_zero,
        },
        "scope": (
            "Deterministic sparse nested search. The coordinate and BA audits "
            "are exact. The low-input audit is not a minimum-distance proof; "
            "the separate full 2^32 enumeration supplies that certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b", type=Path, default=DEFAULT_B)
    parser.add_argument("--seed", type=int, default=0xA3200000)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--p-column-weight", type=int, default=7)
    parser.add_argument("--accumulator-depth", type=int, choices=(1, 2), default=2)
    parser.add_argument("--no-require-p-coverage", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--generators-output", type=Path)
    args = parser.parse_args()
    if args.p_column_weight < 1:
        parser.error("--p-column-weight must be at least 1")
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    selected = payload["selected"]
    if args.generators_output is not None:
        args.generators_output.parent.mkdir(parents=True, exist_ok=True)
        args.generators_output.write_text(
            "\n".join(selected["a_columns_hex"]) + "\n", encoding="utf-8"
        )
    print(f"seed,{selected['seed']}")
    print(
        "accepted_coordinate_profiles,"
        f"{payload['parameters']['accepted_coordinate_profiles']}"
    )
    print(
        "low_input_minimum_output_weight,"
        f"{selected['low_input_audit']['minimum_output_weight']}"
    )
    print(f"ba_zero,{int(selected['ba_zero'])}")
    print(f"wrote,{args.output}")
    if args.generators_output is not None:
        print(f"wrote_generators,{args.generators_output}")


if __name__ == "__main__":
    main()
