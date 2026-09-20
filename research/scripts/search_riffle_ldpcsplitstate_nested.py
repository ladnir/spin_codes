#!/usr/bin/env python3
"""Search a fixed sparse nested pair for Riffle LDPCSplitState.

The construction uses two nested binary codes C1 <= C0 in F_2^256.
The 64-by-256 matrix B checks C0.  The 192-by-256 stacked matrix
[B; H1] checks C1.  A is a systematic sparse encoder for C1, so BA=0
holds by construction.

This script is a constituent search, not an end-to-end distance proof.
It exhaustively audits encoder inputs of weight at most four and samples
the remaining 64-bit input space.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
from collections import Counter
from pathlib import Path

from explore_riffle_ldpcsplitstate_local import gf2_rank, low_kernel_counts


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_search.json"
)


def sample_unique_masks(
    *, bits: int, weight: int, count: int, rng: random.Random
) -> list[int]:
    result: list[int] = []
    used: set[int] = set()
    while len(result) < count:
        mask = sum(1 << position for position in rng.sample(range(bits), weight))
        if mask not in used:
            used.add(mask)
            result.append(mask)
    return result


def sample_covering_columns(
    *, bits: int, weight: int, count: int, rng: random.Random
) -> list[int]:
    while True:
        columns = sample_unique_masks(
            bits=bits, weight=weight, count=count, rng=rng
        )
        if (sum(columns)).bit_length() >= bits and (
            __import__("functools").reduce(int.__or__, columns, 0)
            == (1 << bits) - 1
        ):
            return columns


def prefix_parities(value: int, bits: int) -> int:
    mask = (1 << bits) - 1
    shift = 1
    result = value
    while shift < bits:
        result ^= result << shift
        result &= mask
        shift <<= 1
    return result


def permute_bits(value: int, permutation: list[int]) -> int:
    result = 0
    while value:
        bit = value & -value
        source = bit.bit_length() - 1
        result |= 1 << permutation[source]
        value ^= bit
    return result


def xor_selected(columns: list[int], selector: int) -> int:
    result = 0
    while selector:
        bit = selector & -selector
        result ^= columns[bit.bit_length() - 1]
        selector ^= bit
    return result


def identity_assignment_avoiding_packet_columns(
    *, s_columns: list[int], r_columns: list[int], rng: random.Random
) -> list[int]:
    """Map packet slots to identity coordinates outside three local supports."""
    candidates = []
    for slot in range(64):
        forbidden = s_columns[slot] | r_columns[slot] | r_columns[64 + slot]
        allowed = [index for index in range(64) if not ((forbidden >> index) & 1)]
        rng.shuffle(allowed)
        candidates.append(allowed)

    coordinate_to_slot = [-1] * 64

    def augment(slot: int, seen: list[bool]) -> bool:
        for coordinate in candidates[slot]:
            if seen[coordinate]:
                continue
            seen[coordinate] = True
            previous = coordinate_to_slot[coordinate]
            if previous < 0 or augment(previous, seen):
                coordinate_to_slot[coordinate] = slot
                return True
        return False

    slot_order = list(range(64))
    rng.shuffle(slot_order)
    for slot in slot_order:
        if not augment(slot, [False] * 64):
            raise RuntimeError("failed to assign parity-lane identity coordinates")
    slot_to_coordinate = [-1] * 64
    for coordinate, slot in enumerate(coordinate_to_slot):
        slot_to_coordinate[slot] = coordinate
    assert sorted(slot_to_coordinate) == list(range(64))
    return slot_to_coordinate


def syndrome_to_parity_coordinates(
    syndrome: int, slot_to_coordinate: list[int]
) -> int:
    parity = 0
    for slot, coordinate in enumerate(slot_to_coordinate):
        parity |= ((syndrome >> coordinate) & 1) << slot
    return parity


def packet_cancellation_score(
    *,
    s_columns: list[int],
    r_columns: list[int],
    identity_permutation: list[int],
) -> tuple[int, int, int]:
    by_syndrome: dict[int, list[tuple[int, int]]] = {}
    symbols: list[list[int]] = []
    zero_count = 0
    for slot in range(64):
        packet_columns = (
            s_columns[slot],
            r_columns[slot],
            r_columns[64 + slot],
            1 << identity_permutation[slot],
        )
        packet_symbols = []
        for pattern in range(1, 16):
            syndrome = 0
            for lane, column in enumerate(packet_columns):
                if (pattern >> lane) & 1:
                    syndrome ^= column
            zero_count += int(syndrome == 0)
            packet_symbols.append(syndrome)
            by_syndrome.setdefault(syndrome, []).append((slot, pattern))
        symbols.append(packet_symbols)

    pair_count = 0
    for occurrences in by_syndrome.values():
        for left in range(len(occurrences)):
            for right in range(left + 1, len(occurrences)):
                pair_count += int(occurrences[left][0] != occurrences[right][0])

    triple_count = 0
    for left_slot in range(64):
        for right_slot in range(left_slot + 1, 64):
            for left_syndrome in symbols[left_slot]:
                for right_syndrome in symbols[right_slot]:
                    target = left_syndrome ^ right_syndrome
                    triple_count += sum(
                        third_slot > right_slot
                        for third_slot, _ in by_syndrome.get(target, ())
                    )
    return zero_count, pair_count, triple_count


def optimize_identity_wiring(
    pair: dict[str, object], *, trials: int, seed: int
) -> tuple[int, int, int]:
    rng = random.Random(seed ^ 0x4944454E54495459)
    best_score: tuple[int, int, int] | None = None
    best_permutation: list[int] | None = None
    for _ in range(trials):
        permutation = identity_assignment_avoiding_packet_columns(
            s_columns=pair["s_columns"], r_columns=pair["r_columns"], rng=rng
        )
        score = packet_cancellation_score(
            s_columns=pair["s_columns"],
            r_columns=pair["r_columns"],
            identity_permutation=permutation,
        )
        if best_score is None or score < best_score:
            best_score = score
            best_permutation = permutation
        if score == (0, 0, 0):
            break
    assert best_score is not None and best_permutation is not None

    pair["identity_permutation"] = best_permutation
    pair["b_columns"] = (
        pair["s_columns"]
        + pair["r_columns"]
        + [1 << coordinate for coordinate in best_permutation]
    )
    rebuilt_a_columns = []
    for index, codeword in enumerate(pair["a_columns"]):
        auxiliary = (codeword >> 64) & ((1 << 128) - 1)
        parity_syndrome = pair["s_columns"][index] ^ xor_selected(
            pair["r_columns"], auxiliary
        )
        parity = syndrome_to_parity_coordinates(
            parity_syndrome, best_permutation
        )
        rebuilt_a_columns.append(
            (codeword & ((1 << 192) - 1)) | (parity << 192)
        )
    pair["a_columns"] = rebuilt_a_columns
    assert verify_ba(pair)
    return best_score


def make_pair(
    *,
    seed: int,
    p_column_weight: int,
    b_column_weight: int,
    accumulator_depth: int,
    fixed_b: dict[str, object] | None = None,
) -> dict[str, object]:
    rng = random.Random(seed)

    # H1 has the form [P T 0], where T is the invertible accumulator
    # difference matrix.  Its inverse computes prefix parities.
    p_columns = sample_covering_columns(
        bits=128, weight=p_column_weight, count=64, rng=rng
    )
    accumulator_permutations: list[list[int]] = []
    for _ in range(accumulator_depth - 1):
        permutation = list(range(128))
        rng.shuffle(permutation)
        accumulator_permutations.append(permutation)

    # B has the form [S R I].  Sampling S and R from one unique pool
    # removes duplicate columns between those two coordinate classes.
    if fixed_b is None:
        sr_columns = sample_unique_masks(
            bits=64, weight=b_column_weight, count=192, rng=rng
        )
        s_columns = sr_columns[:64]
        r_columns = sr_columns[64:]
        identity_permutation = identity_assignment_avoiding_packet_columns(
            s_columns=s_columns, r_columns=r_columns, rng=rng
        )
    else:
        s_columns = [
            int(value, 16) for value in fixed_b["matrices"]["S_columns_hex_64"]
        ]
        r_columns = [
            int(value, 16) for value in fixed_b["matrices"]["R_columns_hex_64"]
        ]
        identity_permutation = list(
            fixed_b["matrices"]["parity_lane_identity_permutation"]
        )
    identity_columns = [1 << index for index in identity_permutation]
    b_columns = s_columns + r_columns + identity_columns

    # A(q) = (q, a, p), with a=T^-1 Pq and p=S q + R a.
    a_columns: list[int] = []
    for index in range(64):
        auxiliary = p_columns[index]
        for layer in range(accumulator_depth):
            auxiliary = prefix_parities(auxiliary, 128)
            if layer < accumulator_depth - 1:
                auxiliary = permute_bits(
                    auxiliary, accumulator_permutations[layer]
                )
        parity_syndrome = s_columns[index] ^ xor_selected(r_columns, auxiliary)
        parity = syndrome_to_parity_coordinates(
            parity_syndrome, identity_permutation
        )
        codeword = (1 << index) | (auxiliary << 64) | (parity << 192)
        a_columns.append(codeword)

    return {
        "seed": seed,
        "p_columns": p_columns,
        "s_columns": s_columns,
        "r_columns": r_columns,
        "b_columns": b_columns,
        "a_columns": a_columns,
        "accumulator_permutations": accumulator_permutations,
        "identity_permutation": identity_permutation,
    }


def support_audit(columns: list[int], maximum_input_weight: int) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    global_minimum = 257
    low_23 = 0
    low_55 = 0
    minimum_witness: list[int] = []
    for input_weight in range(1, maximum_input_weight + 1):
        histogram: Counter[int] = Counter()
        row_minimum = 257
        row_witness: tuple[int, ...] = ()
        for support in itertools.combinations(range(64), input_weight):
            value = 0
            for index in support:
                value ^= columns[index]
            histogram[value.bit_count()] += 1
            if value.bit_count() < row_minimum:
                row_minimum = value.bit_count()
                row_witness = support
        minimum = min(histogram)
        if minimum < global_minimum:
            minimum_witness = list(row_witness)
        global_minimum = min(global_minimum, minimum)
        local_23 = sum(count for weight, count in histogram.items() if weight <= 23)
        local_55 = sum(count for weight, count in histogram.items() if weight <= 55)
        low_23 += local_23
        low_55 += local_55
        rows.append(
            {
                "input_weight": input_weight,
                "support_count": sum(histogram.values()),
                "minimum_output_weight": minimum,
                "minimum_witness": list(row_witness),
                "maximum_output_weight": max(histogram),
                "mean_output_weight": sum(
                    weight * count for weight, count in histogram.items()
                )
                / sum(histogram.values()),
                "output_weight_le_23": local_23,
                "output_weight_le_55": local_55,
                "histogram": {str(k): histogram[k] for k in sorted(histogram)},
            }
        )
    return {
        "maximum_input_weight": maximum_input_weight,
        "minimum_output_weight": global_minimum,
        "minimum_witness": minimum_witness,
        "output_weight_le_23": low_23,
        "output_weight_le_55": low_55,
        "by_input_weight": rows,
    }


def quick_score(columns: list[int]) -> tuple[int, int, int]:
    """Cheap score used before the exact weight-three/four audit."""
    minimum = 257
    below_56 = 0
    total_weight = 0
    count = 0
    for input_weight in (1, 2, 3):
        for support in itertools.combinations(range(64), input_weight):
            value = 0
            for index in support:
                value ^= columns[index]
            weight = value.bit_count()
            minimum = min(minimum, weight)
            below_56 += int(weight <= 55)
            total_weight += weight
            count += 1
    return minimum, -below_56, total_weight // count


def output_coordinate_rows(columns: list[int]) -> list[int]:
    return [
        sum(((columns[input_index] >> output) & 1) << input_index for input_index in range(64))
        for output in range(256)
    ]


def coordinate_dependency_score(columns: list[int]) -> tuple[int, int, int]:
    rows = output_coordinate_rows(columns)
    zero_rows = sum(row == 0 for row in rows)
    duplicate_rows = len(rows) - len(set(rows))
    row_set = set(rows)
    triples: set[tuple[int, int, int]] = set()
    row_index = {row: index for index, row in enumerate(rows)}
    for left in range(256):
        for right in range(left + 1, 256):
            third = row_index.get(rows[left] ^ rows[right])
            if third is not None and third not in (left, right):
                triples.add(tuple(sorted((left, right, third))))
    return zero_rows, duplicate_rows, len(triples)


def verify_ba(pair: dict[str, object]) -> bool:
    b_columns = pair["b_columns"]
    for codeword in pair["a_columns"]:
        if xor_selected(b_columns, codeword) != 0:
            return False
    return True


def sampled_audit(
    columns: list[int], *, samples: int, seed: int
) -> dict[str, object]:
    rng = random.Random(seed ^ 0x53414D504C45)
    histogram: Counter[int] = Counter()
    minimum = 257
    for _ in range(samples):
        selector = rng.getrandbits(64)
        if selector == 0:
            selector = 1
        weight = xor_selected(columns, selector).bit_count()
        histogram[weight] += 1
        minimum = min(minimum, weight)
    return {
        "samples": samples,
        "minimum_output_weight": minimum,
        "mean_output_weight": sum(
            weight * count for weight, count in histogram.items()
        )
        / samples,
        "histogram": {str(k): histogram[k] for k in sorted(histogram)},
    }


def hill_climb_audit(
    columns: list[int], *, restarts: int, seed: int
) -> dict[str, object]:
    rng = random.Random(seed ^ 0x48494C4C434C494D42)
    best_weight = 257
    best_selector = 0
    local_minima: Counter[int] = Counter()
    for _ in range(restarts):
        selector = rng.getrandbits(64)
        if selector == 0:
            selector = 1
        codeword = xor_selected(columns, selector)
        weight = codeword.bit_count()
        while True:
            following_weight = weight
            following_index = -1
            for index, column in enumerate(columns):
                candidate_weight = (codeword ^ column).bit_count()
                if candidate_weight < following_weight:
                    following_weight = candidate_weight
                    following_index = index
            if following_index < 0:
                break
            selector ^= 1 << following_index
            codeword ^= columns[following_index]
            weight = following_weight
        if selector == 0:
            # The zero word is not an admissible distance witness.
            continue
        local_minima[weight] += 1
        if weight < best_weight:
            best_weight = weight
            best_selector = selector
    return {
        "restarts": restarts,
        "minimum_output_weight": best_weight,
        "information_hex": f"{best_selector:016x}",
        "information_weight": best_selector.bit_count(),
        "local_minimum_histogram": {
            str(k): local_minima[k] for k in sorted(local_minima)
        },
        "scope": "heuristic local search; not a distance certificate",
    }


def xor_costs(
    pair: dict[str, object],
    *,
    p_column_weight: int,
    b_column_weight: int,
    accumulator_depth: int,
) -> dict[str, object]:
    p_edges = 64 * p_column_weight
    b_sr_edges = 192 * b_column_weight
    return {
        "Pq": p_edges - 128,
        "auxiliary_accumulators": accumulator_depth * 127,
        "final_parity_Sq_plus_Ra": b_sr_edges - 64,
        "A_total": (
            p_edges - 128
            + accumulator_depth * 127
            + b_sr_edges - 64
        ),
        "A_xor_per_256_output_bits": (
            p_edges - 128
            + accumulator_depth * 127
            + b_sr_edges - 64
        )
        / 256,
        "B_total": b_sr_edges,
        "B_xor_per_256_input_bits": b_sr_edges / 256,
        "apply_A_to_input": 256,
        "apply_A_to_input_xor_per_output": 1.0,
        "notes": (
            "Logical XOR count before circuit sharing, fusion, SIMD packing, "
            "or the GF(2^64) state randomizer. B includes the XOR with the "
            "identity parity coordinate in each of its 64 rows."
        ),
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    fixed_b = (
        json.loads(args.fixed_b_source.read_text(encoding="utf-8"))
        if args.fixed_b_source is not None
        else None
    )
    candidates: list[tuple[tuple[int, int, int], dict[str, object]]] = []
    for candidate_index in range(args.candidates):
        seed = args.seed + candidate_index
        pair = make_pair(
            seed=seed,
            p_column_weight=args.p_column_weight,
            b_column_weight=args.b_column_weight,
            accumulator_depth=args.accumulator_depth,
            fixed_b=fixed_b,
        )
        assert verify_ba(pair)
        dependency_score = coordinate_dependency_score(pair["a_columns"])
        expansion_score = quick_score(pair["a_columns"])
        combined_score = (
            -dependency_score[0],
            -dependency_score[1],
            -dependency_score[2],
            *expansion_score,
        )
        candidates.append((combined_score, pair))
    candidates.sort(key=lambda item: item[0], reverse=True)
    score, best = candidates[0]
    packet_score = (
        packet_cancellation_score(
            s_columns=best["s_columns"],
            r_columns=best["r_columns"],
            identity_permutation=best["identity_permutation"],
        )
        if fixed_b is not None
        else optimize_identity_wiring(
            best, trials=args.identity_wiring_trials, seed=args.seed
        )
    )

    support = support_audit(best["a_columns"], args.maximum_input_weight)
    b_counts = low_kernel_counts(best["b_columns"])
    a_dual_counts = low_kernel_counts(output_coordinate_rows(best["a_columns"]))
    zero_output_coordinates = [
        output
        for output in range(256)
        if not any((column >> output) & 1 for column in best["a_columns"])
    ]
    payload = {
        "schema": "riffle-ldpcsplitstate-fixed-nested-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "construction": {
            "coordinates": "(q[64], auxiliary[128], parity[64])",
            "B": "[S R I_64]",
            "H1": (
                "[P T 0], with T the inverse of the layered accumulator "
                "transform"
            ),
            "A": "q -> (q, L(Pq), Sq + R L(Pq))",
            "L": (
                f"{args.accumulator_depth} accumulator layer(s), with a fixed "
                "permutation between adjacent layers"
            ),
            "nesting": "B A = 0 exactly",
        },
        "parameters": {
            "candidates": args.candidates,
            "seed": args.seed,
            "chosen_seed": best["seed"],
            "P_column_weight": args.p_column_weight,
            "S_and_R_column_weight": args.b_column_weight,
            "accumulator_depth": args.accumulator_depth,
            "maximum_exhaustive_input_weight": args.maximum_input_weight,
            "random_samples": args.random_samples,
            "hill_climb_restarts": args.hill_climb_restarts,
            "identity_wiring_trials": args.identity_wiring_trials,
            "fixed_B_source": (
                str(args.fixed_b_source) if args.fixed_b_source is not None else None
            ),
        },
        "ranks": {
            "B": gf2_rank(best["b_columns"], 64),
            "A": gf2_rank(best["a_columns"], 256),
            "stacked_B_H1": 192,
        },
        "verification": {
            "BA_is_zero": verify_ba(best),
            "quick_search_score": list(score),
            "packet_cancellation_score": list(packet_score),
            "A_full_output_support": not zero_output_coordinates,
            "A_zero_output_coordinates": zero_output_coordinates,
        },
        "A_low_input_support_audit": support,
        "A_random_input_audit": sampled_audit(
            best["a_columns"], samples=args.random_samples, seed=args.seed
        ),
        "A_hill_climb_audit": hill_climb_audit(
            best["a_columns"], restarts=args.hill_climb_restarts, seed=args.seed
        ),
        "B_kernel_counts": b_counts,
        "A_dual_low_weight_counts": a_dual_counts,
        "xor_costs": xor_costs(
            best,
            p_column_weight=args.p_column_weight,
            b_column_weight=args.b_column_weight,
            accumulator_depth=args.accumulator_depth,
        ),
        "matrices": {
            "P_columns_hex_128": [f"{x:032x}" for x in best["p_columns"]],
            "S_columns_hex_64": [f"{x:016x}" for x in best["s_columns"]],
            "R_columns_hex_64": [f"{x:016x}" for x in best["r_columns"]],
            "B_columns_hex_64": [f"{x:016x}" for x in best["b_columns"]],
            "A_columns_hex_256": [f"{x:064x}" for x in best["a_columns"]],
            "accumulator_permutations": best["accumulator_permutations"],
            "parity_lane_identity_permutation": best["identity_permutation"],
        },
        "scope": (
            "Fixed constituent audit only. The exhaustive statement covers "
            "state inputs of Hamming weight at most the recorded maximum. "
            "The random-input histogram is diagnostic, not a proof. No affine "
            "input translate, outer spectrum, placement, or end-to-end union "
            "bound is included."
        ),
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0x4E455354)
    parser.add_argument("--p-column-weight", type=int, default=6)
    parser.add_argument("--b-column-weight", type=int, default=3)
    parser.add_argument("--accumulator-depth", type=int, default=1)
    parser.add_argument("--maximum-input-weight", type=int, default=4)
    parser.add_argument("--random-samples", type=int, default=200_000)
    parser.add_argument("--hill-climb-restarts", type=int, default=20_000)
    parser.add_argument("--identity-wiring-trials", type=int, default=64)
    parser.add_argument("--fixed-b-source", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    support = payload["A_low_input_support_audit"]
    sampled = payload["A_random_input_audit"]
    hill = payload["A_hill_climb_audit"]
    print(f"chosen_seed,{payload['parameters']['chosen_seed']}")
    print(f"BA_is_zero,{payload['verification']['BA_is_zero']}")
    print(f"A_low_support_minimum,{support['minimum_output_weight']}")
    print(f"A_random_sample_minimum,{sampled['minimum_output_weight']}")
    print(f"A_random_sample_mean,{sampled['mean_output_weight']:.6f}")
    print(
        f"A_hill_climb_minimum,{hill['minimum_output_weight']},"
        f"input_weight,{hill['information_weight']}"
    )
    print(
        "B_kernel_counts,"
        f"w2,{payload['B_kernel_counts']['weight_2']},"
        f"w3,{payload['B_kernel_counts']['weight_3']},"
        f"w4,{payload['B_kernel_counts']['weight_4']}"
    )
    print(
        "logical_xor_per_output,"
        f"A,{payload['xor_costs']['A_xor_per_256_output_bits']:.6f},"
        f"B,{payload['xor_costs']['B_xor_per_256_input_bits']:.6f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
