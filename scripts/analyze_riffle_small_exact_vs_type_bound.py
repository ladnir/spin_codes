#!/usr/bin/env python3
"""Exact small-instance Riffle spectrum and coefficient-bound comparison.

This is the rational oracle for the scaled extended-BCH [8,4,4]
construction.  It enumerates all GF(16)^B messages, averages the independent
eight-bit block permutations exactly, applies the exact packet-permutation
accumulator enumerator, and compares cumulative expected counts with the
Goal 18 coefficient bound.

The direct message census is intentionally limited to B <= 4.  Larger
instances will use a syndrome dynamic program after this oracle validates the
construction and normalization conventions.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_riffle_packetmul_wrapmul_2lap_bch_family import (  # noqa: E402
    apply_columns,
    build_instance,
)


DEFAULT_OUTPUT_DIRECTORY = (
    ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
)

PACKET_WIDTH = 4
FIELD_SIZE = 16
FIELD_MASK = FIELD_SIZE - 1
FIELD_POLYNOMIAL = 0x13  # x^4 + x + 1
ALPHA_SHIFT = 4


def gf16_multiply(left: int, right: int) -> int:
    result = 0
    a = left
    b = right
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & FIELD_SIZE:
            a ^= FIELD_POLYNOMIAL
    return result & FIELD_MASK


def gf16_power(base: int, exponent: int) -> int:
    result = 1
    value = base
    while exponent:
        if exponent & 1:
            result = gf16_multiply(result, value)
        value = gf16_multiply(value, value)
        exponent >>= 1
    return result


def shifted_coefficients(data_blocks: int) -> tuple[int, ...]:
    coefficients = tuple(
        gf16_power(2, ALPHA_SHIFT + index) for index in range(data_blocks)
    )
    if len(set(coefficients)) != data_blocks or any(value == 0 for value in coefficients):
        raise ValueError("the shifted GF(16) coefficients are not distinct and nonzero")
    return coefficients


def physical_symbols(message: tuple[int, ...], parity_symbols: int) -> tuple[int, ...]:
    if parity_symbols not in (0, 1, 2):
        raise ValueError("parity_symbols must be zero, one, or two")
    result = list(message)
    if parity_symbols:
        p0 = 0
        for value in message:
            p0 ^= value
        result.append(p0)
    if parity_symbols == 2:
        coefficients = shifted_coefficients(len(message))
        p1 = 0
        for coefficient, value in zip(coefficients, message):
            p1 ^= gf16_multiply(coefficient, value)
        result.append(p1)
    return tuple(result)


def symbol_class(value: int) -> int:
    """Return the extended-BCH word weight: 0, 4, or 8."""
    if value == 0:
        return 0
    if value == FIELD_MASK:
        return 8
    return 4


def symbol_class_index(value: int) -> int:
    if value == 0:
        return 0
    if value == FIELD_MASK:
        return 2
    return 1


# Packet-histogram laws after a uniform permutation of one [8,4,4] word.
LOCAL_HISTOGRAM_LAWS: dict[int, tuple[tuple[tuple[int, ...], Fraction], ...]] = {
    0: (((2, 0, 0, 0, 0), Fraction(1)),),
    4: (
        ((1, 0, 0, 0, 1), Fraction(1, 35)),
        ((0, 1, 0, 1, 0), Fraction(16, 35)),
        ((0, 0, 2, 0, 0), Fraction(18, 35)),
    ),
    8: (((0, 0, 0, 0, 2), Fraction(1)),),
}


def add_histograms(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right))


def signature_polynomial(signature: tuple[int, int, int]) -> dict[tuple[int, ...], Fraction]:
    polynomial: dict[tuple[int, ...], Fraction] = {(0, 0, 0, 0, 0): Fraction(1)}
    for weight, multiplicity in zip((0, 4, 8), signature):
        for _ in range(multiplicity):
            following: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
            for current, current_mass in polynomial.items():
                for local, local_mass in LOCAL_HISTOGRAM_LAWS[weight]:
                    following[add_histograms(current, local)] += current_mass * local_mass
            polynomial = dict(following)
    return polynomial


def direct_signature_counts(
    data_blocks: int, parity_symbols: int
) -> dict[tuple[int, int, int], int]:
    signatures: Counter[tuple[int, int, int]] = Counter()
    for message in itertools.product(range(FIELD_SIZE), repeat=data_blocks):
        if not any(message):
            continue
        weights = Counter(symbol_class(value) for value in physical_symbols(message, parity_symbols))
        signatures[(weights[0], weights[4], weights[8])] += 1
    return dict(signatures)


def syndrome_signature_counts(
    data_blocks: int, parity_symbols: int
) -> dict[tuple[int, int, int], int]:
    """Count symbol-class signatures without enumerating 16^B messages."""
    coefficients = shifted_coefficients(data_blocks)
    # State: p0, p1, number of weight-0, weight-4, and weight-8 data blocks.
    current: dict[tuple[int, int, int, int, int], int] = {(0, 0, 0, 0, 0): 1}
    for coefficient in coefficients:
        following: defaultdict[tuple[int, int, int, int, int], int] = defaultdict(int)
        for (p0, p1, zero, ordinary, full), ways in current.items():
            for value in range(FIELD_SIZE):
                counts = [zero, ordinary, full]
                counts[symbol_class_index(value)] += 1
                following[
                    (
                        p0 ^ value,
                        p1 ^ gf16_multiply(coefficient, value),
                        counts[0],
                        counts[1],
                        counts[2],
                    )
                ] += ways
        current = dict(following)

    signatures: Counter[tuple[int, int, int]] = Counter()
    for (p0, p1, zero, ordinary, full), ways in current.items():
        if zero == data_blocks:
            # The unique all-zero message has p0=p1=0.
            continue
        counts = [zero, ordinary, full]
        if parity_symbols:
            counts[symbol_class_index(p0)] += 1
        if parity_symbols == 2:
            counts[symbol_class_index(p1)] += 1
        signatures[tuple(counts)] += ways
    return dict(signatures)


def outer_histogram_distribution(
    data_blocks: int, parity_symbols: int
) -> tuple[dict[tuple[int, ...], Fraction], dict[tuple[int, int, int], int]]:
    signatures = syndrome_signature_counts(data_blocks, parity_symbols)

    result: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
    for signature, messages in signatures.items():
        for histogram, probability in signature_polynomial(signature).items():
            result[histogram] += messages * probability
    return dict(result), dict(signatures)


def transition_multiplicity(old: int, new: int, packet: int) -> int:
    numerator = old + new - packet
    if numerator & 1:
        return 0
    overlap = numerator // 2
    added = new - overlap
    if not 0 <= overlap <= old or not 0 <= added <= PACKET_WIDTH - old:
        return 0
    return math.comb(old, overlap) * math.comb(PACKET_WIDTH - old, added)


TRANSITIONS = tuple(
    tuple(
        tuple(transition_multiplicity(old, new, packet) for new in range(5))
        for old in range(5)
    )
    for packet in range(5)
)


def target_accumulator_enumerators(
    packet_positions: int,
    targets: set[tuple[int, ...]],
    output_weight_maximum: int | None = None,
) -> dict[tuple[int, ...], Counter[int]]:
    """Return A_N(h,w) for targets while pruning impossible partial types."""
    if any(sum(histogram) != packet_positions for histogram in targets):
        raise ValueError("an accumulator target has the wrong packet count")
    allowed = [set() for _ in range(packet_positions + 1)]
    allowed[packet_positions] = set(targets)
    for position in range(packet_positions, 0, -1):
        previous = allowed[position - 1]
        for histogram in allowed[position]:
            for packet, count in enumerate(histogram):
                if count:
                    parent = list(histogram)
                    parent[packet] -= 1
                    previous.add(tuple(parent))

    current: dict[tuple[tuple[int, ...], int, int], int] = {
        ((0, 0, 0, 0, 0), 0, 0): 1
    }
    for position in range(packet_positions):
        following: defaultdict[tuple[tuple[int, ...], int, int], int] = defaultdict(int)
        next_allowed = allowed[position + 1]
        for (counts, old, output_weight), ways in current.items():
            for packet in range(5):
                next_counts = list(counts)
                next_counts[packet] += 1
                next_histogram = tuple(next_counts)
                if next_histogram not in next_allowed:
                    continue
                for new, multiplicity in enumerate(TRANSITIONS[packet][old]):
                    if not multiplicity:
                        continue
                    next_output_weight = output_weight + new
                    if (
                        output_weight_maximum is not None
                        and next_output_weight > output_weight_maximum
                    ):
                        continue
                    following[
                        (next_histogram, new, next_output_weight)
                    ] += ways * multiplicity
        current = dict(following)

    result: defaultdict[tuple[int, ...], Counter[int]] = defaultdict(Counter)
    for (histogram, _state, output_weight), ways in current.items():
        result[histogram][output_weight] += ways
    if not set(result).issubset(targets):
        raise RuntimeError("the targeted accumulator DP produced an untargeted histogram")
    for histogram in targets:
        result[histogram]
    return dict(result)


def histogram_sequence_count(histogram: tuple[int, ...]) -> int:
    packet_positions = sum(histogram)
    result = math.factorial(packet_positions)
    for count in histogram:
        result //= math.factorial(count)
    for packet, count in enumerate(histogram):
        result *= math.comb(PACKET_WIDTH, packet) ** count
    return result


def exact_output_spectrum(
    outer: dict[tuple[int, ...], Fraction],
    accumulator: dict[tuple[int, ...], Counter[int]],
    output_weight_maximum: int | None,
) -> list[Fraction]:
    packet_positions = sum(next(iter(outer)))
    maximum = (
        PACKET_WIDTH * packet_positions
        if output_weight_maximum is None
        else output_weight_maximum
    )
    result = [Fraction(0) for _ in range(maximum + 1)]
    for histogram, outer_count in outer.items():
        denominator = histogram_sequence_count(histogram)
        inner = accumulator[histogram]
        if output_weight_maximum is None and sum(inner.values()) != denominator:
            raise RuntimeError("accumulator histogram mass mismatch")
        for output_weight, inner_count in inner.items():
            result[output_weight] += outer_count * inner_count / denominator
    return result


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        return -math.inf
    return math.log2(value.numerator) - math.log2(value.denominator)


def fraction_payload(value: Fraction) -> dict[str, object]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": format(float(value), ".17g"),
        "log2": log2_fraction(value),
    }


def cumulative(values: list[Fraction]) -> list[Fraction]:
    result = []
    total = Fraction(0)
    for value in values:
        total += value
        result.append(total)
    return result


def first_at_least_one(values: list[Fraction]) -> int | None:
    return next((weight for weight, value in enumerate(values) if value >= 1), None)


def scaled_log_matrix_power_mass(matrix: np.ndarray, exponent: int) -> float:
    vector = np.zeros(5, dtype=np.float64)
    vector[0] = 1.0
    vector_scale = 0.0
    power = matrix.copy()
    power_norm = float(np.max(power))
    power /= power_norm
    power_scale = math.log(power_norm)
    remaining = exponent
    while remaining:
        if remaining & 1:
            vector = vector @ power
            norm = float(np.max(vector))
            vector /= norm
            vector_scale += power_scale + math.log(norm)
        remaining >>= 1
        if remaining:
            power = power @ power
            norm = float(np.max(power))
            power /= norm
            power_scale = 2 * power_scale + math.log(norm)
    return vector_scale + math.log(float(np.sum(vector)))


def log_fraction(value: Fraction) -> float:
    return math.log(value.numerator) - math.log(value.denominator)


def coefficient_bound_objective(
    parameters: np.ndarray,
    distance: int,
    packet_positions: int,
    type_rows: tuple[tuple[tuple[int, ...], float, float], ...],
) -> float:
    log_z = -math.exp(float(parameters[0]))
    log_x = np.empty(5, dtype=np.float64)
    log_x[0] = 0.0
    log_x[1:] = parameters[1:]
    matrix = np.zeros((5, 5), dtype=np.float64)
    for packet in range(5):
        matrix += math.exp(float(log_x[packet])) * np.asarray(TRANSITIONS[packet])
    matrix *= np.exp(log_z * np.arange(5, dtype=np.float64))[None, :]
    log_transfer = scaled_log_matrix_power_mass(matrix, packet_positions)
    terms = [
        log_outer - sum(histogram[k] * log_x[k] for k in range(5)) - log_total
        for histogram, log_outer, log_total in type_rows
    ]
    return -distance * log_z + log_transfer + float(logsumexp(terms))


def optimize_coefficient_bound(
    distance: int,
    packet_positions: int,
    type_rows: tuple[tuple[tuple[int, ...], float, float], ...],
    warm_start: np.ndarray | None,
) -> tuple[float, np.ndarray, bool]:
    bounds = [(-12.0, 3.0)] + [(-14.0, 14.0)] * 4
    starts = [
        np.zeros(5, dtype=np.float64),
        np.array([-2.0, -2.0, -2.0, -2.0, -2.0]),
    ]
    if warm_start is not None:
        starts.insert(0, warm_start)
    results = [
        minimize(
            coefficient_bound_objective,
            start,
            args=(distance, packet_positions, type_rows),
            method="Nelder-Mead",
            bounds=bounds,
            options={"maxiter": 3000, "xatol": 1e-10, "fatol": 1e-10},
        )
        for start in starts
    ]
    best = min(results, key=lambda result: float(result.fun))
    return float(best.fun), best.x, bool(best.success)


def parity_result(
    data_blocks: int,
    parity_symbols: int,
    output_weight_maximum: int | None,
) -> dict[str, object]:
    outer, signatures = outer_histogram_distribution(data_blocks, parity_symbols)
    accumulator = target_accumulator_enumerators(
        2 * (data_blocks + parity_symbols), set(outer), output_weight_maximum
    )
    spectrum = exact_output_spectrum(outer, accumulator, output_weight_maximum)
    tails = cumulative(spectrum)
    crossing = first_at_least_one(tails)
    first_nonzero = next(weight for weight, value in enumerate(tails) if value)
    last_comparison = crossing if crossing is not None else len(tails) - 1

    packet_positions = 2 * (data_blocks + parity_symbols)
    type_rows = tuple(
        (
            histogram,
            log_fraction(count),
            math.log(histogram_sequence_count(histogram)),
        )
        for histogram, count in outer.items()
    )
    comparisons = []
    warm_start = None
    for distance in range(first_nonzero, last_comparison + 1):
        log_bound, warm_start, success = optimize_coefficient_bound(
            distance, packet_positions, type_rows, warm_start
        )
        exact_log2 = log2_fraction(tails[distance])
        bound_log2 = log_bound / math.log(2.0)
        if bound_log2 + 1e-8 < exact_log2:
            raise RuntimeError("numerical coefficient bound fell below the exact tail")
        comparisons.append(
            {
                "distance": distance,
                "exact_cumulative_log2": exact_log2,
                "coefficient_bound_log2": bound_log2,
                "loss_bits": bound_log2 - exact_log2,
                "optimizer_success": success,
            }
        )

    expected_words = (1 << (4 * data_blocks)) - 1
    if output_weight_maximum is None and sum(spectrum, Fraction(0)) != expected_words:
        raise RuntimeError("the exact output spectrum has the wrong total mass")
    if sum(outer.values(), Fraction(0)) != expected_words:
        raise RuntimeError("the outer histogram distribution has the wrong total mass")

    return {
        "parity_symbols": parity_symbols,
        "physical_blocks": data_blocks + parity_symbols,
        "packet_positions": packet_positions,
        "binary_output_length": 4 * packet_positions,
        "exact_output_weight_maximum": (
            4 * packet_positions
            if output_weight_maximum is None
            else output_weight_maximum
        ),
        "complete_output_spectrum": output_weight_maximum is None,
        "message_signatures": len(signatures),
        "outer_packet_histograms": len(outer),
        "exact_expected_count_crossing": crossing,
        "exact_minimum_nonzero_output_weight": first_nonzero,
        "exact_spectrum": [fraction_payload(value) for value in spectrum],
        "exact_cumulative": [fraction_payload(value) for value in tails],
        "coefficient_bound_comparison": comparisons,
    }


def audit_bch8() -> dict[str, object]:
    instance = build_instance(3, 3)
    if instance["columns"] != (0xB, 0xE, 0x7, 0xD):
        raise RuntimeError("the canonical [8,4,4] parity map changed")
    spectrum = Counter()
    for message in range(16):
        parity = apply_columns(instance["columns"], message)
        codeword = message | (parity << 4)
        spectrum[codeword.bit_count()] += 1
        if codeword.bit_count() != symbol_class(message):
            raise RuntimeError("the symbol-class shortcut disagrees with the BCH encoder")
    if spectrum != Counter({4: 14, 0: 1, 8: 1}):
        raise RuntimeError("the [8,4,4] spectrum changed")
    coefficients = shifted_coefficients(4)
    return {
        "systematic_parity_columns_hex": [hex(value) for value in instance["columns"]],
        "ordinary_spectrum": {str(weight): count for weight, count in sorted(spectrum.items())},
        "gf16_polynomial_hex": hex(FIELD_POLYNOMIAL),
        "shifted_coefficients_hex": [hex(value) for value in coefficients],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument(
        "--h-max",
        type=int,
        help="Compute the exact low tail only through this output weight.",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 1 <= args.data_blocks <= 15:
        raise ValueError("the GF(16) coefficient schedule requires 1 <= data_blocks <= 15")
    if args.h_max is not None and args.h_max < 0:
        raise ValueError("h-max must be nonnegative")

    bch_audit = audit_bch8()
    direct_signature_audit = args.data_blocks <= 4
    if direct_signature_audit:
        for parity in range(3):
            direct = direct_signature_counts(args.data_blocks, parity)
            dynamic = syndrome_signature_counts(args.data_blocks, parity)
            if direct != dynamic:
                raise RuntimeError(f"syndrome signature DP failed for parity={parity}")
    results = [
        parity_result(args.data_blocks, parity, args.h_max)
        for parity in range(3)
    ]

    payload = {
        "schema": "riffle-small-exact-vs-type-bound-v1",
        "evidence_label": (
            "EXACT_RATIONAL_FULL_SPECTRUM_AND_NUMERICAL_UPPER_BOUND"
            if args.h_max is None
            else "EXACT_RATIONAL_LOW_TAIL_AND_NUMERICAL_UPPER_BOUND"
        ),
        "construction": "Riffle ShiftAlpha4-eBCH8BlockPerm-ParallelAcc g=4",
        "data_blocks": args.data_blocks,
        "message_bits": 4 * args.data_blocks,
        "bch_audit": bch_audit,
        "parity_results": results,
        "validation": {
            "nonzero_message_mass": (1 << (4 * args.data_blocks)) - 1,
            "direct_message_signatures_match_syndrome_dp_at_all_parity_levels": (
                direct_signature_audit
            ),
            "syndrome_signature_dp_used_for_outer_count": True,
            "exact_output_clipping_has_zero_low_tail_error": args.h_max is not None,
            "local_block_permutation_law": {
                "weight_4_denominator": 70,
                "histogram_counts": [2, 32, 36],
            },
            "all_exact_spectrum_masses_match_nonzero_message_count": True,
            "all_numerical_bounds_dominate_exact_tails": True,
            "all_checks": "PASS",
        },
        "scope": (
            "The exact values average the independent local bit permutations "
            "and the global packet permutation over every nonzero GF(16)^B "
            "message. Output clipping, when requested, is exact through h-max "
            "because accumulator weight never decreases. The coefficient "
            "comparison is a floating-point numerical upper bound, not a proof "
            "certificate."
        ),
    }
    output = args.output or (
        DEFAULT_OUTPUT_DIRECTORY
        / f"goal19_small_exact_vs_type_bound_b{args.data_blocks}.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(
        json.dumps(
            {
                "output": str(output),
                "source_sha256": source_hash,
                "parity_summary": [
                    {
                        "parity_symbols": row["parity_symbols"],
                        "n": row["binary_output_length"],
                        "exact_crossing": row["exact_expected_count_crossing"],
                        "histograms": row["outer_packet_histograms"],
                    }
                    for row in results
                ],
                "status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
