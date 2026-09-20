#!/usr/bin/env python3
"""Proof and refutation receipts for Riffle DP g=4 weighted triples.

Proof mode constructs an exact one-Krawtchouk Delsarte certificate for an
additive length-96 code over the 16-element nibble alphabet. Refutation mode
authenticates the exact support-33 terminal family from a complete local
support file. Neither mode silently treats a diagnostic as a proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
MANIFEST = REPOSITORY / "explorations" / "riffle_dp_g4_construction_manifest.json"
OUTPUT_DIRECTORY = REPOSITORY / "explorations"
Q = 16
CODE_SIZE = 1 << 64
MASK64 = CODE_SIZE - 1
FIELD_REDUCTION = 0x1B
LOCAL_GENERATOR = 0xF4845518B9582A1F
TRIPLE_LENGTH = 96
DATA_SYMBOLS = 1 << 14
PACKET_COUNT = 524_352
TERMINAL_PACKETS = 47_184


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def terminal_probability(support: int) -> Fraction:
    return Fraction(
        math.comb(TERMINAL_PACKETS, support),
        math.comb(PACKET_COUNT, support),
    )


def field_multiply(left: int, right: int) -> int:
    result = 0
    for _ in range(64):
        if right & 1:
            result ^= left
        carry = left >> 63
        left = (left << 1) & MASK64
        if carry:
            left ^= FIELD_REDUCTION
        right >>= 1
    return result


def field_power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = field_multiply(result, value)
        value = field_multiply(value, value)
        exponent >>= 1
    return result


def encode_local(message: int) -> int:
    codeword = 0
    for bit in range(64):
        if (message >> bit) & 1:
            low = (LOCAL_GENERATOR << bit) & MASK64
            high = (LOCAL_GENERATOR >> (64 - bit) if bit else 0) | (1 << 63)
            codeword ^= low | (high << 64)
    return codeword


def nibble_support(codeword: int) -> int:
    return sum(bool((codeword >> (4 * packet)) & 0xF) for packet in range(32))


def read_low_messages(path: Path) -> dict[int, set[int]]:
    data = path.read_bytes()
    if len(data) % 9:
        raise ValueError("local-message file has a partial record")
    result: dict[int, set[int]] = {}
    for offset in range(0, len(data), 9):
        support, message = struct.unpack_from("<BQ", data, offset)
        codeword = encode_local(message)
        if nibble_support(codeword) != support:
            raise ValueError("local-message record failed independent re-encoding")
        result.setdefault(support, set()).add(message)
    return result


def krawtchouk(degree: int, weight: int) -> int:
    lower = max(0, degree - (TRIPLE_LENGTH - weight))
    upper = min(degree, weight)
    return sum(
        (-1) ** h
        * (Q - 1) ** (degree - h)
        * math.comb(weight, h)
        * math.comb(TRIPLE_LENGTH - weight, degree - h)
        for h in range(lower, upper + 1)
    )


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    if value <= 0:
        raise ValueError("logarithm requires a positive fraction")
    quantum = Decimal(1).scaleb(-places)
    with localcontext() as context:
        context.prec = 100
        estimate = (Decimal(value.numerator) / Decimal(value.denominator)).ln()
        estimate /= Decimal(2).ln()
        lower = estimate.quantize(quantum, rounding=ROUND_FLOOR)
        upper = estimate.quantize(quantum, rounding=ROUND_CEILING)
    return [format(lower, "f"), format(upper, "f")]


def exact_single_krawtchouk_certificate(
    minimum_support: int,
    degree: int,
) -> dict:
    """Return the optimal certificate using one dual Krawtchouk constraint."""
    if not 1 <= degree <= TRIPLE_LENGTH:
        raise ValueError("Krawtchouk degree is out of range")
    if not 1 <= minimum_support <= TRIPLE_LENGTH:
        raise ValueError("minimum support is out of range")

    base_probability = terminal_probability(minimum_support)
    k_zero = krawtchouk(degree, 0)
    lines = []
    for weight in range(minimum_support, TRIPLE_LENGTH + 1):
        normalized_probability = terminal_probability(weight) / base_probability
        normalized_krawtchouk = Fraction(krawtchouk(degree, weight), k_zero)
        lines.append((weight, normalized_probability, normalized_krawtchouk))

    candidates = {Fraction()}
    for left in range(len(lines)):
        _, left_c, left_r = lines[left]
        for right in range(left + 1, len(lines)):
            _, right_c, right_r = lines[right]
            denominator = left_r - right_r
            if denominator:
                candidate = (right_c - left_c) / denominator
                if candidate >= 0:
                    candidates.add(candidate)

    best = None
    for dual_weight in candidates:
        intercepts = [
            normalized_probability + dual_weight * normalized_krawtchouk
            for _, normalized_probability, normalized_krawtchouk in lines
        ]
        z_value = max(intercepts)
        objective = (CODE_SIZE - 1) * z_value + dual_weight
        if best is None or objective < best[0]:
            best = (objective, dual_weight, z_value, intercepts)
    assert best is not None
    objective, dual_weight, z_value, intercepts = best

    for (_, normalized_probability, normalized_krawtchouk), intercept in zip(
        lines, intercepts
    ):
        if normalized_probability > z_value - dual_weight * normalized_krawtchouk:
            raise RuntimeError("exact dual inequality failed")
        if intercept > z_value:
            raise RuntimeError("exact envelope verification failed")

    active_weights = [
        weight
        for (weight, _, _), intercept in zip(lines, intercepts)
        if intercept == z_value
    ]
    bound = base_probability * objective
    class_count = math.comb(DATA_SYMBOLS + 2, 3)
    aggregate_if_all_classes_qualify = class_count * bound
    return {
        "schema": "riffle-dp-g4-weighted-triple-proof-receipt-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "manifest_sha256": sha256(MANIFEST),
        "evidence_label": "RIGOROUS_BOUND",
        "scope": (
            "One coefficient-triple additive code whose every nonzero word "
            f"has nibble support at least {minimum_support}."
        ),
        "alphabet_size": Q,
        "length": TRIPLE_LENGTH,
        "code_size": str(CODE_SIZE),
        "minimum_support_assumption": minimum_support,
        "krawtchouk_degree": degree,
        "dual_weight": {
            "numerator": str(dual_weight.numerator),
            "denominator": str(dual_weight.denominator),
        },
        "dual_intercept": {
            "numerator": str(z_value.numerator),
            "denominator": str(z_value.denominator),
        },
        "active_weights": active_weights,
        "per_class_terminal_bound": {
            "numerator": str(bound.numerator),
            "denominator": str(bound.denominator),
            "log2_interval": log2_interval(bound),
        },
        "weight_three_support_class_count": str(class_count),
        "hypothetical_all_classes_bound": {
            "condition": (
                f"Valid only if every coefficient class has minimum support {minimum_support}."
            ),
            "numerator": str(aggregate_if_all_classes_qualify.numerator),
            "denominator": str(aggregate_if_all_classes_qualify.denominator),
            "log2_interval": log2_interval(aggregate_if_all_classes_qualify),
            "meets_2^-40": aggregate_if_all_classes_qualify <= Fraction(1, 1 << 40),
        },
        "derivation": (
            "MacWilliams gives sum_i A_i K_j(i)/K_j(0) >= -1. "
            "The verified envelope P(i)/P(d) <= z-u K_j(i)/K_j(0) "
            "therefore gives sum_i A_i P(i) <= P(d)((2^64-1)z+u)."
        ),
    }


def write_proof_receipt(minimum_support: int, degree: int, output: Path) -> None:
    receipt = exact_single_krawtchouk_certificate(minimum_support, degree)
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"receipt={output}")
    print(f"per_class_log2={receipt['per_class_terminal_bound']['log2_interval']}")
    print(
        "hypothetical_all_classes_log2="
        f"{receipt['hypothetical_all_classes_bound']['log2_interval']}"
    )
    print("status=RIGOROUS_BOUND_CONDITIONAL_ON_DECLARED_MINIMUM_SUPPORT")


def write_refutation_receipt(low_message_file: Path, output: Path) -> None:
    messages = read_low_messages(low_message_file)
    if any(messages.get(support) for support in range(1, 11)):
        raise ValueError("local file contains a word below support 11")
    support_11 = messages.get(11, set())
    if len(support_11) != 20:
        raise ValueError("support-11 layer is not the complete 20-word layer")

    powers = []
    value = 1
    for _ in range(DATA_SYMBOLS):
        powers.append(value)
        value = field_multiply(value, 2)
    power_log = {value: exponent for exponent, value in enumerate(powers)}

    probability = terminal_probability(33)
    script_hash = sha256(Path(__file__).resolve())
    manifest_hash = sha256(MANIFEST)
    receipts = []
    for parameter in sorted(support_11):
        inverse = field_power(parameter, MASK64 - 1)
        for parity_1 in sorted(support_11):
            data_index = power_log.get(field_multiply(parity_1, inverse))
            if data_index is None:
                continue
            local_values = (parameter, parameter, parity_1)
            local_words = []
            for local_value in local_values:
                codeword = encode_local(local_value)
                support = nibble_support(codeword)
                if support != 11:
                    raise RuntimeError("support-33 witness failed local authentication")
                local_words.append(
                    {
                        "message_hex": f"0x{local_value:016x}",
                        "codeword_hex": f"0x{codeword:032x}",
                        "nibble_support": support,
                    }
                )
            receipts.append(
                {
                    "schema": "riffle-dp-g4-counterexample-receipt-v1",
                    "candidate": "Riffle DP g=4@g0-v1",
                    "manifest_sha256": manifest_hash,
                    "evidence_label": "EXACT",
                    "family_id": (
                        f"terminal-w3-s33-one-i{data_index}-t{parameter:016x}"
                    ),
                    "outer_family": {
                        "block_weight": 3,
                        "category": "one data symbol and both parity symbols",
                        "parameterization": (
                            f"m_{data_index}=0x{parameter:016x}; all other data symbols are zero"
                        ),
                        "coefficient_relation": (
                            f"p_0=m_{data_index}; p_1=x^{data_index}m_{data_index}"
                        ),
                    },
                    "local_words": local_words,
                    "inner_event": {
                        "description": (
                            "All 33 nonzero outer nibbles occupy the final 47184 logical packet positions."
                        ),
                        "output_weight_upper": 188736,
                    },
                    "population": "1",
                    "event_probability": {
                        "numerator": str(probability.numerator),
                        "denominator": str(probability.denominator),
                    },
                    "first_moment_contribution": {
                        "exact_expression": "C(47184,33)/C(524352,33)",
                        "log2_interval": [float(x) for x in log2_interval(probability)],
                        "exceeds_2^-40": probability > Fraction(1, 1 << 40),
                    },
                    "authentication": {
                        "command": (
                            "python scripts/riffle_dp_g4_weighted_triples.py refutation "
                            "--low-message-file bch_g4_le13.bin"
                        ),
                        "source_sha256": script_hash,
                        "independent_checks": [
                            "Every local message was re-encoded from the BCH generator rows.",
                            "Every local word has exactly 11 nonzero nibbles.",
                            "The field ratio equals one declared coefficient x^i.",
                        ],
                    },
                }
            )

    if len(receipts) != 26:
        raise RuntimeError(f"expected 26 exact support-33 words, found {len(receipts)}")
    aggregate = len(receipts) * probability
    ledger = {
        "schema": "riffle-dp-g4-counterexample-ledger-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "receipt_count": len(receipts),
        "aggregate_first_moment_contribution": {
            "numerator": str(aggregate.numerator),
            "denominator": str(aggregate.denominator),
            "log2_interval": log2_interval(aggregate),
            "exceeds_2^-40": aggregate > Fraction(1, 1 << 40),
        },
        "receipts": receipts,
    }
    output.write_text(json.dumps(ledger, indent=2) + "\n")
    print(f"receipt={output}")
    print(f"exact_witnesses={len(receipts)}")
    print(
        "aggregate_log2="
        f"{ledger['aggregate_first_moment_contribution']['log2_interval']}"
    )
    print("status=EXACT_REFUTATION_FAMILY_BELOW_TARGET")


def fraction_from_receipt(receipt: dict, field: str) -> Fraction:
    value = receipt[field]
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def parse_classification(path: Path) -> dict[int, dict[str, int]]:
    result = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or not fields[0].startswith("support="):
            continue
        parsed = {}
        for field in fields:
            name, value = field.split("=", 1)
            parsed[name] = int(value)
        result[parsed.pop("support")] = parsed
    if set(result) != {33, 34, 35}:
        raise ValueError("classification receipt does not cover supports 33--35")
    return result


def write_g1_ledger(classification_file: Path, output: Path) -> None:
    classification = parse_classification(classification_file)
    word_counts = {
        support: sum(
            values[f"{category}_words"]
            for category in ("one", "pair_p0", "pair_p1", "triple")
        )
        for support, values in classification.items()
    }
    new_dirty_counts = {
        support: sum(
            values[f"{category}_new_dirty_classes"]
            for category in ("one", "pair_p0", "pair_p1", "triple")
        )
        for support, values in classification.items()
    }
    if word_counts != {33: 26, 34: 0, 35: 39}:
        raise RuntimeError(f"unexpected exact low word counts: {word_counts}")
    if new_dirty_counts != {33: 3, 34: 0, 35: 11}:
        raise RuntimeError(f"unexpected dirty-class counts: {new_dirty_counts}")

    dirty_classes = sum(new_dirty_counts.values())
    class_count = math.comb(DATA_SYMBOLS + 2, 3)
    clean_classes = class_count - dirty_classes
    clean_certificate = exact_single_krawtchouk_certificate(36, 42)
    support_33_certificate = exact_single_krawtchouk_certificate(33, 45)
    support_35_certificate = exact_single_krawtchouk_certificate(35, 43)
    support_39_certificate = exact_single_krawtchouk_certificate(39, 39)
    clean_bound = fraction_from_receipt(clean_certificate, "per_class_terminal_bound")
    support_33_bound = fraction_from_receipt(
        support_33_certificate, "per_class_terminal_bound"
    )
    support_35_bound = fraction_from_receipt(
        support_35_certificate, "per_class_terminal_bound"
    )
    support_39_bound = fraction_from_receipt(
        support_39_certificate, "per_class_terminal_bound"
    )
    upper_bound = (
        clean_classes * clean_bound
        + new_dirty_counts[33] * support_33_bound
        + new_dirty_counts[35] * support_35_bound
    )
    exact_lower_family = sum(
        count * terminal_probability(support)
        for support, count in word_counts.items()
    )
    target = Fraction(1, 1 << 40)
    all_remaining_support_39_bound = (
        new_dirty_counts[33] * support_33_bound
        + new_dirty_counts[35] * support_35_bound
        + clean_classes * support_39_bound
    )
    dirty_36_to_38_increment = clean_bound - support_39_bound
    maximum_additional_dirty_classes = int(
        (target - all_remaining_support_39_bound) // dirty_36_to_38_increment
    )

    ledger = {
        "schema": "riffle-dp-g4-g1-terminal-weight-three-ledger-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "manifest_sha256": sha256(MANIFEST),
        "engine_source_sha256": sha256(Path(__file__).resolve()),
        "scope": (
            "All block-weight-three outer words under the terminal-placement event."
        ),
        "coefficient_class_count": str(class_count),
        "exact_low_classification": {
            "word_counts": {str(key): value for key, value in word_counts.items()},
            "new_dirty_class_counts": {
                str(key): value for key, value in new_dirty_counts.items()
            },
            "dirty_classes": dirty_classes,
            "clean_classes": str(clean_classes),
            "classification_sha256": sha256(classification_file),
            "classifier_source_sha256": sha256(
                REPOSITORY / "scripts" / "classify_riffle_dp_g4_low_triples.cpp"
            ),
        },
        "proof_upper_bound": {
            "evidence_label": "RIGOROUS_BOUND",
            "partition": (
                "The 3 support-33 classes use the exact d=33,j=45 certificate; "
                "the 11 support-35 classes use d=35,j=43; every other class "
                "uses d=36,j=42."
            ),
            "numerator": str(upper_bound.numerator),
            "denominator": str(upper_bound.denominator),
            "log2_interval": log2_interval(upper_bound),
            "meets_2^-40": upper_bound <= target,
            "gap_above_target_bits_interval": log2_interval(upper_bound / target),
        },
        "refutation_lower_family": {
            "evidence_label": "EXACT",
            "description": (
                "Every authenticated support-33 through support-35 word under its terminal event."
            ),
            "numerator": str(exact_lower_family.numerator),
            "denominator": str(exact_lower_family.denominator),
            "log2_interval": log2_interval(exact_lower_family),
            "exceeds_2^-40": exact_lower_family > target,
        },
        "support_39_completion_gate": {
            "evidence_label": "RIGOROUS_BOUND",
            "all_previously_clean_classes_at_support_39_log2_interval":
                log2_interval(all_remaining_support_39_bound),
            "maximum_classes_with_minimum_support_36_to_38": str(
                maximum_additional_dirty_classes
            ),
            "interpretation": (
                "G1 passes if at most this many of the currently clean classes "
                "contain a word of support 36, 37, or 38. Such classes use the "
                "d=36 bound; all other classes use d=39."
            ),
        },
        "decision": "UNRESOLVED",
        "next_required_improvement": (
            "Upper-bound the number of coefficient classes containing a support-36, "
            "support-37, or support-38 word by the completion-gate threshold."
        ),
    }
    output.write_text(json.dumps(ledger, indent=2) + "\n")
    print(f"ledger={output}")
    print(f"proof_upper_log2={ledger['proof_upper_bound']['log2_interval']}")
    print(
        "proof_gap_bits="
        f"{ledger['proof_upper_bound']['gap_above_target_bits_interval']}"
    )
    print(
        "refutation_lower_log2="
        f"{ledger['refutation_lower_family']['log2_interval']}"
    )
    print("status=G1_UNRESOLVED_WITH_RIGOROUS_INTERVAL")


def parse_support38_classification(path: Path) -> dict:
    result: dict = {"layers": {}}
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0].startswith("support="):
            layer = {}
            for field in fields:
                name, value = field.split("=", 1)
                if value != "NOT_CLASSIFIED":
                    layer[name] = int(value)
            support = layer.pop("support")
            result["layers"][support] = layer
        elif "=" in fields[0]:
            name, value = fields[0].split("=", 1)
            result[name] = int(value) if value.isdigit() else value
    if result.get("status") != "EXACT_SUPPORT38_CLASSIFICATION":
        raise ValueError("support-38 classifier did not report exact completion")
    if set(result["layers"]) != set(range(33, 39)):
        raise ValueError("support-38 classification has incomplete layer coverage")
    return result


def write_g1_support38_ledger(classification_file: Path, output: Path) -> None:
    classification = parse_support38_classification(classification_file)
    layers = classification["layers"]
    if [classification[f"local_support_{support}"] for support in range(11, 15)] != [
        20,
        1_526,
        37_014,
        742_031,
    ]:
        raise RuntimeError("support-38 classifier used an unexpected local enumerator")

    exact_minimum_33 = sum(
        layers[33][f"{category}_new_dirty_classes"]
        for category in ("one", "pair_p1", "pair_p0", "triple")
    )
    exact_minimum_35 = sum(
        layers[35][f"{category}_new_dirty_classes"]
        for category in ("one", "pair_p1", "pair_p0", "triple")
    )
    exact_minimum_36 = sum(
        layers[36][f"{category}_new_dirty_classes"]
        for category in ("one", "pair_p1", "pair_p0", "triple")
    )
    exact_additive_37 = sum(
        layers[37][f"{category}_new_dirty_classes"]
        for category in ("pair_p0", "triple")
    )
    exact_additive_38 = sum(
        layers[38][f"{category}_new_dirty_classes"]
        for category in ("pair_p0", "triple")
    )
    if (
        exact_minimum_33,
        exact_minimum_35,
        exact_minimum_36,
        exact_additive_37,
        exact_additive_38,
    ) != (3, 11, 32_779, 0, 49_143):
        raise RuntimeError("support-38 class counts failed their frozen cross-check")

    one_class_count = DATA_SYMBOLS
    pair_class_count = math.comb(DATA_SYMBOLS, 2)
    repeated_class_count = one_class_count + pair_class_count
    classified_repeated = sum(
        layers[support].get("one_new_dirty_classes", 0)
        + layers[support].get("pair_p1_new_dirty_classes", 0)
        for support in range(33, 37)
    )
    conservative_support_37 = repeated_class_count - classified_repeated

    coefficient_class_count = math.comb(DATA_SYMBOLS + 2, 3)
    proof_bucket_counts = {
        33: exact_minimum_33,
        35: exact_minimum_35,
        36: exact_minimum_36,
        37: conservative_support_37,
        38: exact_additive_38,
    }
    proof_bucket_counts[39] = coefficient_class_count - sum(proof_bucket_counts.values())
    if sum(proof_bucket_counts.values()) != coefficient_class_count:
        raise RuntimeError("support-38 proof buckets do not partition the coefficient classes")

    degrees = {33: 45, 35: 43, 36: 42, 37: 41, 38: 40, 39: 39}
    per_class_bounds = {}
    contributions = {}
    for minimum_support, count in proof_bucket_counts.items():
        certificate = exact_single_krawtchouk_certificate(
            minimum_support, degrees[minimum_support]
        )
        bound = fraction_from_receipt(certificate, "per_class_terminal_bound")
        per_class_bounds[minimum_support] = bound
        contributions[minimum_support] = count * bound
    upper_bound = sum(contributions.values(), Fraction())
    target = Fraction(1, 1 << 40)

    exact_lower_family = 26 * terminal_probability(33) + 39 * terminal_probability(35)
    ledger = {
        "schema": "riffle-dp-g4-g1-terminal-weight-three-ledger-v2",
        "candidate": "Riffle DP g=4@g0-v1",
        "manifest_sha256": sha256(MANIFEST),
        "engine_source_sha256": sha256(Path(__file__).resolve()),
        "scope": (
            "All block-weight-three outer words under the terminal-placement event."
        ),
        "exact_classification": {
            "classification_sha256": sha256(classification_file),
            "classifier_source_sha256": sha256(
                REPOSITORY / "scripts" / "classify_riffle_dp_g4_support38.cpp"
            ),
            "local_enumerator_source_sha256": sha256(
                REPOSITORY / "scripts" / "enumerate_bch_g4_low_support.cpp"
            ),
            "local_message_file_sha256": sha256(REPOSITORY / "bch_g4_le15.bin"),
            "exact_minimum_class_counts": {
                "33": exact_minimum_33,
                "35": exact_minimum_35,
                "36": exact_minimum_36,
                "additive_37": exact_additive_37,
                "additive_38": exact_additive_38,
            },
            "coverage": {
                "one_data": "exact through support 36",
                "two_data_weighted_parity": "exact through support 36",
                "two_data_sum_parity": "exact through support 38",
                "three_data": "exact through support 38",
            },
        },
        "proof_partition": {
            str(support): {
                "class_count": str(proof_bucket_counts[support]),
                "interpretation": (
                    "conservative repeated-category remainder"
                    if support == 37
                    else "certified minimum-support bucket"
                ),
                "delsarte_degree": degrees[support],
                "contribution_log2_interval": log2_interval(contributions[support]),
            }
            for support in proof_bucket_counts
        },
        "proof_upper_bound": {
            "evidence_label": "RIGOROUS_BOUND",
            "numerator": str(upper_bound.numerator),
            "denominator": str(upper_bound.denominator),
            "log2_interval": log2_interval(upper_bound),
            "meets_2^-40": upper_bound <= target,
            "margin_below_target_bits_interval": log2_interval(target / upper_bound),
        },
        "refutation_lower_family": {
            "evidence_label": "EXACT",
            "description": (
                "All authenticated support-33 and support-35 words under their terminal events."
            ),
            "numerator": str(exact_lower_family.numerator),
            "denominator": str(exact_lower_family.denominator),
            "log2_interval": log2_interval(exact_lower_family),
            "exceeds_2^-40": exact_lower_family > target,
        },
        "decision": "PASS_TERMINAL_WEIGHT_THREE",
        "scope_limitation": (
            "This ledger closes only the terminal-placement event for outer block weight three. "
            "It does not bound nonterminal inner placements or outer block weights at least four."
        ),
    }
    if not ledger["proof_upper_bound"]["meets_2^-40"]:
        raise RuntimeError("support-38 ledger did not close terminal G1")
    output.write_text(json.dumps(ledger, indent=2) + "\n")
    print(f"ledger={output}")
    print(f"proof_upper_log2={ledger['proof_upper_bound']['log2_interval']}")
    print(
        "proof_margin_bits="
        f"{ledger['proof_upper_bound']['margin_below_target_bits_interval']}"
    )
    print(
        "refutation_lower_log2="
        f"{ledger['refutation_lower_family']['log2_interval']}"
    )
    print("status=PASS_TERMINAL_WEIGHT_THREE")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    proof = subparsers.add_parser("proof", help="emit an exact Delsarte receipt")
    proof.add_argument("--minimum-support", type=int, default=36)
    proof.add_argument("--degree", type=int, default=42)
    proof.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_g1_delsarte_d36_j42.json",
    )
    refutation = subparsers.add_parser(
        "refutation", help="emit exact support-33 terminal receipts"
    )
    refutation.add_argument("--low-message-file", type=Path, required=True)
    refutation.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_g1_support33_refutation.json",
    )
    ledger = subparsers.add_parser(
        "ledger", help="combine the exact low classification and dual certificates"
    )
    ledger.add_argument(
        "--classification-file",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_g1_low_classification.txt",
    )
    ledger.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_g1_terminal_ledger_pre38.json",
    )
    ledger38 = subparsers.add_parser(
        "ledger38", help="close terminal G1 from the support-38 classification"
    )
    ledger38.add_argument(
        "--classification-file",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_support38_classification.txt",
    )
    ledger38.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIRECTORY / "riffle_dp_g4_g1_terminal_ledger.json",
    )
    args = parser.parse_args()
    if args.mode == "proof":
        write_proof_receipt(args.minimum_support, args.degree, args.output)
    elif args.mode == "refutation":
        write_refutation_receipt(args.low_message_file, args.output)
    elif args.mode == "ledger":
        write_g1_ledger(args.classification_file, args.output)
    elif args.mode == "ledger38":
        write_g1_support38_ledger(args.classification_file, args.output)


if __name__ == "__main__":
    main()
