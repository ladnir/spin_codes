#!/usr/bin/env python3
"""Evaluate the transfer-weighted Q=1 spectrum of one sparse-EA constituent.

The outer constituent is C=A E, where every one of the B rows of E is an
independent uniform r-subset of the K input coordinates and A is the
zero-initialized accumulator.  One sampled full-rank constituent is reused in
all L outer rows.

For Q=1, linearity of expectation permits averaging the actual SPIN transfer
functional over the one sampled constituent.  This is different from
requiring simultaneous upper caps on every realized weight shell.

The transfer factors are recovered from the existing B=512 RandomStepConv-M22
Q=1 diagnostic.  The arithmetic here is nearest numpy.longdouble and is a
parameter diagnostic, not an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
DEFAULT_TRANSFER = WORKSTREAM / "single_random_constituent_B512_q1_s22_d109.json"
DEFAULT_OUTPUT = HERE / "ea_K256_B512_combined_q1_degree_sweep.json"
LOG2 = np.log(np.longdouble(2))


def krawtchouk(length: int, degree: int, weight: int) -> int:
    return sum(
        (-1) ** j
        * math.comb(weight, j)
        * math.comb(length - weight, degree - j)
        for j in range(
            max(0, degree - (length - weight)), min(degree, weight) + 1
        )
    )


def activation_probabilities(message_bits: int, right_degree: int) -> np.ndarray:
    denominator = np.longdouble(math.comb(message_bits, right_degree))
    return np.asarray(
        [
            (np.longdouble(1) - np.longdouble(
                krawtchouk(message_bits, right_degree, weight)
            ) / denominator) / np.longdouble(2)
            for weight in range(1, message_bits + 1)
        ],
        dtype=np.longdouble,
    )


def expected_accumulated_spectrum(
    message_bits: int, output_bits: int, right_degree: int
) -> np.ndarray:
    """Return E[A_w(AE)] for w=0,...,B in long-double arithmetic."""
    q = activation_probabilities(message_bits, right_degree)[:, None]
    one_minus_q = np.longdouble(1) - q
    zero = np.zeros((message_bits, output_bits + 1), dtype=np.longdouble)
    one = np.zeros_like(zero)
    zero[:, 0] = np.asarray(
        [math.comb(message_bits, weight) for weight in range(1, message_bits + 1)],
        dtype=np.longdouble,
    )
    for completed in range(output_bits):
        following_zero = np.zeros_like(zero)
        following_one = np.zeros_like(one)
        active = slice(0, completed + 1)
        following_zero[:, active] = (
            one_minus_q * zero[:, active] + q * one[:, active]
        )
        following_one[:, 1 : completed + 2] = (
            q * zero[:, active] + one_minus_q * one[:, active]
        )
        zero, one = following_zero, following_one
    return np.sum(zero + one, axis=0)


def expected_accumulated_spectrum_truncated(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    maximum_weight: int,
) -> np.ndarray:
    """Return exact low-weight recurrence through ``maximum_weight``.

    Accumulated weight is monotone during the recursion, so discarding states
    above the cutoff cannot later change a retained coefficient.
    """
    if maximum_weight >= output_bits:
        return expected_accumulated_spectrum(
            message_bits, output_bits, right_degree
        )
    q = activation_probabilities(message_bits, right_degree)[:, None]
    one_minus_q = np.longdouble(1) - q
    zero = np.zeros((message_bits, maximum_weight + 1), dtype=np.longdouble)
    one = np.zeros_like(zero)
    zero[:, 0] = np.asarray(
        [math.comb(message_bits, weight) for weight in range(1, message_bits + 1)],
        dtype=np.longdouble,
    )
    for _ in range(output_bits):
        following_zero = one_minus_q * zero + q * one
        following_one = np.zeros_like(one)
        following_one[:, 1:] = q * zero[:, :-1] + one_minus_q * one[:, :-1]
        zero, one = following_zero, following_one
    return np.sum(zero + one, axis=0)


def expected_accumulated_spectrum_log2_truncated(
    message_bits: int,
    output_bits: int,
    right_degree: int,
    maximum_weight: int,
    rescale_interval: int = 16,
) -> np.ndarray:
    """Return log2 E[A_w] without materializing central binomial counts.

    Each input-weight row carries a probability law initialized at one.  The
    rows are rescaled independently during the accumulator recursion.  The
    final aggregation adds the binomial message multiplicities in log space.
    """
    q = np.asarray(
        activation_probabilities(message_bits, right_degree), dtype=np.float64
    )[:, None]
    one_minus_q = 1.0 - q
    zero = np.zeros((message_bits, maximum_weight + 1), dtype=np.float64)
    one = np.zeros_like(zero)
    zero[:, 0] = 1.0
    log_scales = np.zeros(message_bits, dtype=np.float64)
    for completed in range(output_bits):
        following_zero = one_minus_q * zero + q * one
        following_one = np.zeros_like(one)
        following_one[:, 1:] = q * zero[:, :-1] + one_minus_q * one[:, :-1]
        zero, one = following_zero, following_one
        if (completed + 1) % rescale_interval == 0 or completed + 1 == output_bits:
            scales = np.maximum(np.max(zero, axis=1), np.max(one, axis=1))
            positive = scales > 0
            zero[positive] /= scales[positive, None]
            one[positive] /= scales[positive, None]
            log_scales[positive] += np.log(scales[positive])

    input_log_counts = np.asarray(
        [
            math.lgamma(message_bits + 1)
            - math.lgamma(weight + 1)
            - math.lgamma(message_bits - weight + 1)
            for weight in range(1, message_bits + 1)
        ],
        dtype=np.float64,
    )
    result = np.full(maximum_weight + 1, -math.inf, dtype=np.float64)
    for output_weight in range(maximum_weight + 1):
        masses = zero[:, output_weight] + one[:, output_weight]
        positive = masses > 0
        if not np.any(positive):
            continue
        terms = input_log_counts[positive] + log_scales[positive] + np.log(masses[positive])
        maximum = np.max(terms)
        result[output_weight] = (
            maximum + math.log(float(np.sum(np.exp(terms - maximum))))
        ) / math.log(2.0)
    return result


def log2_sum(values: np.ndarray) -> float:
    positive = values[values > 0]
    if positive.size == 0:
        return -math.inf
    maximum = np.max(np.log2(positive))
    return float(maximum + np.log2(np.sum(np.exp2(np.log2(positive) - maximum))))


def load_transfer(path: Path, block_bits: int) -> tuple[np.ndarray, dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates = [
        block for block in payload["blocks"] if int(block["block_bits"]) == block_bits
    ]
    if len(candidates) != 1:
        raise ValueError(f"transfer receipt does not contain unique B={block_bits} block")
    block = candidates[0]
    outer_rows = int(block["outer_rows"])
    factors = np.zeros(block_bits + 1, dtype=np.longdouble)
    for row in block["weight_rows"]:
        weight = int(row["weight"])
        # pointwise = log2(L) + log2(E A_w for a random injection) + log2 transfer
        inner_log2 = (
            np.longdouble(row["pointwise_log2_upper"])
            - np.log2(np.longdouble(outer_rows))
            - np.longdouble(row["expected_multiplicity_log2"])
        )
        factors[weight] = np.exp2(inner_log2)
    return factors, {
        "file": path.name,
        "schema": payload.get("schema"),
        "status": payload.get("status"),
        "outer_rows": outer_rows,
        "output_bits": int(block["output_bits"]),
        "distance_cutoff": int(block["distance_cutoff"]),
        "memory_bits": int(block["memory_bits"]),
    }


def evaluate_degree(
    *,
    message_bits: int,
    output_bits: int,
    right_degree: int,
    rank_attempts: int,
    transfer: np.ndarray,
    outer_rows: int,
) -> dict[str, object]:
    spectrum = expected_accumulated_spectrum(
        message_bits, output_bits, right_degree
    )
    kernel_mean = spectrum[0]
    rank_success_lower = max(np.longdouble(0), np.longdouble(1) - kernel_mean)
    abort_upper = min(np.longdouble(1), kernel_mean) ** rank_attempts
    weighted_terms = spectrum[1:] * transfer[1:]
    unconditional_q1 = np.longdouble(outer_rows) * np.sum(weighted_terms)
    conditional_q1 = (
        np.longdouble(np.inf)
        if rank_success_lower == 0
        else unconditional_q1 / rank_success_lower
    )
    setup_or_q1 = abort_upper + conditional_q1
    dominant = int(np.argmax(weighted_terms)) + 1

    def bits(value: np.longdouble) -> float | None:
        if not np.isfinite(value) or value <= 0:
            return None
        return float(-np.log2(value))

    return {
        "right_degree": right_degree,
        "expected_left_degree": output_bits * right_degree / message_bits,
        "sparse_map_xors": output_bits * max(0, right_degree - 1),
        "accumulator_xors": output_bits - 1,
        "total_xors_per_constituent": output_bits * max(0, right_degree - 1)
        + output_bits
        - 1,
        "unconditional_kernel_expected_count": float(kernel_mean),
        "rank_success_lower": float(rank_success_lower),
        "rank_test_abort_upper": float(abort_upper),
        "unconditional_q1_upper": float(unconditional_q1),
        "unconditional_q1_margin_bits": bits(unconditional_q1),
        "conditional_q1_upper": float(conditional_q1),
        "conditional_q1_margin_bits": bits(conditional_q1),
        "setup_abort_or_q1_upper": float(setup_or_q1),
        "setup_abort_or_q1_margin_bits": bits(setup_or_q1),
        "dominant_outer_weight": dominant,
        "dominant_weighted_term_log2": (
            None
            if weighted_terms[dominant - 1] <= 0
            else float(np.log2(weighted_terms[dominant - 1]))
        ),
        "expected_spectrum_log2": [
            None if value <= 0 else float(np.log2(value)) for value in spectrum
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        default=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 29, 33],
    )
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--transfer", type=Path, default=DEFAULT_TRANSFER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if any(not 1 <= degree <= args.message_bits for degree in args.degrees):
        parser.error("every degree must lie in [1,message_bits]")

    transfer, transfer_meta = load_transfer(args.transfer, args.output_bits)
    rows = []
    for degree in args.degrees:
        print(f"degree,{degree}", flush=True)
        rows.append(
            evaluate_degree(
                message_bits=args.message_bits,
                output_bits=args.output_bits,
                right_degree=degree,
                rank_attempts=args.rank_attempts,
                transfer=transfer,
                outer_rows=int(transfer_meta["outer_rows"]),
            )
        )
    passing = [
        row for row in rows
        if row["setup_abort_or_q1_margin_bits"] is not None
        and float(row["setup_abort_or_q1_margin_bits"]) >= 40
    ]
    payload = {
        "schema": "ea-transfer-weighted-combined-q1-sweep-v1",
        "status": "BINARY_LONGDOUBLE_DIAGNOSTIC",
        "probability_space": {
            "outer_attempt": "sample B independent uniform weight-r rows E, then set C=A E",
            "outer_setup": "sample independent attempts, accept the first full-rank C, and abort after the stated number of failures",
            "outer_reuse": "reuse the one accepted constituent in every outer row",
            "routing_and_inner": "the frozen Q=1 transfer averages uniform row and region permutations and RandomStepConv-M22 setup",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "rank_attempts": args.rank_attempts,
            "degrees": args.degrees,
        },
        "transfer_source": transfer_meta,
        "first_tested_degree_closing_q1_at_40_bits": (
            None if not passing else passing[0]["right_degree"]
        ),
        "rows": rows,
        "proved_vs_diagnostic": [
            "The one-word sparse-map law and accumulator recursion are exact identities.",
            "Linearity of expectation is sufficient for occupation Q=1 even though one constituent is reused.",
            "Conditioning divides by the lower bound 1-E[nonzero kernel size]; sixteen-attempt abort is bounded by the sixteenth power of that first moment.",
            "The numerical values use nearest long-double arithmetic and transfer factors recovered from a binary64 diagnostic; they are not outward certificates.",
            "No occupation Q>=2 statement follows from this receipt.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "first_tested_degree_closing_q1_at_40_bits": payload["first_tested_degree_closing_q1_at_40_bits"],
        "rows": [{
            "degree": row["right_degree"],
            "xors": row["total_xors_per_constituent"],
            "margin_bits": row["setup_abort_or_q1_margin_bits"],
            "dominant_weight": row["dominant_outer_weight"],
            "kernel_mean": row["unconditional_kernel_expected_count"],
        } for row in rows],
    }, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
