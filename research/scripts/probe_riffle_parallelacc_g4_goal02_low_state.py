#!/usr/bin/env python3
"""Diagnostic lower families for the worst ParallelAcc support-33 multiset."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
MOMENT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal02_support33_moment.json"
OUTPUT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal02_low_state_probe.json"
PACKET_COUNT = 524352
SUPPORT = 33
DISTANCE = 188766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    return [
        value
        for packet in range(32)
        if (value := (word >> (4 * packet)) & 15) != 0
    ]


def selected_profile() -> tuple[str, Counter[int]]:
    source = json.loads(SOURCE.read_text())
    for receipt in source["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        final_state = 0
        for value in values:
            final_state ^= value
        if final_state == 0:
            return receipt["family_id"], Counter(values)
    raise RuntimeError("support-33 source no longer has a final-XOR-zero word")


def order_counts_by_zero_prefix(counts: Counter[int], cap: int) -> np.ndarray:
    values = tuple(sorted(counts))
    maxima = tuple(counts[value] for value in values)
    strides = []
    state_count = 1
    for maximum in maxima:
        strides.append(state_count)
        state_count *= maximum + 1

    dynamic = np.zeros((state_count, SUPPORT + 1), dtype=np.float64)
    dynamic[0, 0] = 1.0
    for index in range(state_count):
        if not dynamic[index].any():
            continue
        state = 0
        for value, maximum, stride in zip(values, maxima, strides):
            used = (index // stride) % (maximum + 1)
            if used & 1:
                state ^= value
        for value, maximum, stride in zip(values, maxima, strides):
            used = (index // stride) % (maximum + 1)
            next_state = state ^ value
            if used == maximum or next_state.bit_count() > cap:
                continue
            shift = int(next_state == 0)
            dynamic[index + stride, shift:] += (
                dynamic[index, : SUPPORT + 1 - shift] * (maximum - used)
            )
    return dynamic[-1]


def placement_log_probability(*, zero_prefixes: int, cap: int) -> float:
    nonzero_prefixes = SUPPORT - zero_prefixes
    zero_slots = zero_prefixes + 1
    extra_positions = PACKET_COUNT - SUPPORT
    budget = (DISTANCE - cap * nonzero_prefixes) // cap
    if budget < 0:
        return -math.inf

    allocated = np.arange(budget + 1, dtype=np.float64)
    log_nonzero_compositions = (
        gammaln(allocated + nonzero_prefixes)
        - gammaln(nonzero_prefixes)
        - gammaln(allocated + 1)
    )
    zero_positions = extra_positions - allocated
    log_zero_compositions = (
        gammaln(zero_positions + zero_slots)
        - gammaln(zero_slots)
        - gammaln(zero_positions + 1)
    )
    log_bad_placements = logsumexp(log_nonzero_compositions + log_zero_compositions)
    log_all_placements = (
        gammaln(PACKET_COUNT + 1)
        - gammaln(SUPPORT + 1)
        - gammaln(PACKET_COUNT - SUPPORT + 1)
    )
    return float(log_bad_placements - log_all_placements)


def main() -> None:
    family_id, counts = selected_profile()
    rows = []
    for cap in (2, 3, 4):
        order_counts = order_counts_by_zero_prefix(counts, cap)
        terms = []
        nonzero_rows = []
        for zero_prefixes, count in enumerate(order_counts):
            if not count:
                continue
            log_order_probability = math.log(count) - gammaln(SUPPORT + 1)
            log_placement_probability = placement_log_probability(
                zero_prefixes=zero_prefixes,
                cap=cap,
            )
            term = log_order_probability + log_placement_probability
            terms.append(term)
            nonzero_rows.append(
                {
                    "zero_prefixes": zero_prefixes,
                    "log2_order_probability": log_order_probability / math.log(2),
                    "log2_placement_probability": log_placement_probability / math.log(2),
                    "log2_contribution": term / math.log(2),
                }
            )
        aggregate = float(logsumexp(terms) / math.log(2)) if terms else None
        rows.append(
            {
                "maximum_allowed_nonzero_prefix_weight": cap,
                "aggregate_log2_lower_bound": aggregate,
                "zero_prefix_rows": nonzero_rows,
            }
        )

    payload = {
        "schema": "riffle-parallelacc-g4-goal02-low-state-probe-v1",
        "candidate": "Riffle ParallelAcc g=4",
        "evidence_label": "DIAGNOSTIC_LOWER_BOUND_FAMILY",
        "source_sha256": digest(SOURCE),
        "moment_receipt_sha256": digest(MOMENT),
        "outer_family": family_id,
        "packet_value_counts": {str(value): count for value, count in sorted(counts.items())},
        "rows": rows,
        "event": (
            "Every active-order prefix has weight at most c. The total extra gap "
            "length after nonzero prefixes is at most floor((d-c(h-m0))/c)."
        ),
        "scope_limitation": (
            "The order counts use binary64 dynamic programming and placement sums "
            "use binary64 logarithms. The event is sufficient but not necessary."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in rows:
        print(
            f"cap={row['maximum_allowed_nonzero_prefix_weight']} "
            f"aggregate_log2_lower_bound={row['aggregate_log2_lower_bound']:.12f}"
        )
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_LOWER_BOUND_FAMILY")


if __name__ == "__main__":
    main()
