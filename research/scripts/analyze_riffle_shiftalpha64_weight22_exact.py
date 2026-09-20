#!/usr/bin/env python3
"""Exactly scan shifted alpha lags on all minimum EBCH messages.

The input is the sorted little-endian uint64 list produced by
``enumerate_ebch128_weight22_affine.py``.  For every requested exponent e,
this script computes

    |L_22 intersection gamma^e L_22|,

where L_22 is the complete set of messages whose extended BCH [128,64]
encoding has binary weight 22.  It prints a JSON receipt and does not modify
the workspace.

Two 24-bit occupancy filters make exact membership checks sparse.  They have
no false negatives; every surviving candidate is then checked by binary
search in the complete sorted list.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as dp  # noqa: E402
from analyze_riffle_bch_alpha_joint_spectrum import encode_weights  # noqa: E402


EXPECTED_COUNT = 243_840
EXPECTED_WEIGHT = 22
FILTER_BITS = 24
FILTER_SIZE = 1 << FILTER_BITS
FILTER_MASK = np.uint64(FILTER_SIZE - 1)
DEFAULT_LAG_START = 64
DEFAULT_LAG_COUNT = dp.USER_BLOCKS


def read_messages(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if len(data) % 8:
        raise ValueError("minimum-message file has a partial uint64 record")
    messages = np.frombuffer(data, dtype="<u8").astype(np.uint64, copy=True)
    if messages.size != EXPECTED_COUNT:
        raise ValueError(
            f"expected {EXPECTED_COUNT} minimum messages, found {messages.size}"
        )
    if np.any(messages[1:] <= messages[:-1]):
        raise ValueError("minimum-message file is not strictly increasing")
    if np.any(encode_weights(messages) != EXPECTED_WEIGHT):
        raise ValueError("minimum-message file contains a non-weight-22 BCH word")
    return messages


def multiply_x(values: np.ndarray) -> None:
    """Multiply a uint64 array by gamma=X in place."""
    top = values >> np.uint64(63)
    values <<= np.uint64(1)
    values ^= top * np.uint64(dp.FIELD_REDUCTION)


def validate_vector_multiply(messages: np.ndarray, exponent: int) -> None:
    probes = messages[np.linspace(0, messages.size - 1, 17, dtype=np.int64)].copy()
    expected = probes.copy()
    for _ in range(exponent):
        multiply_x(expected)
    for source, target in zip(probes.tolist(), expected.tolist()):
        scalar = dp.field_multiply(dp.field_power(2, exponent), int(source))
        if int(target) != scalar:
            raise RuntimeError("vectorized gamma multiplication failed scalar check")


def build_filters(messages: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    low_filter = np.zeros(FILTER_SIZE, dtype=np.bool_)
    middle_filter = np.zeros(FILTER_SIZE, dtype=np.bool_)
    low_filter[(messages & FILTER_MASK).astype(np.int64)] = True
    middle_filter[((messages >> np.uint64(FILTER_BITS)) & FILTER_MASK).astype(np.int64)] = True
    if not np.all(low_filter[(messages & FILTER_MASK).astype(np.int64)]):
        raise RuntimeError("low occupancy filter has a false negative")
    if not np.all(
        middle_filter[
            ((messages >> np.uint64(FILTER_BITS)) & FILTER_MASK).astype(np.int64)
        ]
    ):
        raise RuntimeError("middle occupancy filter has a false negative")
    return low_filter, middle_filter


def exact_membership_count(candidates: np.ndarray, sorted_messages: np.ndarray) -> int:
    if not candidates.size:
        return 0
    locations = np.searchsorted(sorted_messages, candidates)
    inside = locations < sorted_messages.size
    return int(
        np.count_nonzero(
            inside
            & (sorted_messages[np.minimum(locations, sorted_messages.size - 1)] == candidates)
        )
    )


def scan(
    messages: np.ndarray,
    lag_start: int,
    lag_count: int,
    progress_every: int,
) -> dict[str, object]:
    low_filter, middle_filter = build_filters(messages)
    values = messages.copy()
    for _ in range(lag_start):
        multiply_x(values)

    total_pairs = 0
    filter_candidates = 0
    hit_rows: list[dict[str, int]] = []
    maximum = 0
    started = time.perf_counter()
    for lag_offset in range(lag_count):
        low = (values & FILTER_MASK).astype(np.int64)
        middle = (
            (values >> np.uint64(FILTER_BITS)) & FILTER_MASK
        ).astype(np.int64)
        survivors = values[low_filter[low] & middle_filter[middle]]
        hits = exact_membership_count(survivors, messages)
        filter_candidates += int(survivors.size)
        total_pairs += hits
        maximum = max(maximum, hits)
        if hits:
            hit_rows.append({"exponent": lag_start + lag_offset, "ordered_pairs": hits})
        multiply_x(values)
        if progress_every and (lag_offset + 1) % progress_every == 0:
            elapsed = time.perf_counter() - started
            print(
                f"scanned {lag_offset + 1}/{lag_count} exponents in {elapsed:.1f}s",
                file=sys.stderr,
                flush=True,
            )

    elapsed = time.perf_counter() - started
    return {
        "exponent_start_inclusive": lag_start,
        "exponent_end_exclusive": lag_start + lag_count,
        "exponents_scanned": lag_count,
        "sources_per_exponent": int(messages.size),
        "source_exponent_pairs": int(messages.size) * lag_count,
        "occupancy_filter_survivors": filter_candidates,
        "exact_weight22_to_weight22_ordered_pairs": total_pairs,
        "exponents_with_hits": len(hit_rows),
        "maximum_pairs_at_one_exponent": maximum,
        "hit_rows": hit_rows,
        "elapsed_seconds": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--messages", type=Path, required=True)
    parser.add_argument("--lag-start", type=int, default=DEFAULT_LAG_START)
    parser.add_argument("--lag-count", type=int, default=DEFAULT_LAG_COUNT)
    parser.add_argument("--progress-every", type=int, default=1024)
    args = parser.parse_args()
    if args.lag_start < 0 or args.lag_count < 1:
        raise ValueError("lag interval must be nonnegative and nonempty")

    messages = read_messages(args.messages)
    validate_vector_multiply(messages, args.lag_start)
    result = scan(messages, args.lag_start, args.lag_count, args.progress_every)
    payload = {
        "schema": "riffle-shiftalpha64-weight22-exact-v1",
        "evidence_label": "EXACT_COMPLETE_MINIMUM_SHELL_INTERSECTION",
        "construction": {
            "name": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
            "alpha_schedule": "alpha_i = gamma^(64+i), 0 <= i < 16384",
            "field": "GF(2^64) with X^64+X^4+X^3+X+1",
            "outer": "extended primitive binary BCH [128,64,22]",
        },
        "validation": {
            "minimum_message_count": int(messages.size),
            "expected_A_22": EXPECTED_COUNT,
            "strictly_sorted_and_unique": True,
            "all_bch_weights_equal_22": True,
            "vector_field_arithmetic_matches_scalar_probes": True,
            "occupancy_filters_have_no_false_negatives": True,
        },
        "scan": result,
        "scope": (
            "Exact for every minimum-weight BCH message and every exponent in "
            "the stated interval. Occupancy filters only discard values that "
            "cannot be members; all survivors receive an exact binary-search check."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
