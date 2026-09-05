#!/usr/bin/env python3
"""Write native-product sections for all L22 x L24 x L24 ratios.

Each section stores one inverse denominator vector followed by every numerator
vector for that affine orbit.  A native PCLMUL converter performs the final
field products.  This layout avoids repeated Python field multiplication and
stores each denominator vector only once.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
import time
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

from analyze_ebch128_adjacent_triples import (  # noqa: E402
    affine_message_actions,
    transformed_messages_fast,
)
from analyze_ebch128_weight22_triples import batch_inverse  # noqa: E402
from analyze_ebch128_weight24_relation_core import (  # noqa: E402
    CONSTRUCTION,
    DEFAULT_WEIGHT24,
    EXPECTED_SHELL_SIZES,
    classify_affine_orbits,
    encode_halves,
    read_messages,
)
from enumerate_ebch128_weight22_affine import coordinate_maps  # noqa: E402


DEFAULT_RELATIONS = CONSTRUCTION / "receipts" / "goal14_weight24_relation_core.json"


def partner_indices(
    rep_low: int,
    rep_high: int,
    low: np.ndarray,
    high: np.ndarray,
    *,
    chunk_size: int,
) -> np.ndarray:
    """Return L24 words whose XOR with the representative has weight 22."""
    pieces = []
    rep_low64 = np.uint64(rep_low)
    rep_high64 = np.uint64(rep_high)
    for begin in range(0, low.size, chunk_size):
        end = min(begin + chunk_size, low.size)
        weights = np.bitwise_count(low[begin:end] ^ rep_low64)
        weights += np.bitwise_count(high[begin:end] ^ rep_high64)
        local = np.flatnonzero(weights == 22)
        if local.size:
            pieces.append(local.astype(np.int64) + begin)
    if not pieces:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(pieces)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--weight24", type=Path, default=DEFAULT_WEIGHT24)
    parser.add_argument("--relations", type=Path, default=DEFAULT_RELATIONS)
    parser.add_argument("--chunk-size", type=int, default=1 << 20)
    args = parser.parse_args()
    if args.chunk_size <= 0:
        raise ValueError("chunk size must be positive")

    relation_receipt = json.loads(args.relations.read_text(encoding="utf-8-sig"))
    expected_mass = int(
        relation_receipt["relations"]["L22_x_L24^2"]["ordered_relations"]
    )
    expected_rows = {
        row["representative_message_hex"]: (
            int(row["orbit_size"]),
            int(row["partner_count"]),
        )
        for row in relation_receipt["relations"]["L22_x_L24^2"][
            "from_L24_orbit_rows"
        ]
    }

    started = time.perf_counter()
    messages = read_messages(
        args.weight24, expected=EXPECTED_SHELL_SIZES[24]
    )
    low, high = encode_halves(messages)
    coordinate_to_point, point_to_coordinate, multiplication = coordinate_maps()
    actions = affine_message_actions(
        coordinate_to_point, point_to_coordinate, multiplication
    )
    representatives = classify_affine_orbits(
        messages, low, high, weight=24, actions=actions
    )

    written = 0
    inverse_values = 0
    pair_seeds = 0
    with args.output.open("wb") as destination:
        for representative_index, representative in enumerate(representatives):
            first_all = transformed_messages_fast(
                representative["codeword_low"]
                | (representative["codeword_high"] << 64),
                actions,
            )
            _, action_indices = np.unique(first_all, return_index=True)
            first = first_all[action_indices]
            if first.size != representative["orbit_size"]:
                raise RuntimeError("affine section has the wrong orbit size")
            inverses = np.asarray(batch_inverse(first.tolist()), dtype=np.uint64)
            partners = partner_indices(
                representative["codeword_low"],
                representative["codeword_high"],
                low,
                high,
                chunk_size=args.chunk_size,
            )
            expected = expected_rows[hex(representative["message"])]
            if expected != (representative["orbit_size"], int(partners.size)):
                raise RuntimeError("mixed relation receipt changed")

            destination.write(struct.pack("<II", int(first.size), int(partners.size)))
            inverses.astype("<u8", copy=False).tofile(destination)
            inverse_values += int(inverses.size)

            for partner_index in partners:
                partner_codeword = int(low[partner_index]) | (
                    int(high[partner_index]) << 64
                )
                second = transformed_messages_fast(partner_codeword, actions)[
                    action_indices
                ]
                if np.any(second == 0):
                    raise RuntimeError("mixed relation produced a zero numerator")
                second.astype("<u8", copy=False).tofile(destination)
                written += int(second.size)
                pair_seeds += 1

            if (representative_index + 1) % 16 == 0 or (
                representative_index + 1 == len(representatives)
            ):
                print(
                    f"mixed ratios: {representative_index + 1}/"
                    f"{len(representatives)} orbits, {written}/{expected_mass} records",
                    file=sys.stderr,
                    flush=True,
                )

    if written != expected_mass:
        raise RuntimeError(f"raw mixed ratio mass changed: {written}")
    expected_bytes = 8 * len(representatives) + 8 * (
        expected_mass + inverse_values
    )
    if args.output.stat().st_size != expected_bytes:
        raise RuntimeError("mixed ratio-product stream has the wrong byte length")
    payload = {
        "schema": "ebch128-goal14-mixed-ratio-products-v1",
        "profile": "22_24_24",
        "ratio": "second L24 message divided by first L24 message",
        "binary_format": (
            "per orbit: uint32 orbit_size, uint32 partner_count, "
            "orbit_size uint64 denominator inverses, then partner_count "
            "vectors of orbit_size uint64 numerators"
        ),
        "ordered_relations": written,
        "stored_inverse_values": inverse_values,
        "affine_pair_seeds": pair_seeds,
        "affine_orbits": len(representatives),
        "output": str(args.output),
        "output_bytes": args.output.stat().st_size,
        "elapsed_seconds": time.perf_counter() - started,
        "validation": {
            "orbit_sections_remove_stabilizer_duplicates": True,
            "orbit_partner_counts_match_relation_receipt": True,
            "ratio_mass_matches_exact_relation_count": True,
            "all_numerators_and_denominators_are_nonzero": True,
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
