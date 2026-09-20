#!/usr/bin/env python3
"""Count the first two Goal 13 low-shell/complement relations exactly."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as bch  # noqa: E402
from analyze_ebch128_adjacent_triples import (  # noqa: E402
    ALL_ONE_MESSAGE,
    DEFAULT_MESSAGES,
    classify_minimum_shell,
    write_raw_ratio_stream_222224,
)
from analyze_ebch128_weight22_triples import read_messages  # noqa: E402
from analyze_riffle_bch_alpha_joint_spectrum import encode_weights  # noqa: E402
from enumerate_ebch128_weight22_affine import AFFINE_GROUP_SIZE  # noqa: E402


ALL_ONE_CODEWORD = (1 << 128) - 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--messages", type=Path, default=DEFAULT_MESSAGES)
    parser.add_argument(
        "--raw-ratio-2226",
        type=Path,
        help="write one little-endian uint64 ratio per ordered (22,22,26) relation",
    )
    args = parser.parse_args()

    messages = read_messages(args.messages)
    messages_array = np.asarray(messages, dtype=np.uint64)
    (
        representatives,
        message_by_codeword,
        coordinate_to_point,
        point_to_coordinate,
        multiplication,
    ) = classify_minimum_shell(messages)
    complete = set(message_by_codeword)

    twisted_total = 0
    weight26_total = 0
    rows = []
    weight26_partners_by_representative: list[list[int]] = []
    witnesses: dict[str, list[dict[str, str]]] = {
        "106_106_106": [],
        "22_22_26": [],
    }
    for representative in representatives:
        representative_message = message_by_codeword[representative]
        twisted_partners = [
            other
            for other in complete
            if representative ^ other ^ ALL_ONE_CODEWORD in complete
        ]
        target_weights = encode_weights(
            messages_array ^ np.uint64(representative_message)
        )
        weight26_partners = messages_array[target_weights == 26]
        weight26_partners_by_representative.append(
            [int(value) for value in weight26_partners]
        )
        twisted_total += AFFINE_GROUP_SIZE * len(twisted_partners)
        weight26_total += AFFINE_GROUP_SIZE * int(weight26_partners.size)
        rows.append(
            {
                "representative_message_hex": hex(representative_message),
                "twisted_minimum_partners": len(twisted_partners),
                "weight26_xor_partners": int(weight26_partners.size),
            }
        )

        for other in twisted_partners:
            if len(witnesses["106_106_106"]) >= 4:
                break
            third = representative ^ other ^ ALL_ONE_CODEWORD
            witnesses["106_106_106"].append(
                {
                    "low_first_message_hex": hex(representative_message),
                    "low_second_message_hex": hex(message_by_codeword[other]),
                    "low_third_message_hex": hex(message_by_codeword[third]),
                }
            )
        for other_message in weight26_partners[:4]:
            if len(witnesses["22_22_26"]) >= 4:
                break
            witnesses["22_22_26"].append(
                {
                    "first_weight22_message_hex": hex(representative_message),
                    "second_weight22_message_hex": hex(int(other_message)),
                    "weight26_message_hex": hex(
                        representative_message ^ int(other_message)
                    ),
                }
            )

    if twisted_total % 6:
        raise RuntimeError("twisted ordered triple count is not divisible by six")
    if weight26_total % 2:
        raise RuntimeError("(22,22,26) ordered relation count is not even")
    if bch.encode_message(ALL_ONE_MESSAGE) != ALL_ONE_CODEWORD:
        raise RuntimeError("all-one BCH message changed")

    raw_ratio_count = None
    if args.raw_ratio_2226:
        raw_ratio_count = write_raw_ratio_stream_222224(
            args.raw_ratio_2226,
            representatives=representatives,
            partners_by_representative=weight26_partners_by_representative,
            message_by_codeword=message_by_codeword,
            coordinate_to_point=coordinate_to_point,
            point_to_coordinate=point_to_coordinate,
            multiplication=multiplication,
        )
        if raw_ratio_count != weight26_total:
            raise RuntimeError("raw (22,22,26) ratio stream has the wrong size")

    payload = {
        "schema": "ebch128-goal13-low-complement-relations-v1",
        "code": "extended primitive binary BCH [128,64,22]",
        "minimum_shell_size": len(messages),
        "affine_orbits": len(representatives),
        "affine_orbit_size": AFFINE_GROUP_SIZE,
        "profiles": {
            "106_106_106": {
                "low_shell_identity": "x+y+z=e with x,y,z in L22",
                "ordered_relations": twisted_total,
                "unordered_relations": twisted_total // 6,
                "log2_ordered_relations": (
                    math.log2(twisted_total) if twisted_total else None
                ),
            },
            "22_22_26": {
                "low_shell_identity": "x,y in L22 and wt(B(x+y))=26",
                "ordered_relations": weight26_total,
                "unordered_light_pairs": weight26_total // 2,
                "log2_ordered_relations": (
                    math.log2(weight26_total) if weight26_total else None
                ),
                "raw_ratio_stream": (
                    {
                        "output": str(args.raw_ratio_2226),
                        "records": raw_ratio_count,
                        "binary_record": "little-endian uint64 ratio",
                    }
                    if args.raw_ratio_2226
                    else None
                ),
            },
        },
        "representative_rows": rows,
        "witnesses": witnesses,
        "validation": {
            "minimum_shell_partitioned_into_affine_orbits": True,
            "all_one_message_reencoded": True,
            "target_weights_computed_by_linear_reencoding": True,
            "symmetry_divisibility_checks": True,
            "raw_ratio_stream_matches_relation_count": (
                True if args.raw_ratio_2226 else "NOT_RUN"
            ),
        },
        "scope": (
            "Exact global relation counts. Outer coefficient schedules and "
            "inner accumulator probabilities are not included."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
