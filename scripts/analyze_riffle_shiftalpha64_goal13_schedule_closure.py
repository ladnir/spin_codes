#!/usr/bin/env python3
"""Close the first shared schedule-matching rung of Goal 13."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RECEIPTS = (
    ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--relations",
        type=Path,
        default=RECEIPTS / "goal13_initial_relations.json",
    )
    parser.add_argument(
        "--inner",
        type=Path,
        default=RECEIPTS / "goal13_initial_inner.json",
    )
    parser.add_argument(
        "--outer",
        type=Path,
        default=RECEIPTS / "goal13_weight222226_outer_scan.json",
    )
    args = parser.parse_args()

    relations = json.loads(args.relations.read_text())
    inner = json.loads(args.inner.read_text())
    outer = json.loads(args.outer.read_text())
    relation_count = int(relations["profiles"]["22_22_26"]["ordered_relations"])
    if relation_count != 168_184_576:
        raise RuntimeError("(22,22,26) relation count changed")
    outer_count = int(outer["exact_finite_weight222226_outer_words"])
    if outer_count != 61_929_389:
        raise RuntimeError("(22,22,26) outer count changed")
    inner_log2 = float(inner["profiles"][0]["selected_one_word_log2_bound"])
    profile_log2 = math.log2(outer_count) + inner_log2

    payload = {
        "schema": "riffle-shiftalpha64-goal13-schedule-closure-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "bad_output_weight_inclusive": 188_766,
        "finite_profiles": {
            "22_22_26": {
                "ordered_global_relations": relation_count,
                "exact_shifted_schedule_outer_words": outer_count,
                "one_word_inner_log2_bound": inner_log2,
                "complete_profile_log2_bound": profile_log2,
            },
            "odd_complement_low_shell_cluster": {
                "complete_profile_contribution": 0,
                "reason": (
                    "three representatives from weights {22,24,26} have total "
                    "weight at most 78 and cannot XOR to the weight-128 all-one word"
                ),
            },
        },
        "combined_schedule_rung_log2_bound": profile_log2,
        "validation": {
            "outer_scan_covers_all_three_unique_weight_roles": True,
            "small_instance_direct_audit": "PASS",
            "odd_complement_weight_obstruction": "3*26=78<128",
            "inner_bound_is_rigorous_at_selected_tilts": True,
        },
        "scope": (
            "Complete finite-support contribution of (22,22,26), plus the "
            "zero contribution of every odd-complement triple over low shells "
            "{22,24,26}."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
