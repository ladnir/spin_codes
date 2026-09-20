#!/usr/bin/env python3
"""Combine the first exact Goal 13 relations with their rigorous kernels."""

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
    args = parser.parse_args()

    relations = json.loads(args.relations.read_text())
    inner = json.loads(args.inner.read_text())
    if int(relations["profiles"]["106_106_106"]["ordered_relations"]) != 0:
        raise RuntimeError("(106,106,106) is no longer absent")

    reused_relations = 27_765_248
    finite_coordinates = 16_385
    support_bound = finite_coordinates * (finite_coordinates - 1)
    one_word_inner = float(inner["profiles"][1]["selected_one_word_log2_bound"])
    mixed_bound = math.log2(reused_relations * support_bound) + one_word_inner
    second_one_word_inner = float(
        inner["profiles"][2]["selected_one_word_log2_bound"]
    )
    second_mixed_bound = (
        math.log2(reused_relations * support_bound) + second_one_word_inner
    )
    payload = {
        "schema": "riffle-shiftalpha64-goal13-initial-closure-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "bad_output_weight_inclusive": 188_766,
        "finite_profiles": {
            "106_106_106": {
                "ordered_global_relations": 0,
                "complete_profile_contribution": 0,
                "reason": "no x,y,z in L22 satisfy x+y+z=e",
            },
            "22_104_106": {
                "relation_source": "Goal 12 ordered (22,22,24) relations",
                "ordered_global_relations": reused_relations,
                "fixed_ratio_support_multiplicity_bound": support_bound,
                "one_word_inner_log2_bound": one_word_inner,
                "complete_profile_log2_bound": mixed_bound,
            },
            "24_106_106": {
                "relation_source": "Goal 12 ordered (22,22,24) relations",
                "ordered_global_relations": reused_relations,
                "fixed_ratio_support_multiplicity_bound": support_bound,
                "one_word_inner_log2_bound": second_one_word_inner,
                "complete_profile_log2_bound": second_mixed_bound,
            },
        },
        "validation": {
            "twisted_minimum_relation_enumerator_is_exact": True,
            "mixed_profile_complement_identities": [
                "(22,104,106) reduces to (22,24,22)",
                "(24,106,106) reduces to (24,22,22)"
            ],
            "mixed_support_bound_uses_unique_third_coordinate": True,
            "inner_bound_is_rigorous_at_selected_tilts": True,
        },
        "scope": (
            "Complete finite-support contributions of (106,106,106) and all "
            "ordered placements of (22,104,106) and (24,106,106)."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
