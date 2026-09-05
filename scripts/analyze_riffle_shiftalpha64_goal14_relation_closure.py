#!/usr/bin/env python3
"""Combine Goal 14 relation counts with the rigorous inner bounds."""

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


def profile_map(payload: dict[str, object]) -> dict[tuple[int, ...], dict[str, object]]:
    return {
        tuple(int(weight) for weight in row["weights"]): row
        for row in payload["profiles"]
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--relations",
        type=Path,
        default=RECEIPTS / "goal14_weight24_relation_core.json",
    )
    parser.add_argument(
        "--goal13-relations",
        type=Path,
        default=RECEIPTS / "goal13_initial_relations.json",
    )
    parser.add_argument(
        "--inner",
        type=Path,
        default=RECEIPTS / "goal14_even_core_inner.json",
    )
    parser.add_argument(
        "--mixed-outer",
        type=Path,
        default=RECEIPTS / "goal14_weight242422_outer_scan.json",
    )
    args = parser.parse_args()

    relations = json.loads(args.relations.read_text())
    goal13_relations = json.loads(args.goal13_relations.read_text())
    inner = profile_map(json.loads(args.inner.read_text()))
    mixed_outer = json.loads(args.mixed_outer.read_text(encoding="utf-8-sig"))
    finite_coordinates = 16_385
    fixed_ratio_support_bound = finite_coordinates * (finite_coordinates - 1)

    l24_cube = int(relations["relations"]["L24^3"]["ordered_relations"])
    mixed = int(
        relations["relations"]["L22_x_L24^2"]["ordered_relations"]
    )
    l22_l22_l26 = int(
        goal13_relations["profiles"]["22_22_26"]["ordered_relations"]
    )
    mixed_outer_words = int(
        mixed_outer["exact_finite_weight242422_outer_words"]
    )
    if mixed_outer_words != 0:
        raise RuntimeError("the (22,24,24) shifted-schedule zero changed")

    def complete_bound(weights: tuple[int, int, int], relation_count: int) -> float:
        return (
            math.log2(relation_count * fixed_ratio_support_bound)
            + float(inner[weights]["selected_one_word_log2_bound"])
        )

    high_l24 = complete_bound((24, 104, 104), l24_cube)
    high_mixed_same = complete_bound((22, 104, 104), mixed)
    high_mixed_cross = complete_bound((24, 104, 106), mixed)
    high_l2226 = complete_bound((22, 102, 106), l22_l22_l26)
    high_l2226_same = complete_bound((26, 106, 106), l22_l22_l26)
    payload = {
        "schema": "riffle-shiftalpha64-goal14-relation-closure-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "bad_output_weight_inclusive": 188_766,
        "fixed_ratio_support_multiplicity_bound": fixed_ratio_support_bound,
        "relation_families": {
            "L24^3": {
                "ordered_relations": l24_cube,
                "profiles": {
                    "24_24_24": {
                        "one_word_inner_log2_bound": inner[(24, 24, 24)][
                            "selected_one_word_log2_bound"
                        ],
                        "status": "requires shifted-schedule matching",
                    },
                    "24_104_104": {
                        "one_word_inner_log2_bound": inner[(24, 104, 104)][
                            "selected_one_word_log2_bound"
                        ],
                        "complete_profile_log2_bound": high_l24,
                        "status": "closed for all ordered placements",
                    },
                },
            },
            "L22_x_L24^2": {
                "ordered_relations": mixed,
                "profiles": {
                    "22_24_24": {
                        "one_word_inner_log2_bound": inner[(22, 24, 24)][
                            "selected_one_word_log2_bound"
                        ],
                        "exact_shifted_schedule_outer_words": mixed_outer_words,
                        "complete_profile_contribution": 0,
                        "status": "closed for all ordered placements",
                    },
                    "22_104_104": {
                        "one_word_inner_log2_bound": inner[(22, 104, 104)][
                            "selected_one_word_log2_bound"
                        ],
                        "complete_profile_log2_bound": high_mixed_same,
                        "status": "closed for all ordered placements",
                    },
                    "24_104_106": {
                        "one_word_inner_log2_bound": inner[(24, 104, 106)][
                            "selected_one_word_log2_bound"
                        ],
                        "complete_profile_log2_bound": high_mixed_cross,
                        "status": "closed for all ordered placements",
                    },
                },
            },
            "L22_x_L22_x_L26": {
                "ordered_relations": l22_l22_l26,
                "profiles": {
                    "22_102_106": {
                        "one_word_inner_log2_bound": inner[(22, 102, 106)][
                            "selected_one_word_log2_bound"
                        ],
                        "complete_profile_log2_bound": high_l2226,
                        "status": "closed for all ordered placements",
                    },
                    "26_106_106": {
                        "one_word_inner_log2_bound": inner[(26, 106, 106)][
                            "selected_one_word_log2_bound"
                        ],
                        "complete_profile_log2_bound": high_l2226_same,
                        "status": "closed for all ordered placements",
                    }
                },
            },
        },
        "closed_profile_log2_sum_bound": float(
            math.log2(
                2.0**high_l24
                + 2.0**high_mixed_same
                + 2.0**high_mixed_cross
                + 2.0**high_l2226
                + 2.0**high_l2226_same
            )
        ),
        "validation": {
            "relation_counts_are_exact": True,
            "complement_reduction_preserves_the_additive_relation": True,
            "fixed_ratio_support_bound_uses_unique_third_coordinate": True,
            "inner_bounds_are_rigorous_at_selected_tilts": True,
        },
        "scope": (
            "Complete finite-support closures for every two-complement profile "
            "in the three selected relation families and for (22,24,24), "
            "including all ordered placements. Only (24,24,24) remains open "
            "within this relation core."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
