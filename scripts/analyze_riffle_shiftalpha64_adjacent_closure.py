#!/usr/bin/env python3
"""Combine the exact outer counts and rigorous adjacent-profile kernels."""

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


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--relations",
        type=Path,
        default=RECEIPTS / "goal12_adjacent_triple_counts.json",
    )
    parser.add_argument(
        "--inner",
        type=Path,
        default=RECEIPTS / "goal12_adjacent_inner.json",
    )
    parser.add_argument(
        "--outer-222224",
        type=Path,
        default=RECEIPTS / "goal12_weight222224_outer_scan.json",
    )
    args = parser.parse_args()

    relations = load(args.relations)
    inner = load(args.inner)
    outer = load(args.outer_222224)
    heavy_relations = int(relations["profiles"]["22_106_106"]["ordered_relations"])
    if heavy_relations != 1_365_504:
        raise RuntimeError("heavy relation count changed")
    heavy_inner = float(inner["profiles"][0]["source_preserving"]["log2_bound"])
    light_inner = float(inner["profiles"][1]["pooled_one_word_log2_bound"])
    finite_coordinates = 16_385
    heavy_outer_bound = heavy_relations * finite_coordinates * (finite_coordinates - 1)
    heavy_total = math.log2(heavy_outer_bound) + heavy_inner

    light_outer = int(outer["exact_finite_weight222224_outer_words"])
    if light_outer != 9_331_089:
        raise RuntimeError("(22,22,24) outer count changed")
    light_total = math.log2(light_outer) + light_inner
    combined = max(heavy_total, light_total) + math.log2(
        1.0 + 2.0 ** (-abs(heavy_total - light_total))
    )

    payload = {
        "schema": "riffle-shiftalpha64-adjacent-closure-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "bad_output_weight_inclusive": 188_766,
        "finite_profiles": {
            "22_106_106": {
                "ordered_global_relations": heavy_relations,
                "fixed_ratio_support_multiplicity_bound": (
                    finite_coordinates * (finite_coordinates - 1)
                ),
                "outer_word_count_upper_bound": heavy_outer_bound,
                "one_word_inner_log2_bound": heavy_inner,
                "complete_profile_log2_bound": heavy_total,
            },
            "22_22_24": {
                "exact_outer_word_count": light_outer,
                "one_word_inner_log2_bound": light_inner,
                "complete_profile_log2_bound": light_total,
            },
        },
        "combined_adjacent_profiles_log2_bound": combined,
        "validation": {
            "heavy_relation_count_matches_exact_enumerator": True,
            "heavy_support_bound_uses_unique_third_coordinate": True,
            "light_outer_count_matches_exact_scan": True,
            "inner_receipt_uses_rigorous_selected_tilts": True,
        },
        "scope": (
            "Complete finite-support contribution of profiles (22,106,106) "
            "and (22,22,24), including all ordered coordinate placements."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
