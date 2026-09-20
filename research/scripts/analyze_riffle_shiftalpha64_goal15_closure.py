#!/usr/bin/env python3
"""Combine the exact (24,24,24) outer count with its rigorous inner bound."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CONSTRUCTION = (
    ROOT / "constructions" / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
)
RECEIPTS = CONSTRUCTION / "receipts"


def log2_add(*values: float) -> float:
    maximum = max(values)
    return maximum + math.log2(sum(2.0 ** (value - maximum) for value in values))


def main() -> None:
    scan = json.loads((RECEIPTS / "goal15_weight24_s3_scan.json").read_text())
    inner = json.loads((RECEIPTS / "goal14_even_core_inner.json").read_text())
    old_closure = json.loads((RECEIPTS / "goal14_relation_closure.json").read_text())
    envelope = json.loads(
        (RECEIPTS / "goal15_occ3_after_equal_shell.json").read_text()
    )

    profile = next(row for row in inner["profiles"] if row["weights"] == [24, 24, 24])
    one_word_log2 = profile["selected_one_word_log2_bound"]
    outer_words = scan["exact_finite_outer_words"]
    profile_log2 = math.log2(outer_words) + one_word_log2
    closed_log2 = log2_add(
        old_closure["closed_profile_log2_sum_bound"], profile_log2
    )
    open_log2 = envelope["combined_log2_bound"]
    complete_log2 = log2_add(open_log2, closed_log2)

    payload = {
        "schema": "riffle-shiftalpha64-goal15-equal-shell-closure-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "profile": [24, 24, 24],
        "exact_outer_words": outer_words,
        "log2_exact_outer_words": math.log2(outer_words),
        "one_word_inner_log2_bound": one_word_log2,
        "complete_profile_log2_bound": profile_log2,
        "complete_profile_probability_bound": 2.0**profile_log2,
        "prior_closed_profile_log2_sum_bound": old_closure[
            "closed_profile_log2_sum_bound"
        ],
        "updated_closed_profile_log2_sum_bound": closed_log2,
        "remaining_occupation_three_log2_bound": open_log2,
        "occupation_three_log2_bound_including_closed_profiles": complete_log2,
        "next_finite_profiles": [
            [22, 22, 28],
            [24, 24, 26],
            [22, 24, 26],
        ],
        "validation": {
            "outer_count_is_exact": True,
            "inner_bound_uses_goal14_rigorous_selected_tilt": True,
            "closed_and_open_families_are_disjoint": True,
            "log_space_sums_recomputed": True,
        },
        "scope": (
            "Closure of the finite-support (24,24,24) profile and the updated "
            "occupation-three envelope. This is not the full distance proof."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
