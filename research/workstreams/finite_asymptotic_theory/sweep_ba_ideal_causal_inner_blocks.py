#!/usr/bin/env python3
"""Sweep block size for the ideal random-convolution diagnostic."""

from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.special import logsumexp

from analyze_ba_ideal_causal_inner_prefix import (
    LN2,
    ba_log_local_prefix_counts,
    independent_rows_result,
    random_injection_log_local_prefix_counts,
)
from evaluate_ba_ideal_causal_inner_q1 import (
    log_average_for_weight,
    log_capped_row_union,
    region_average_tables,
)
from analyze_golay_ba_rm2sub_joint import expected_ba_log_spectrum


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ba_ideal_causal_g1_block_sweep.json"
OUTPUT_BITS = 2_119_680
BLOCK_SIZES = [240, 360, 480, 720, 960, 1440]


def reuse_aware_q1_margin(
    outer_bits: int,
    outer_rows: int,
    distance: int,
) -> float:
    exact_region, _ = region_average_tables(
        outer_bits,
        outer_rows,
        distance,
    )
    log_spectrum = expected_ba_log_spectrum(outer_bits)
    terms = []
    for weight in range(1, outer_bits + 1):
        log_multiplicity = float(log_spectrum[weight])
        if not math.isfinite(log_multiplicity):
            continue
        log_single = log_average_for_weight(
            weight,
            outer_bits,
            exact_region,
        )
        terms.append(
            log_multiplicity
            + log_capped_row_union(log_single, outer_rows)
        )
    return -float(logsumexp(terms)) / LN2


def main() -> None:
    rows = []
    for outer_bits in BLOCK_SIZES:
        if OUTPUT_BITS % outer_bits:
            raise ArithmeticError(f"{outer_bits} does not divide {OUTPUT_BITS}")
        outer_rows = OUTPUT_BITS // outer_bits
        distance = math.floor(0.11 * OUTPUT_BITS)
        random_result = independent_rows_result(
            "independent uniform random injections",
            random_injection_log_local_prefix_counts(
                outer_bits,
                outer_bits // 2,
            ),
            outer_rows,
            distance,
        )
        ba_result = independent_rows_result(
            "independently sampled Golay--BA-3 rows",
            ba_log_local_prefix_counts(outer_bits),
            outer_rows,
            distance,
        )
        rows.append(
            {
                "outer_bits": outer_bits,
                "outer_rows": outer_rows,
                "output_bits": OUTPUT_BITS,
                "parent_message_bits": OUTPUT_BITS // 2,
                "independent_random_all_message_margin_bits": (
                    random_result["margin_bits"]
                ),
                "independent_ba_all_message_margin_bits": (
                    ba_result["margin_bits"]
                ),
                "repeated_ba_reuse_aware_q1_margin_bits": (
                    reuse_aware_q1_margin(outer_bits, outer_rows, distance)
                ),
                "independent_ba_dominant_prefix_fraction": (
                    ba_result["dominant_prefix_fraction"]
                ),
                "independent_ba_dominant_suffix_fraction": (
                    ba_result["dominant_suffix_fraction"]
                ),
            }
        )

    result = {
        "schema": "ba-ideal-causal-g1-block-sweep-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "output_bits": OUTPUT_BITS,
            "bad_output_weight_at_most": math.floor(0.11 * OUTPUT_BITS),
            "inner": "invertible random lower-triangular Toeplitz convolution",
        },
        "rows": rows,
        "limitations": [
            "All arithmetic is nearest binary64.",
            "The all-message BA column uses independently sampled BA codes across outer rows.",
            "The reuse-aware column covers occupation one for one repeated BA code.",
            "Neither BA column alone is a full theorem for one randomly sampled and repeated code.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
