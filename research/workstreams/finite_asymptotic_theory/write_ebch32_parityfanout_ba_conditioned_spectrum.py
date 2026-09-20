#!/usr/bin/env python3
"""Write the analyzer-facing conditioned EBCH32--PF31x33--BA spectrum."""

from __future__ import annotations

import json
import math
from pathlib import Path

from certify_ebch32_parityfanout_ba_setup import (
    B,
    LOWER_WEIGHT,
    UPPER_WEIGHT,
    conditioned_spectrum_upper,
)


OUTPUT = Path(__file__).resolve().parent / (
    "ebch32_parityfanout31x33_ba3_B256_conditioned_spectrum_upper.json"
)


def main() -> None:
    spectrum, tail, good = conditioned_spectrum_upper()
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-b256-conditioned-spectrum-upper-v1",
        "status": "OUTWARD_BINARY64_EXPECTATION_UPPER_FROM_MARKOV_CONDITIONING",
        "outer_bits": B,
        "outer_dimension": B // 2,
        "conditioning_event": (
            f"no nonzero final BA word outside weights {LOWER_WEIGHT}..{UPPER_WEIGHT}"
        ),
        "good_event_probability_lower_hex": good.hex(),
        "tail_expected_word_count_upper_hex": tail.hex(),
        "spectrum": [
            {
                "weight": weight,
                "log2_expected_multiplicity": (
                    0.0
                    if weight == 0
                    else math.log2(spectrum[weight])
                    if spectrum[weight] > 0.0
                    else None
                ),
                "expected_multiplicity_upper_hex": (
                    spectrum[weight].hex() if spectrum[weight] > 0.0 else None
                ),
            }
            for weight in range(B + 1)
        ],
        "scope": (
            "The hexadecimal fields are outward expected-spectrum bounds. "
            "The log2 fields are nearest-binary64 analyzer inputs and serve "
            "only to select witnesses."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
