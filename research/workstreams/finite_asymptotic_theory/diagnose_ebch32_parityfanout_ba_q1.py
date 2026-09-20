#!/usr/bin/env python3
"""Nearest-binary64 tilt search for the EBCH32--PF31x33--BA q=1 proof."""

from __future__ import annotations

import json
import math
from pathlib import Path

from compare_ebch32_ba_k20 import expected_ba_log_spectrum
from evaluate_golay_ba_rm2sub_finite import one_active
from evaluate_random_outer_rm2sub_one_active import (
    DEFAULT_ACTIVATION,
    DEFAULT_LIVE_SPECTRUM,
)


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_q1_d11_diagnostic.json"
B = 256
L = 8192
DISTANCE = 230_686
LOCAL_SPECTRUM = {
    0: 1,
    8: 620,
    12: 13_888,
    16: 36_518,
    20: 13_888,
    24: 620,
    32: 1,
}


def main() -> None:
    spectrum = expected_ba_log_spectrum(
        B, 32, 16, LOCAL_SPECTRUM, parity_fanout=(31, 33)
    )
    active = one_active(
        spectrum=spectrum,
        outer_bits=B,
        outer_blocks=L,
        distance=DISTANCE,
        step_bits=128,
        state_bits=19,
        constituent_distance=48,
        activation=DEFAULT_ACTIVATION,
        live_spectrum_path=DEFAULT_LIVE_SPECTRUM,
        tilt_minimum=-10.0,
        tilt_maximum=-6.0,
        tilt_spacing=0.01,
    )
    tail_logs = [
        float(spectrum[weight])
        for weight in list(range(1, 24)) + list(range(233, B + 1))
        if math.isfinite(float(spectrum[weight]))
    ]
    maximum = max(tail_logs)
    tail = maximum + math.log(sum(math.exp(value - maximum) for value in tail_logs))
    good = 1.0 - math.exp(tail)
    permitted_rows = [
        row for row in active["weight_rows"] if 24 <= int(row["outer_weight"]) <= 232
    ]
    permitted_terms = [
        float(row["pointwise_log2_upper"]) * math.log(2.0)
        for row in permitted_rows
    ]
    maximum_term = max(permitted_terms)
    aggregate = maximum_term + math.log(
        sum(math.exp(value - maximum_term) for value in permitted_terms)
    ) - math.log(good)
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-b256-q1-d11-diagnostic-v1",
        "status": "NEAREST_BINARY64_WITNESS_SEARCH",
        "parameters": {
            "row_length": B,
            "row_count": L,
            "output_length": B * L,
            "distance": DISTANCE,
            "permitted_outer_weights": [24, 232],
        },
        "conditioned_q1": {
            "aggregate_log2_upper": aggregate / math.log(2.0),
            "aggregate_margin_bits": -aggregate / math.log(2.0),
            "weight_rows": permitted_rows,
        },
        "limitations": [
            "Nearest-binary64 values select witnesses only.",
            "No numerical bound in this receipt is a certificate.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["conditioned_q1"] | {"weight_rows": len(permitted_rows)}, indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
