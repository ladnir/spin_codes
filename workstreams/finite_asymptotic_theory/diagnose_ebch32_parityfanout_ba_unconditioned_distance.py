#!/usr/bin/env python3
"""Search the unconditioned EBCH32--PF31x33--BA-3 q=1 distance margin.

This is a nearest-binary64 diagnostic.  It does not provide an outward
certificate.  The search identifies a finite-distance target with enough
slack to justify translating the complete occupation cover.
"""

from __future__ import annotations

import argparse
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
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_unconditioned_distance_diagnostic.json"
B = 256
L = 8192
N = B * L
LOCAL_SPECTRUM = {
    0: 1,
    8: 620,
    12: 13_888,
    16: 36_518,
    20: 13_888,
    24: 620,
    32: 1,
}


def evaluate(spectrum, distance: int) -> dict[str, object]:
    active = one_active(
        spectrum=spectrum,
        outer_bits=B,
        outer_blocks=L,
        distance=distance,
        step_bits=128,
        state_bits=19,
        constituent_distance=48,
        activation=DEFAULT_ACTIVATION,
        live_spectrum_path=DEFAULT_LIVE_SPECTRUM,
        tilt_minimum=-12.0,
        tilt_maximum=-5.0,
        tilt_spacing=0.01,
    )
    return {
        "distance": distance,
        "relative_bad_weight": distance / N,
        "guaranteed_relative_distance": (distance + 1) / N,
        "q1_aggregate_log2_upper": active["aggregate_log2_upper"],
        "q1_margin_bits": active["aggregate_margin_bits"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--relative-distances",
        default="0.09,0.092,0.094,0.095,0.096,0.098,0.10",
        help="comma-separated bad-weight fractions",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    relative_distances = [
        float(value) for value in args.relative_distances.split(",") if value.strip()
    ]
    if any(not 0.0 < value < 0.5 for value in relative_distances):
        raise ValueError("every relative distance must lie strictly between 0 and 1/2")

    spectrum = expected_ba_log_spectrum(
        B, 32, 16, LOCAL_SPECTRUM, parity_fanout=(31, 33)
    )
    rows = [evaluate(spectrum, math.floor(value * N)) for value in relative_distances]
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-unconditioned-distance-diagnostic-v1",
        "status": "NEAREST_BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": 1 << 20,
            "row_length": B,
            "row_count": L,
            "output_length": N,
            "outer_sampling": "independent unconditioned rows",
        },
        "q1_results": rows,
        "limitations": [
            "Nearest-binary64 arithmetic selects a target only.",
            "The receipt covers occupation one, not the complete occupation union.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
