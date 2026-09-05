#!/usr/bin/env python3
"""Sweep RM2Sub state dimension at t=64 on two exact BCH spectra."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_inner_calibration import (  # noqa: E402
    evaluate_randomstep_case,
    tilt_grid,
)
from small_k_replay.evaluate_rm2sub_epoch_calibration import (  # noqa: E402
    outer_spectrum_path,
    rm2sub_case,
)


GENERATED = HERE / "rm2sub_calibration_constituents"
DEFAULT_JSON = HERE / "rm2sub_t64_state_calibration_bch_q1_d100.json"
DEFAULT_CSV = HERE / "rm2sub_t64_state_calibration_bch_q1_d100.csv"


def configuration(state_bits: int) -> dict[str, object]:
    stem = f"t64_s{state_bits}"
    selection_path = GENERATED / f"{stem}_selection.json"
    selection = json.loads(selection_path.read_text(encoding="utf-8"))["selected"]
    return {
        "name": stem,
        "step_bits": 64,
        "state_bits": state_bits,
        "a_distance": selection["minimum_A_distance"],
        "kernel_distance": selection["minimum_kernel_distance"],
        "kernel_weight_four": selection["weight_four_kernel_words"],
        "a": GENERATED / f"{stem}_a_spectrum.json",
        "b": GENERATED / f"{stem}_b_kernel_spectrum.json",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponents", type=int, nargs="+", default=[13, 15])
    parser.add_argument("--constituents", nargs="+", default=["ebch32", "ebch128"])
    parser.add_argument("--state-bits", type=int, nargs="+", default=[15, 16, 18, 19, 20])
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()
    distance = args.distance_numerator / args.distance_denominator
    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)
    configs = [configuration(state_bits) for state_bits in args.state_bits]

    cases: list[dict[str, object]] = []
    for key in args.constituents:
        constituent = CONSTITUENTS[key]
        spectrum = load_spectrum(constituent)
        spectrum_path = outer_spectrum_path(key, spectrum)
        for exponent in args.message_exponents:
            message_bits = 1 << exponent
            outer_rows = message_bits // constituent.dimension
            output_bits = outer_rows * constituent.block_bits
            bad_weight = math.floor(distance * output_bits)
            common = {
                "constituent": constituent.name,
                "block_bits": constituent.block_bits,
                "dimension": constituent.dimension,
                "message_exponent": exponent,
                "message_bits": message_bits,
                "outer_rows": outer_rows,
                "output_bits": output_bits,
                "bad_weight": bad_weight,
            }
            for config in configs:
                state_bits = int(config["state_bits"])
                for comparison, memory_bits in (
                    ("same-stored-state-bits", state_bits),
                    ("same-live-state-reset-rate", state_bits + 6),
                ):
                    print(
                        f"case,{key},kexp,{exponent},RSC,M,{memory_bits},"
                        f"reference_s,{state_bits}",
                        flush=True,
                    )
                    row = evaluate_randomstep_case(
                        constituent=constituent,
                        spectrum=spectrum,
                        message_exponent=exponent,
                        memory_bits=memory_bits,
                        bad_weight=bad_weight,
                        tilts=tilts,
                    )
                    row.update(
                        {
                            "configuration": f"m{memory_bits}",
                            "comparison": comparison,
                            "reference_rm2sub_state_bits": state_bits,
                            "persistence_exponent": float(memory_bits),
                            **common,
                        }
                    )
                    cases.append(row)
                print(
                    f"case,{key},kexp,{exponent},RM2Sub,{config['name']}",
                    flush=True,
                )
                row = rm2sub_case(
                    key=key,
                    spectrum_path=spectrum_path,
                    exponent=exponent,
                    distance=distance,
                    tilts=tilts,
                    config=config,
                )
                row.update(
                    {
                        "comparison": "rm2sub",
                        "reference_rm2sub_state_bits": state_bits,
                        **common,
                    }
                )
                cases.append(row)

    payload = {
        "schema": "rm2sub-t64-state-calibration-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "shell-separated occupation-one union-bound margin in bits",
        "integer_event": (
            "output weight at most floor((distance_numerator/"
            "distance_denominator)*N)"
        ),
        "probability_space": {
            "outer": "one fixed exact-spectrum BCH constituent reused in every row",
            "routing": (
                "independent uniform row-coordinate and transposed-region "
                "permutations"
            ),
            "randomstep": "independent random linear transitions at every bit",
            "rm2sub": (
                "one fixed, exactly audited A/B pair with BA=0 for each state "
                "dimension; independent nonzero state multipliers per epoch"
            ),
        },
        "parameters": {
            "step_bits": 64,
            "message_exponents": args.message_exponents,
            "constituents": args.constituents,
            "state_bits": args.state_bits,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "persistence_matching_rule": "RandomStepConv M=s+log2(64)=s+6",
            "tilt_count": len(tilts),
        },
        "cases": cases,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "Only occupation one is covered.",
            "The selected A/B samples are exactly audited but not proven optimal.",
            "Persistence matching equates a first-order reset exponent, not the full transition law.",
            "The comparison does not measure implementation cost or randomness cost.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "constituent", "block_bits", "dimension", "message_exponent",
        "message_bits", "outer_rows", "output_bits", "bad_weight", "inner",
        "configuration", "comparison", "reference_rm2sub_state_bits",
        "step_bits", "state_bits", "persistence_exponent", "epochs_per_region",
        "a_minimum_distance", "kernel_minimum_distance", "kernel_weight_four",
        "margin_bits", "dominant_weight", "dominant_log_surprisal",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cases)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
