#!/usr/bin/env python3
"""Evaluate one nested t=128 RM2Sub chain with the fixed RM(4,9) outer."""

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
from small_k_replay.evaluate_inner_calibration import tilt_grid  # noqa: E402
from small_k_replay.evaluate_rm2sub_epoch_calibration import (  # noqa: E402
    outer_spectrum_path,
    rm2sub_case as q1_case,
)
from small_k_replay.evaluate_rm2sub_q2_calibration import (  # noqa: E402
    rm2sub_case as q2_case,
)


def configuration(chain_dir: Path, state_bits: int) -> dict[str, object]:
    stem = f"t128_s{state_bits}"
    selection_path = chain_dir / f"{stem}_selection.json"
    selected = json.loads(selection_path.read_text(encoding="utf-8"))["selected"]
    return {
        "name": stem,
        "step_bits": 128,
        "state_bits": state_bits,
        "a_distance": selected["minimum_A_distance"],
        "kernel_distance": selected["minimum_kernel_distance"],
        "kernel_weight_four": selected["weight_four_kernel_words"],
        "selection": selection_path,
        "a": chain_dir / f"{stem}_a_spectrum.json",
        "b": chain_dir / f"{stem}_b_kernel_spectrum.json",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chain-dir", type=Path, required=True)
    parser.add_argument("--state-bits", type=int, nargs="+", default=[12, 13, 14, 15, 16])
    parser.add_argument("--message-exponent", type=int, default=16)
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument("--skip-q2", action="store_true")
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--csv-output", type=Path, required=True)
    args = parser.parse_args()

    constituent = CONSTITUENTS["rm49"]
    spectrum = load_spectrum(constituent)
    spectrum_path = outer_spectrum_path("rm49", spectrum)
    message_bits = 1 << args.message_exponent
    outer_rows = message_bits // constituent.dimension
    output_bits = outer_rows * constituent.block_bits
    distance = args.distance_numerator / args.distance_denominator
    bad_weight = math.floor(distance * output_bits)
    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)
    common = {
        "constituent": constituent.name,
        "message_exponent": args.message_exponent,
        "message_bits": message_bits,
        "outer_rows": outer_rows,
        "output_bits": output_bits,
        "bad_weight": bad_weight,
    }
    cases = []
    for state_bits in args.state_bits:
        config = configuration(args.chain_dir, state_bits)
        if outer_rows % 128:
            raise ValueError("128-bit epoch does not divide the transposed region")
        print(f"case,s{state_bits},Q1", flush=True)
        row = q1_case(
            key="rm49",
            spectrum_path=spectrum_path,
            exponent=args.message_exponent,
            distance=distance,
            tilts=tilts,
            config=config,
        )
        row.update({**common, "occupation": 1})
        cases.append(row)
        if not args.skip_q2:
            print(f"case,s{state_bits},Q2", flush=True)
            row = q2_case(
                key="rm49",
                exponent=args.message_exponent,
                config=config,
                bad_weight=bad_weight,
                tilts=tilts,
            )
            row.update({**common, "occupation": 2})
            cases.append(row)

    manifest_path = args.chain_dir / "NESTED_FAMILY_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload = {
        "schema": "rm2sub-nested-t128-rm49-screen-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "question": "How do Q1 and Q2 change under literal nested extension of A_s?",
        "construction": (
            "one fixed RM(4,9) outer repeated in every row and one nested "
            "unselected RM2Sub chain"
        ),
        "chain_manifest": str(manifest_path),
        "chain_seed": manifest["parameters"]["seed"],
        "probability_space": {
            "outer": "one fixed authenticated RM(4,9) constituent repeated in every row",
            "routing": "independent uniform row and region permutations",
            "inner": (
                "one deterministic nested A/B chain from the recorded seed; "
                "independent nonzero epoch scalars"
            ),
        },
        "parameters": {
            "state_bits": args.state_bits,
            "message_exponent": args.message_exponent,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "tilt_minimum": args.tilt_minimum,
            "tilt_maximum": args.tilt_maximum,
            "tilt_spacing": args.tilt_spacing,
            "tilt_count": len(tilts),
            "q2_enabled": not args.skip_q2,
        },
        "cases": cases,
        "limitations": [
            "The numerical bounds use nearest binary64 arithmetic.",
            "A finite witness grid is used.",
            "The result covers at most occupations one and two.",
            "One nested chain does not establish typical behavior over chain seeds.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "configuration",
        "state_bits",
        "step_bits",
        "a_minimum_distance",
        "kernel_minimum_distance",
        "kernel_weight_four",
        "occupation",
        "margin_bits",
        "dominant_weight",
        "dominant_first_weight",
        "dominant_second_weight",
        "dominant_log_surprisal",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cases)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
