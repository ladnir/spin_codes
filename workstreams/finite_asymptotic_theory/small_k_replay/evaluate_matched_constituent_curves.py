#!/usr/bin/env python3
"""Compare matched fixed-size outer constituents under one memory schedule.

For k=2^e, set M=e+2.  This is the schedule M=ceil(log2(k))+2 and reproduces
M=22 at k=2^20.  At each block size, compare the authenticated structured
constituent with one uniform full-rank random constituent sampled once and
reused in every row.  The target distance is 10 percent.

Every value is the occupation-one transfer margin.  It is not a complete
all-occupation distance certificate.  Arithmetic is nearest binary64.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

from small_k_replay.evaluate_exact_spectra_q1_phase import (
    CONSTITUENTS,
    evaluate_case,
    load_spectrum,
    tilt_grid,
)
from small_k_replay.evaluate_reused_random_constituent_q1 import (
    expected_spectrum,
)


DEFAULT_JSON = HERE / "matched_constituent_k_margin_d100.json"
DEFAULT_CSV = HERE / "matched_constituent_k_margin_d100.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--message-exponents", type=int, nargs="+", default=list(range(8, 21))
    )
    parser.add_argument("--memory-offset", type=int, default=2)
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-5.0)
    parser.add_argument("--fine-step", type=float, default=0.1)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    tilts = tilt_grid(args)
    families = []
    for key in ("ebch32", "ebch128", "rm49"):
        structured = CONSTITUENTS[key]
        families.append(
            (
                f"{structured.name} exact",
                "structured",
                structured,
                load_spectrum(structured),
            )
        )
        random_name = (
            f"random full-rank [{structured.block_bits},{structured.dimension}]"
        )
        random_constituent = type(structured)(
            name=random_name,
            block_bits=structured.block_bits,
            dimension=structured.dimension,
            minimum_distance=1,
            spectrum_path=Path(),
        )
        families.append(
            (
                f"{random_name} reused",
                "random",
                random_constituent,
                expected_spectrum(structured.block_bits),
            )
        )

    rows = []
    for label, model, constituent, spectrum in families:
        for exponent in args.message_exponents:
            message_bits = 1 << exponent
            if message_bits < constituent.dimension:
                continue
            memory_bits = exponent + args.memory_offset
            print(
                f"case,{label},kexp,{exponent},memory,{memory_bits}",
                flush=True,
            )
            row = evaluate_case(
                constituent=constituent,
                spectrum=spectrum,
                message_exponent=exponent,
                memory_bits=memory_bits,
                distance_numerator=args.distance_numerator,
                distance_denominator=args.distance_denominator,
                tilts=tilts,
            )
            row["series"] = label
            row["outer_model"] = model
            row["block_bits"] = constituent.block_bits
            row["dimension"] = constituent.dimension
            rows.append(row)

    payload = {
        "schema": "matched-fixed-constituent-q1-k-margin-curve-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "occupation-one union-bound margin in bits",
        "memory_schedule": {
            "formula": "M(k)=ceil(log2(k))+2 for the evaluated powers of two",
            "offset": args.memory_offset,
            "anchor": {"message_bits": 1 << 20, "memory_bits": 22},
        },
        "distance": {
            "fraction": args.distance_numerator / args.distance_denominator,
            "cutoff_rule": "D=ceil(fraction*N)",
        },
        "probability_space": {
            "structured_outer": "one fixed authenticated constituent reused in every row",
            "random_outer": "one uniform full-rank matched-size constituent sampled once and reused in every row",
            "routing": "independent uniform row-coordinate and transposed-region permutations",
            "inner": "independent RandomStepConv maps sampled once and shared by all messages",
        },
        "cases": rows,
        "limitations": [
            "These curves cover occupation one only and are not full distance-certificate margins.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "Random constituents are compared only with the structured constituent of the same block size and dimension.",
        ],
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    fieldnames = [
        "series",
        "outer_model",
        "block_bits",
        "dimension",
        "message_exponent",
        "message_bits",
        "memory_bits",
        "distance_cutoff",
        "relative_distance_lower",
        "q1_margin_bits_diagnostic",
        "dominant_weight",
        "dominant_log_surprisal",
    ]
    args.csv_output.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
