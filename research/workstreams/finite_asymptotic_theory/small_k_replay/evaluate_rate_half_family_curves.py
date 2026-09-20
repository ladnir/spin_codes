#!/usr/bin/env python3
"""Compare rate-half BCH-derived, RM, and random outer constituents.

For k=2^e, use M=e+2, output length N=2k, and D=ceil(0.10N).
Every structured constituent is fixed and repeated in every outer row.  Each
random control samples one uniform full-rank [B,B/2] constituent once and
repeats it in every row.  The computation is the shell-separated Q=1 transfer
bound with nearest binary64 arithmetic; it is not an all-occupation distance
certificate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    Constituent,
    evaluate_case,
    load_spectrum,
    tilt_grid,
)
from small_k_replay.evaluate_reused_random_constituent_q1 import (  # noqa: E402
    expected_spectrum,
)


DEFAULT_JSON = HERE / "rate_half_family_k_margin_d100.json"
DEFAULT_CSV = HERE / "rate_half_family_k_margin_d100.csv"
STRUCTURED = {
    "bch": ("ebch8", "ebch32", "xbch64", "ebch128"),
    "rm": ("rm13", "rm25", "rm37", "rm49"),
}
RANDOM_BLOCKS = (8, 16, 32, 64, 128, 256, 512, 1024)


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
    if not 0 < args.distance_numerator < args.distance_denominator:
        parser.error("distance fraction must lie in (0,1)")

    tilts = tilt_grid(args)
    families: list[tuple[str, str, Constituent, dict[int, int | float]]] = []
    spectrum_sources = []
    for family, keys in STRUCTURED.items():
        for key in keys:
            constituent = CONSTITUENTS[key]
            spectrum = load_spectrum(constituent)
            label = f"{constituent.name} exact"
            families.append((label, family, constituent, spectrum))
            spectrum_sources.append(
                {
                    "series": label,
                    "path": str(constituent.spectrum_path),
                    "sha256": hashlib.sha256(
                        constituent.spectrum_path.read_bytes()
                    ).hexdigest(),
                }
            )
    for block_bits in RANDOM_BLOCKS:
        dimension = block_bits // 2
        constituent = Constituent(
            name=f"random full-rank [{block_bits},{dimension}]",
            block_bits=block_bits,
            dimension=dimension,
            minimum_distance=1,
            spectrum_path=Path(),
        )
        label = f"random full-rank [{block_bits},{dimension}] reused"
        families.append((label, "random", constituent, expected_spectrum(block_bits)))

    rows = []
    for label, family, constituent, spectrum in families:
        for exponent in args.message_exponents:
            if (1 << exponent) < constituent.dimension:
                continue
            memory_bits = exponent + args.memory_offset
            print(f"case,{label},kexp,{exponent},memory,{memory_bits}", flush=True)
            row = evaluate_case(
                constituent=constituent,
                spectrum=spectrum,
                message_exponent=exponent,
                memory_bits=memory_bits,
                distance_numerator=args.distance_numerator,
                distance_denominator=args.distance_denominator,
                tilts=tilts,
            )
            row.update(
                {
                    "series": label,
                    "family": family,
                    "outer_model": "random" if family == "random" else "structured",
                    "block_bits": constituent.block_bits,
                    "dimension": constituent.dimension,
                }
            )
            rows.append(row)

    payload = {
        "schema": "rate-half-fixed-constituent-q1-k-margin-curve-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "shell-separated occupation-one union-bound margin in bits",
        "memory_schedule": {
            "formula": "M(k)=ceil(log2(k))+2 for evaluated powers of two",
            "offset": args.memory_offset,
        },
        "distance": {
            "fraction": args.distance_numerator / args.distance_denominator,
            "cutoff_rule": "D=ceil(fraction*N)",
        },
        "probability_space": {
            "structured_outer": "one fixed authenticated constituent reused in every row",
            "random_outer": "one uniform full-rank [B,B/2] constituent sampled once and reused in every row",
            "routing": "independent uniform row-coordinate and transposed-region permutations",
            "inner": "independent RandomStepConv maps sampled once and shared by all messages",
        },
        "structured_spectrum_sources": spectrum_sources,
        "coverage": {
            "bch_derived_block_bits": [8, 32, 64, 128],
            "rm_block_bits": [8, 32, 128, 512],
            "random_block_bits": list(RANDOM_BLOCKS),
            "bch_stop": "no authenticated complete rate-half BCH-derived spectrum is available at block length 256 or above",
            "rm_stop": "the next rate-half RM member after length 512 is RM(5,11) at length 2048, beyond the requested cap",
            "random_stop": "requested block-length cap 1024",
        },
        "cases": rows,
        "limitations": [
            "These curves cover occupation one only and are not full distance-certificate margins.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite tilt grid supplies candidate witnesses and is not asserted optimal.",
            "Random curves average over the one sampled constituent; they do not certify that a particular sampled constituent realizes the average spectrum.",
        ],
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    fieldnames = [
        "series", "family", "outer_model", "block_bits", "dimension",
        "message_exponent", "message_bits", "memory_bits", "distance_cutoff",
        "relative_distance_lower", "q1_margin_bits_diagnostic", "dominant_weight",
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
