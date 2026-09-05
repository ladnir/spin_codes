#!/usr/bin/env python3
"""Replay the rate-half outer-family Q1 scan with calibrated RM2Sub.

The epoch size is t=64.  The state may be fixed at s=16 or scheduled so that
s+log2(t) matches the earlier RandomStepConv memory schedule M=e+2.  A case is
admissible only when the transposed region length is a multiple of 64.  No
padding or shortened epoch is introduced.  Structured constituents are fixed
and repeated.  A random control samples one full-rank constituent once and
repeats it; its Q1 average is evaluated through the exact expected spectrum.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
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
    Constituent,
    load_spectrum,
)
from small_k_replay.evaluate_inner_calibration import tilt_grid  # noqa: E402
from small_k_replay.evaluate_reused_random_constituent_q1 import (  # noqa: E402
    expected_spectrum,
)
from small_k_replay.evaluate_rm2sub_epoch_calibration import (  # noqa: E402
    rm2sub_case,
)
from small_k_replay.evaluate_rm2sub_t64_state_calibration import (  # noqa: E402
    configuration as t64_configuration,
)


DEFAULT_JSON = HERE / "rate_half_family_k_margin_d100_rm2sub_t64_s16.json"
DEFAULT_CSV = HERE / "rate_half_family_k_margin_d100_rm2sub_t64_s16.csv"
MATCHED_JSON = HERE / "rate_half_family_k_margin_d100_rm2sub_t64_matched.json"
MATCHED_CSV = HERE / "rate_half_family_k_margin_d100_rm2sub_t64_matched.csv"
SPECTRUM_DIR = HERE / "spectra" / "rm2sub_family_scan"
STRUCTURED = {
    "bch": ("ebch8", "ebch32", "xbch64", "ebch128"),
    "rm": ("rm13", "rm25", "rm37", "rm49"),
}
RANDOM_BLOCKS = (8, 16, 32, 64, 128, 256, 512, 1024)


def write_exact_spectrum(
    *, label: str, constituent: Constituent, spectrum: dict[int, int | float]
) -> Path:
    SPECTRUM_DIR.mkdir(parents=True, exist_ok=True)
    safe = label.lower().replace(" ", "_").replace("[", "").replace("]", "")
    safe = safe.replace(",", "_").replace("-", "_")
    path = SPECTRUM_DIR / f"{safe}.json"
    path.write_text(
        json.dumps(
            {
                "schema": "exact-binary-weight-counts-v1",
                "constituent": constituent.name,
                "length": constituent.block_bits,
                "dimension": constituent.dimension,
                "weight_counts": {
                    str(weight): int(count) for weight, count in spectrum.items()
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def write_expected_spectrum(
    *, label: str, constituent: Constituent, spectrum: dict[int, int | float]
) -> Path:
    SPECTRUM_DIR.mkdir(parents=True, exist_ok=True)
    path = SPECTRUM_DIR / f"random_b{constituent.block_bits}_expected.json"
    rows = [{"weight": 0, "log2_expected_multiplicity": 0.0}]
    for weight, count in sorted(spectrum.items()):
        if count > 0:
            rows.append(
                {
                    "weight": weight,
                    "log2_expected_multiplicity": math.log2(count),
                }
            )
    path.write_text(
        json.dumps(
            {
                "schema": "expected-binary-weight-spectrum-v1",
                "constituent": label,
                "length": constituent.block_bits,
                "dimension": constituent.dimension,
                "spectrum": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--message-exponents", type=int, nargs="+", default=list(range(8, 21))
    )
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument(
        "--state-schedule",
        choices=("fixed16", "persistence-matched"),
        default="fixed16",
    )
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()
    if args.state_schedule == "persistence-matched":
        if args.json_output == DEFAULT_JSON:
            args.json_output = MATCHED_JSON
        if args.csv_output == DEFAULT_CSV:
            args.csv_output = MATCHED_CSV
    distance = args.distance_numerator / args.distance_denominator
    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)
    def config_for(exponent: int) -> dict[str, object]:
        state_bits = (
            16
            if args.state_schedule == "fixed16"
            else max(7, exponent - 4)
        )
        return t64_configuration(state_bits)

    schedule_configs = {
        int(config_for(exponent)["state_bits"]): config_for(exponent)
        for exponent in args.message_exponents
    }

    families: list[
        tuple[str, str, str | None, Constituent, dict[int, int | float], Path]
    ] = []
    spectrum_sources: list[dict[str, object]] = []
    for family, keys in STRUCTURED.items():
        for key in keys:
            constituent = CONSTITUENTS[key]
            spectrum = load_spectrum(constituent)
            label = f"{constituent.name} exact"
            receipt = write_exact_spectrum(
                label=label, constituent=constituent, spectrum=spectrum
            )
            families.append((label, family, key, constituent, spectrum, receipt))
            spectrum_sources.append(
                {
                    "series": label,
                    "source_path": str(constituent.spectrum_path),
                    "source_sha256": hashlib.sha256(
                        constituent.spectrum_path.read_bytes()
                    ).hexdigest(),
                    "normalized_receipt": str(receipt),
                    "normalized_receipt_sha256": hashlib.sha256(
                        receipt.read_bytes()
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
        spectrum = expected_spectrum(block_bits)
        receipt = write_expected_spectrum(
            label=label, constituent=constituent, spectrum=spectrum
        )
        families.append((label, "random", None, constituent, spectrum, receipt))

    rows: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for label, family, constituent_key, constituent, spectrum, receipt in families:
        for exponent in args.message_exponents:
            message_bits = 1 << exponent
            if message_bits < constituent.dimension:
                skipped.append(
                    {
                        "series": label,
                        "message_exponent": exponent,
                        "reason": "message length is below constituent dimension",
                    }
                )
                continue
            if message_bits % constituent.dimension:
                skipped.append(
                    {
                        "series": label,
                        "message_exponent": exponent,
                        "reason": "constituent dimension does not divide message length",
                    }
                )
                continue
            outer_rows = message_bits // constituent.dimension
            config = config_for(exponent)
            if outer_rows % int(config["step_bits"]):
                skipped.append(
                    {
                        "series": label,
                        "message_exponent": exponent,
                        "outer_rows": outer_rows,
                        "reason": "64-bit epoch does not divide transposed region",
                    }
                )
                continue
            print(
                f"case,{label},kexp,{exponent},t64_s{config['state_bits']}",
                flush=True,
            )
            row = (
                rm2sub_case(
                    key=constituent_key,
                    spectrum_path=receipt,
                    exponent=exponent,
                    distance=distance,
                    tilts=tilts,
                    config=config,
                )
                if constituent_key is not None
                else None
            )

            # Random controls are not entries in CONSTITUENTS.  Reuse an
            # ephemeral key so the common evaluator sees the intended shape.
            if family == "random":
                ephemeral_key = f"_random_{constituent.block_bits}"
                CONSTITUENTS[ephemeral_key] = constituent
                try:
                    row = rm2sub_case(
                        key=ephemeral_key,
                        spectrum_path=receipt,
                        exponent=exponent,
                        distance=distance,
                        tilts=tilts,
                        config=config,
                    )
                finally:
                    del CONSTITUENTS[ephemeral_key]
            assert row is not None
            row.update(
                {
                    "series": label,
                    "family": family,
                    "outer_model": (
                        "random-ensemble-average" if family == "random" else "fixed"
                    ),
                    "block_bits": constituent.block_bits,
                    "dimension": constituent.dimension,
                    "message_exponent": exponent,
                    "message_bits": message_bits,
                    "outer_rows": outer_rows,
                    "output_bits": 2 * message_bits,
                    "bad_weight": math.floor(distance * 2 * message_bits),
                }
            )
            rows.append(row)

    payload = {
        "schema": "rate-half-fixed-constituent-rm2sub-q1-k-margin-curve-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "shell-separated occupation-one union-bound margin in bits",
        "distance": {
            "numerator": args.distance_numerator,
            "denominator": args.distance_denominator,
            "bad_event": "output weight at most floor((numerator/denominator)*N)",
        },
        "inner": {
            "name": "RM2Sub",
            "step_bits": 64,
            "state_schedule": args.state_schedule,
            "schedule_formula": (
                "s(e)=16"
                if args.state_schedule == "fixed16"
                else "s(e)=max(7,e-4), so s+6=e+2 for e>=11"
            ),
            "constituents": [
                {
                    "state_bits": state_bits,
                    "persistence_exponent": state_bits + 6,
                    "a_minimum_distance": selected["a_distance"],
                    "kernel_minimum_distance": selected["kernel_distance"],
                    "kernel_weight_four": selected["kernel_weight_four"],
                    "a_spectrum": str(selected["a"]),
                    "b_activation": str(selected["b"]),
                }
                for state_bits, selected in sorted(schedule_configs.items())
            ],
        },
        "probability_space": {
            "structured_outer": "one fixed authenticated constituent reused in every row",
            "random_outer": "one uniform full-rank constituent sampled once and reused in every row; Q1 averages over this sample",
            "routing": "independent uniform row-coordinate and transposed-region permutations",
            "inner": "one fixed audited A/B pair; independent nonzero state multipliers sampled per 64-bit epoch",
        },
        "admissible_length_rule": (
            "message_bits must be divisible by the constituent dimension, and "
            "outer_rows=message_bits/dimension must be divisible by 64"
        ),
        "structured_spectrum_sources": spectrum_sources,
        "cases": rows,
        "skipped_inadmissible_cases": skipped,
        "limitations": [
            "These curves cover occupation one only and are not distance certificates.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite tilt grid is not proved optimal.",
            "The random curves are joint Q1 ensemble averages and do not certify a realized constituent spectrum.",
            "No arbitrary-length padding or shortened-epoch wrapper is included.",
            "For e<11, the persistence-matched schedule is floored at s=7 because the selected RM(2,6)-subcode family contains the constant and six linear generators.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "series", "family", "outer_model", "block_bits", "dimension",
        "message_exponent", "message_bits", "outer_rows", "output_bits",
        "bad_weight", "state_bits", "step_bits", "epochs_per_region",
        "margin_bits", "dominant_weight", "dominant_log_surprisal",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
