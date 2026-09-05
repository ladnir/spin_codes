#!/usr/bin/env python3
"""Run the first 10%-distance RM2Sub parameter-study tranche.

The tranche uses epoch lengths 128 and 256, message lengths 2^14 and 2^16,
and persistence exponents 20, 22, 24, and 26.  It evaluates occupation one
for exact BCH/RM constituents and matched random-outer expectations.  It also
evaluates occupation two for the exact BCH/RM constituents.  Every structured
outer is one fixed constituent repeated in all outer rows.

The output is a binary64 screening receipt.  It is not an all-occupation or
outward-rounded distance certificate.
"""

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
    Constituent,
    load_spectrum,
)
from small_k_replay.evaluate_inner_calibration import tilt_grid  # noqa: E402
from small_k_replay.evaluate_reused_random_constituent_q1 import (  # noqa: E402
    expected_spectrum,
)
from small_k_replay.evaluate_rm2sub_epoch_calibration import (  # noqa: E402
    outer_spectrum_path,
    rm2sub_case as q1_case,
)
from small_k_replay.evaluate_rm2sub_q2_calibration import (  # noqa: E402
    rm2sub_case as q2_case,
)


GENERATED = HERE / "rm2sub_calibration_constituents"
DEFAULT_JSON = HERE / "rm2sub_primary_tranche_d100.json"
DEFAULT_CSV = HERE / "rm2sub_primary_tranche_d100.csv"
STRUCTURED_KEYS = ("ebch128", "rm49")


def configuration(step_bits: int, persistence_exponent: int) -> dict[str, object]:
    state_bits = persistence_exponent - int(math.log2(step_bits))
    stem = f"t{step_bits}_s{state_bits}"
    selection_path = GENERATED / f"{stem}_selection.json"
    if not selection_path.exists():
        raise FileNotFoundError(f"missing selected RM2Sub map: {selection_path}")
    selected = json.loads(selection_path.read_text(encoding="utf-8"))["selected"]
    return {
        "name": stem,
        "step_bits": step_bits,
        "state_bits": state_bits,
        "a_distance": selected["minimum_A_distance"],
        "kernel_distance": selected["minimum_kernel_distance"],
        "kernel_weight_four": selected["weight_four_kernel_words"],
        "selection": selection_path,
        "a": GENERATED / f"{stem}_a_spectrum.json",
        "b": GENERATED / f"{stem}_b_kernel_spectrum.json",
    }


def common_fields(
    *, constituent: Constituent, exponent: int, distance: float
) -> dict[str, object]:
    message_bits = 1 << exponent
    outer_rows = message_bits // constituent.dimension
    output_bits = outer_rows * constituent.block_bits
    return {
        "constituent": constituent.name,
        "block_bits": constituent.block_bits,
        "dimension": constituent.dimension,
        "message_exponent": exponent,
        "message_bits": message_bits,
        "outer_rows": outer_rows,
        "output_bits": output_bits,
        "bad_weight": math.floor(distance * output_bits),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponents", type=int, nargs="+", default=[14, 16])
    parser.add_argument("--constituents", nargs="+", default=list(STRUCTURED_KEYS))
    parser.add_argument("--step-bits", type=int, nargs="+", default=[128, 256])
    parser.add_argument(
        "--persistence-exponents", type=int, nargs="+", default=[20, 22, 24, 26]
    )
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument("--skip-q2", action="store_true")
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    distance = args.distance_numerator / args.distance_denominator
    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)
    configs = [
        configuration(step_bits, persistence_exponent)
        for step_bits in args.step_bits
        for persistence_exponent in args.persistence_exponents
    ]
    cases: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []

    unknown = sorted(set(args.constituents) - set(STRUCTURED_KEYS))
    if unknown:
        parser.error(f"unsupported tranche constituents: {', '.join(unknown)}")

    for key in args.constituents:
        constituent = CONSTITUENTS[key]
        exact_spectrum = load_spectrum(constituent)
        exact_path = outer_spectrum_path(key, exact_spectrum)
        random_constituent = Constituent(
            name=f"random full-rank [{constituent.block_bits},{constituent.dimension}]",
            block_bits=constituent.block_bits,
            dimension=constituent.dimension,
            minimum_distance=1,
            spectrum_path=Path(),
        )
        random_spectrum = expected_spectrum(constituent.block_bits)
        random_key = f"_primary_random_{constituent.block_bits}"
        random_path = HERE / "spectra" / f"{random_key}_expected.json"
        random_path.write_text(
            json.dumps(
                {
                    "schema": "expected-binary-weight-spectrum-v1",
                    "constituent": random_constituent.name,
                    "length": random_constituent.block_bits,
                    "dimension": random_constituent.dimension,
                    "spectrum": [
                        {"weight": 0, "log2_expected_multiplicity": 0.0},
                        *[
                            {
                                "weight": weight,
                                "log2_expected_multiplicity": math.log2(count),
                            }
                            for weight, count in sorted(random_spectrum.items())
                            if weight > 0 and count > 0
                        ],
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        for exponent in args.message_exponents:
            if (1 << exponent) % constituent.dimension:
                continue
            common = common_fields(
                constituent=constituent, exponent=exponent, distance=distance
            )
            for config in configs:
                step_bits = int(config["step_bits"])
                if int(common["outer_rows"]) % step_bits:
                    skipped.append(
                        {
                            "constituent": constituent.name,
                            "message_exponent": exponent,
                            "configuration": config["name"],
                            "reason": "epoch length does not divide the transposed region",
                        }
                    )
                    continue

                print(f"case,{key},kexp,{exponent},{config['name']},Q1", flush=True)
                row = q1_case(
                    key=key,
                    spectrum_path=exact_path,
                    exponent=exponent,
                    distance=distance,
                    tilts=tilts,
                    config=config,
                )
                row.update(
                    {
                        **common,
                        "family": "structured",
                        "outer_model": "fixed-exact-spectrum",
                        "occupation": 1,
                    }
                )
                cases.append(row)

                CONSTITUENTS[random_key] = random_constituent
                try:
                    print(
                        f"case,random{constituent.block_bits},kexp,{exponent},"
                        f"{config['name']},Q1",
                        flush=True,
                    )
                    row = q1_case(
                        key=random_key,
                        spectrum_path=random_path,
                        exponent=exponent,
                        distance=distance,
                        tilts=tilts,
                        config=config,
                    )
                finally:
                    del CONSTITUENTS[random_key]
                row.update(
                    {
                        **common_fields(
                            constituent=random_constituent,
                            exponent=exponent,
                            distance=distance,
                        ),
                        "family": "random-control",
                        "outer_model": "single-reused-constituent-Q1-ensemble-average",
                        "occupation": 1,
                    }
                )
                cases.append(row)

                if not args.skip_q2:
                    print(f"case,{key},kexp,{exponent},{config['name']},Q2", flush=True)
                    row = q2_case(
                        key=key,
                        exponent=exponent,
                        config=config,
                        bad_weight=int(common["bad_weight"]),
                        tilts=tilts,
                    )
                    row.update(
                        {
                            **common,
                            "family": "structured",
                            "outer_model": "fixed-exact-spectrum",
                            "occupation": 2,
                        }
                    )
                    cases.append(row)

    payload = {
        "schema": "rm2sub-primary-parameter-tranche-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "primary_target": {
            "relative_distance": distance,
            "bad_event": (
                "output weight at most floor((distance_numerator/"
                "distance_denominator)*N)"
            ),
            "target_margin_bits": 40,
        },
        "bonus_target": (
            "10.9% or 11% may be evaluated for surviving configurations, but does "
            "not influence this tranche's ranking"
        ),
        "probability_space": {
            "structured_outer": (
                "one fixed authenticated BCH or RM constituent repeated in every row"
            ),
            "random_outer": (
                "one uniform full-rank constituent sampled once and repeated; only "
                "the joint occupation-one expectation is evaluated"
            ),
            "routing": (
                "independent uniform row-coordinate permutations and independent "
                "uniform position permutations in the transposed regions"
            ),
            "inner": (
                "one fixed exactly audited A/B pair and independent nonzero field "
                "multipliers in each epoch"
            ),
        },
        "admissible_length_rule": (
            "the constituent dimension divides k and the epoch length divides "
            "L=k/dimension; no padding or partial epoch is used"
        ),
        "parameters": {
            "message_exponents": args.message_exponents,
            "constituents": args.constituents,
            "step_bits": args.step_bits,
            "persistence_exponents": args.persistence_exponents,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "tilt_count": len(tilts),
            "q2_enabled": not args.skip_q2,
        },
        "inner_configurations": [
            {
                "name": config["name"],
                "step_bits": config["step_bits"],
                "state_bits": config["state_bits"],
                "persistence_exponent": int(config["state_bits"])
                + math.log2(int(config["step_bits"])),
                "a_minimum_distance": config["a_distance"],
                "kernel_minimum_distance": config["kernel_distance"],
                "kernel_weight_four": config["kernel_weight_four"],
                "selection_receipt": str(config["selection"]),
            }
            for config in configs
        ],
        "cases": cases,
        "skipped_inadmissible_cases": skipped,
        "limitations": [
            "The receipt covers occupations one and two only.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite tilt grid is not proved optimal.",
            "Random-control values do not certify a realized constituent.",
            "BCH plus fanout is deferred until its single-reused-constituent envelope is specified.",
            "No arbitrary-length wrapper or implementation benchmark is included.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "family",
        "outer_model",
        "constituent",
        "block_bits",
        "dimension",
        "message_exponent",
        "message_bits",
        "outer_rows",
        "output_bits",
        "bad_weight",
        "occupation",
        "configuration",
        "step_bits",
        "state_bits",
        "persistence_exponent",
        "epochs_per_region",
        "a_minimum_distance",
        "kernel_minimum_distance",
        "kernel_weight_four",
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
