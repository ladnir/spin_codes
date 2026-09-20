#!/usr/bin/env python3
"""Compare RM2Sub epoch sizes at approximately equal live persistence.

The experiment fixes an exact-spectrum BCH outer constituent and the complete
routing law.  Each RM2Sub configuration satisfies

    state_bits + log2(step_bits) = 26.

Thus, the configurations have the same first-order live-state reset exponent.
The comparison isolates the effect of epoch length and of the authenticated
RM2Sub A/B spectra.  RandomStepConv with 26 memory bits is the common reference.
Only the occupation-one part of the bad-distance union bound is evaluated.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from argparse import Namespace
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from analyze_riffle_splitstate_preaddmul_one_active import (  # noqa: E402
    evaluate as evaluate_rm2sub,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_inner_calibration import (  # noqa: E402
    evaluate_randomstep_case,
    tilt_grid,
)


DEFAULT_JSON = HERE / "rm2sub_epoch_calibration_bch_q1_d100.json"
DEFAULT_CSV = HERE / "rm2sub_epoch_calibration_bch_q1_d100.csv"
GENERATED = HERE / "rm2sub_calibration_constituents"
T128_RECEIPTS = REPOSITORY / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state"
)

CONFIGS = (
    {
        "name": "t64_s20",
        "step_bits": 64,
        "state_bits": 20,
        "a_distance": 16,
        "kernel_distance": 8,
        "kernel_weight_four": 0,
        "a": GENERATED / "t64_s20_a_spectrum.json",
        "b": GENERATED / "t64_s20_b_kernel_spectrum.json",
    },
    {
        "name": "t128_s19",
        "step_bits": 128,
        "state_bits": 19,
        "a_distance": 48,
        "kernel_distance": 6,
        "kernel_weight_four": 0,
        "a": T128_RECEIPTS / "s19_rm2sub_a_spectrum_audit.json",
        "b": T128_RECEIPTS / "s19_rm2sub_b_kernel_spectrum.json",
    },
    {
        "name": "t256_s18",
        "step_bits": 256,
        "state_bits": 18,
        "a_distance": 96,
        "kernel_distance": 4,
        "kernel_weight_four": 576,
        "a": GENERATED / "t256_s18_a_spectrum.json",
        "b": GENERATED / "t256_s18_b_kernel_spectrum.json",
    },
)


def outer_spectrum_path(key: str, spectrum: dict[int, int]) -> Path:
    constituent = CONSTITUENTS[key]
    path = HERE / "spectra" / f"{key}_weight_counts.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(
            json.dumps(
                {
                    "schema": "exact-binary-weight-counts-v1",
                    "constituent": constituent.name,
                    "length": constituent.block_bits,
                    "dimension": constituent.dimension,
                    "weight_counts": {
                        str(weight): count for weight, count in spectrum.items()
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return path


def rm2sub_case(
    *, key: str, spectrum_path: Path, exponent: int, distance: float,
    tilts: list[float], config: dict[str, object]
) -> dict[str, object]:
    constituent = CONSTITUENTS[key]
    message_bits = 1 << exponent
    outer_rows = message_bits // constituent.dimension
    args = Namespace(
        message_bits=message_bits,
        outer_bits=constituent.block_bits,
        outer_dimension=constituent.dimension,
        modeled_minimum_distance=constituent.minimum_distance,
        outer_spectrum=spectrum_path,
        random_linear_extension=False,
        random_linear_extension_bits=0,
        relative_distance=distance,
        step_bits=config["step_bits"],
        state_bits=config["state_bits"],
        constituent_distance=config["a_distance"],
        live_moment_order=3,
        log_surprisals=tilts,
        activation=config["b"],
        live_spectrum=config["a"],
        ideal_random_activation=False,
        three_state_refresh=False,
    )
    result = evaluate_rm2sub(args)
    return {
        "inner": "RM2Sub",
        "configuration": config["name"],
        "step_bits": config["step_bits"],
        "state_bits": config["state_bits"],
        "persistence_exponent": (
            config["state_bits"] + math.log2(config["step_bits"])
        ),
        "epochs_per_region": outer_rows // config["step_bits"],
        "a_minimum_distance": config["a_distance"],
        "kernel_minimum_distance": config["kernel_distance"],
        "kernel_weight_four": config["kernel_weight_four"],
        "margin_bits": result["aggregate_margin_bits"],
        "dominant_weight": result["dominant_row"]["outer_weight"],
        "dominant_log_surprisal": result["dominant_row"][
            "best_log_surprisal"
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponents", type=int, nargs="+", default=[13, 15])
    parser.add_argument("--constituents", nargs="+", default=["ebch32", "ebch128"])
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()
    distance = args.distance_numerator / args.distance_denominator
    if not 0.0 < distance < 1.0:
        parser.error("distance fraction must lie in (0,1)")
    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)

    cases: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for key in args.constituents:
        constituent = CONSTITUENTS[key]
        spectrum = load_spectrum(constituent)
        spectrum_path = outer_spectrum_path(key, spectrum)
        for exponent in args.message_exponents:
            message_bits = 1 << exponent
            if message_bits % constituent.dimension:
                continue
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

            print(f"case,{key},kexp,{exponent},RandomStepConv,M,26", flush=True)
            reference = evaluate_randomstep_case(
                constituent=constituent,
                spectrum=spectrum,
                message_exponent=exponent,
                memory_bits=26,
                bad_weight=bad_weight,
                tilts=tilts,
            )
            reference.update(
                {
                    "configuration": "m26",
                    "persistence_exponent": 26.0,
                    **common,
                }
            )
            cases.append(reference)

            for config in CONFIGS:
                step_bits = int(config["step_bits"])
                if outer_rows % step_bits:
                    skipped.append(
                        {
                            **common,
                            "configuration": config["name"],
                            "reason": (
                                f"step_bits={step_bits} does not divide "
                                f"outer_rows={outer_rows}"
                            ),
                        }
                    )
                    continue
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
                row.update(common)
                cases.append(row)

    payload = {
        "schema": "rm2sub-epoch-calibration-v1",
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
            "randomstep": (
                "independent random linear state transitions sampled at every "
                "bit position"
            ),
            "rm2sub": (
                "one fixed, exactly audited A/B pair with BA=0 for each "
                "configuration; independent nonzero state multipliers per epoch"
            ),
        },
        "parameters": {
            "message_exponents": args.message_exponents,
            "constituents": args.constituents,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "tilt_count": len(tilts),
            "matched_persistence_exponent": 26,
        },
        "cases": cases,
        "skipped_inadmissible_cases": skipped,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "Only occupation one is covered.",
            "Equal first-order persistence does not imply equal transition laws.",
            "The t=64 and t=256 A/B pairs are selected calibration samples, not proven-optimal constituents.",
            "The comparison does not measure implementation cost or randomness cost.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "constituent", "block_bits", "dimension", "message_exponent",
        "message_bits", "outer_rows", "output_bits", "bad_weight", "inner",
        "configuration", "step_bits", "state_bits", "persistence_exponent",
        "epochs_per_region", "a_minimum_distance", "kernel_minimum_distance",
        "kernel_weight_four", "margin_bits", "dominant_weight",
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
