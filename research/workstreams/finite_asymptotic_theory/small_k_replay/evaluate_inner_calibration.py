#!/usr/bin/env python3
"""Calibrate RM2Sub against RandomStepConv on exact BCH spectra.

The experiment fixes the outer constituent and the complete routing law.  It
compares occupation-one bounds for the same integer bad-weight threshold.
RM2Sub uses an authenticated t=128 constituent.  RandomStepConv uses a state
with the same number of bits.  All calculations use nearest binary64.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

import evaluate_ebch128_randomstepconv_g1 as randomstep  # noqa: E402
from analyze_riffle_splitstate_preaddmul_one_active import (  # noqa: E402
    evaluate as evaluate_rm2sub,
)
from evaluate_ebch128_randomstepconv_q1_exact import (  # noqa: E402
    uniform_coefficients,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    Constituent,
    load_spectrum,
    one_active_region_coefficients,
)


LOG2 = math.log(2.0)
DEFAULT_JSON = HERE / "inner_calibration_bch_q1_d100.json"
DEFAULT_CSV = HERE / "inner_calibration_bch_q1_d100.csv"

RM2SUB_ROOT = REPOSITORY / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts"
)
SMALLER_ROOT = REPOSITORY / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_t128_s32/receipts/smaller_state"
)

RM2SUB_CONFIGS = {
    15: {
        "a": SMALLER_ROOT / "s15_rm2sub_a_spectrum_audit.json",
        "b": SMALLER_ROOT / "s15_rm2sub_b_kernel_spectrum.json",
        "distance": 48,
    },
    16: {
        "a": SMALLER_ROOT / "s16_rm2sub_a_spectrum_audit.json",
        "b": SMALLER_ROOT / "s16_rm2sub_b_kernel_spectrum.json",
        "distance": 48,
    },
    18: {
        "a": RM2SUB_ROOT / "min_state/s18_rm2sub_a_spectrum_audit.json",
        "b": RM2SUB_ROOT / "min_state/s18_rm2sub_b_kernel_spectrum.json",
        "distance": 48,
    },
    19: {
        "a": RM2SUB_ROOT / "min_state/s19_rm2sub_a_spectrum_audit.json",
        "b": RM2SUB_ROOT / "min_state/s19_rm2sub_b_kernel_spectrum.json",
        "distance": 48,
    },
    20: {
        "a": RM2SUB_ROOT / "s20_rm2sub_a_spectrum_audit.json",
        "b": RM2SUB_ROOT / "s20_rm2sub_b_kernel_spectrum.json",
        "distance": 32,
    },
}


def tilt_grid(minimum: float, maximum: float, spacing: float) -> list[float]:
    count = math.floor((maximum - minimum) / spacing + 1e-12)
    return [float(minimum + spacing * index) for index in range(count + 1)]


def evaluate_randomstep_case(
    *,
    constituent: Constituent,
    spectrum: dict[int, int],
    message_exponent: int,
    memory_bits: int,
    bad_weight: int,
    tilts: list[float],
) -> dict[str, object]:
    message_bits = 1 << message_exponent
    outer_rows = message_bits // constituent.dimension
    weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    best = {weight: math.inf for weight in weights}
    witness = {weight: math.nan for weight in weights}
    for log_surprisal in tilts:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = randomstep.step_matrices(z, memory_bits)
        zero_region, one_region = one_active_region_coefficients(
            zero, active, outer_rows
        )
        coordinates = uniform_coefficients(
            zero_region,
            one_region,
            constituent.block_bits,
            constituent.block_bits,
        )
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        correction = bad_weight * surprisal
        for weight in weights:
            value = (
                math.log(outer_rows)
                + math.log(spectrum[weight])
                + min(0.0, float(moments[weight]) + correction)
            )
            if value < best[weight]:
                best[weight] = value
                witness[weight] = log_surprisal
    aggregate = float(logsumexp(list(best.values())))
    dominant = max(weights, key=lambda weight: best[weight])
    return {
        "inner": "RandomStepConv",
        "state_bits": memory_bits,
        "step_bits": 1,
        "epochs_per_region": outer_rows,
        "margin_bits": -aggregate / LOG2,
        "dominant_weight": dominant,
        "dominant_log_surprisal": witness[dominant],
    }


def evaluate_rm2sub_case(
    *,
    constituent: Constituent,
    rm2sub_spectrum_path: Path,
    message_exponent: int,
    state_bits: int,
    relative_distance: float,
    tilts: list[float],
) -> dict[str, object]:
    message_bits = 1 << message_exponent
    outer_rows = message_bits // constituent.dimension
    step_bits = 128
    if outer_rows % step_bits:
        raise ValueError(
            f"RM2Sub t={step_bits} does not divide region length {outer_rows}"
        )
    config = RM2SUB_CONFIGS[state_bits]
    args = Namespace(
        message_bits=message_bits,
        outer_bits=constituent.block_bits,
        outer_dimension=constituent.dimension,
        modeled_minimum_distance=constituent.minimum_distance,
        outer_spectrum=rm2sub_spectrum_path,
        random_linear_extension=False,
        random_linear_extension_bits=0,
        relative_distance=relative_distance,
        step_bits=step_bits,
        state_bits=state_bits,
        constituent_distance=config["distance"],
        live_moment_order=3,
        log_surprisals=tilts,
        activation=config["b"],
        live_spectrum=config["a"],
        ideal_random_activation=False,
        three_state_refresh=False,
    )
    payload = evaluate_rm2sub(args)
    return {
        "inner": "RM2Sub",
        "state_bits": state_bits,
        "step_bits": step_bits,
        "epochs_per_region": outer_rows // step_bits,
        "margin_bits": payload["aggregate_margin_bits"],
        "dominant_weight": payload["dominant_row"]["outer_weight"],
        "dominant_log_surprisal": payload["dominant_row"]["best_log_surprisal"],
        "a_minimum_distance": config["distance"],
        "a_spectrum": str(config["a"]),
        "b_activation": str(config["b"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponents", type=int, nargs="+", default=[13, 15])
    parser.add_argument("--constituents", nargs="+", default=["ebch32", "ebch128"])
    parser.add_argument("--state-bits", type=int, nargs="+", default=sorted(RM2SUB_CONFIGS))
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--tilt-minimum", type=float, default=-12.0)
    parser.add_argument("--tilt-maximum", type=float, default=0.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.1)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv-output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()
    if not 0 < args.distance_numerator < args.distance_denominator:
        parser.error("distance fraction must lie in (0,1)")
    unknown_states = set(args.state_bits) - set(RM2SUB_CONFIGS)
    if unknown_states:
        parser.error(f"no authenticated RM2Sub receipt for states {sorted(unknown_states)}")

    tilts = tilt_grid(args.tilt_minimum, args.tilt_maximum, args.tilt_spacing)
    distance_fraction = args.distance_numerator / args.distance_denominator
    rows = []
    for key in args.constituents:
        constituent = CONSTITUENTS[key]
        spectrum = load_spectrum(constituent)
        rm2sub_spectrum_path = HERE / "spectra" / f"{key}_weight_counts.json"
        rm2sub_spectrum_path.parent.mkdir(parents=True, exist_ok=True)
        rm2sub_spectrum_path.write_text(
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
        for exponent in args.message_exponents:
            message_bits = 1 << exponent
            if message_bits % constituent.dimension:
                continue
            outer_rows = message_bits // constituent.dimension
            output_bits = outer_rows * constituent.block_bits
            bad_weight = math.floor(distance_fraction * output_bits)
            for state_bits in args.state_bits:
                print(
                    f"case,{constituent.name},kexp,{exponent},s,{state_bits}",
                    flush=True,
                )
                common = {
                    "constituent": constituent.name,
                    "block_bits": constituent.block_bits,
                    "dimension": constituent.dimension,
                    "message_exponent": exponent,
                    "message_bits": message_bits,
                    "output_bits": output_bits,
                    "outer_rows": outer_rows,
                    "bad_weight": bad_weight,
                    "distance_fraction": distance_fraction,
                }
                rsc = evaluate_randomstep_case(
                    constituent=constituent,
                    spectrum=spectrum,
                    message_exponent=exponent,
                    memory_bits=state_bits,
                    bad_weight=bad_weight,
                    tilts=tilts,
                )
                rsc["comparison"] = "same-stored-state-bits"
                rsc["reference_rm2sub_state_bits"] = state_bits
                rsc.update(common)
                rows.append(rsc)
                persistence_memory = state_bits + 7
                rsc_persistence = evaluate_randomstep_case(
                    constituent=constituent,
                    spectrum=spectrum,
                    message_exponent=exponent,
                    memory_bits=persistence_memory,
                    bad_weight=bad_weight,
                    tilts=tilts,
                )
                rsc_persistence["comparison"] = "same-live-state-reset-rate"
                rsc_persistence["reference_rm2sub_state_bits"] = state_bits
                rsc_persistence.update(common)
                rows.append(rsc_persistence)
                rm = evaluate_rm2sub_case(
                    constituent=constituent,
                    rm2sub_spectrum_path=rm2sub_spectrum_path,
                    message_exponent=exponent,
                    state_bits=state_bits,
                    relative_distance=distance_fraction,
                    tilts=tilts,
                )
                rm["comparison"] = "rm2sub"
                rm["reference_rm2sub_state_bits"] = state_bits
                rm.update(common)
                rows.append(rm)

    payload = {
        "schema": "bch-inner-q1-calibration-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "shell-separated occupation-one union-bound margin in bits",
        "probability_space": {
            "outer": "one fixed exact-spectrum BCH constituent reused in every row",
            "routing": "independent uniform row-coordinate and transposed-region permutations",
            "randomstep": "independent random linear state transitions sampled at every bit position",
            "rm2sub": "one authenticated A,B pair with BA=0; independent nonzero state multipliers sampled per t-bit epoch",
        },
        "integer_event": "output weight at most floor(0.10*N)",
        "parameters": {
            "message_exponents": args.message_exponents,
            "constituents": args.constituents,
            "state_bits": args.state_bits,
            "rm2sub_step_bits": 128,
            "persistence_matching_rule": "RandomStepConv M=s+log2(t)=s+7",
            "tilt_count": len(tilts),
        },
        "cases": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "Only occupation one is covered.",
            "The RM2Sub t=128 sweep changes the number of epochs per region through k and B; it does not yet compare independently constructed t=64 or t=256 constituents.",
            "Equal state-bit counts do not imply equal implementation cost or equal randomness cost.",
            "The persistence-matched RandomStepConv comparison equates first-order live-state reset rates, not complete transition laws.",
            "The s=20 A constituent has minimum distance 32; the s=15,16,18,19 constituents have minimum distance 48.",
        ],
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "constituent", "block_bits", "dimension", "message_exponent",
        "message_bits", "output_bits", "outer_rows", "bad_weight", "inner",
        "state_bits", "step_bits", "epochs_per_region", "comparison",
        "reference_rm2sub_state_bits", "margin_bits",
        "dominant_weight", "dominant_log_surprisal", "a_minimum_distance",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
