#!/usr/bin/env python3
"""Exact-spectrum Q=2 calibration for RM2Sub and RandomStepConv.

Two distinct outer rows are active.  Their constituent codewords are counted
by the exact repeated-code spectrum.  Conditional on their weights, the two
row-coordinate permutations make their routed supports independent.  The
region calculation retains the without-replacement placement law.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
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
from analyze_riffle_bitshuffle_splitstate_regular_bulk import (  # noqa: E402
    load_nonzero_spectrum,
    load_uniform_nonactivation,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from evaluate_ebch128_randomstepconv_q1_exact import (  # noqa: E402
    uniform_coefficients,
)
from evaluate_single_random_constituent_q2 import (  # noqa: E402
    pair_support_coefficients,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_inner_calibration import tilt_grid  # noqa: E402
from small_k_replay.evaluate_rm2sub_epoch_calibration import (  # noqa: E402
    outer_spectrum_path,
)
from small_k_replay.evaluate_rm2sub_t64_state_calibration import (  # noqa: E402
    configuration as t64_configuration,
)


LOG2 = math.log(2.0)
DEFAULT_JSON = HERE / "rm2sub_q2_calibration_bch_k13_d100.json"
DEFAULT_CSV = HERE / "rm2sub_q2_calibration_bch_k13_d100.csv"


def aggregate_exact_spectrum(
    *, best: np.ndarray, spectrum: dict[int, int], outer_rows: int
) -> dict[str, object]:
    weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    log_pairs = math.log(math.comb(outer_rows, 2))
    terms: list[float] = []
    rows: list[tuple[float, int, int]] = []
    for first in weights:
        for second in weights:
            value = (
                log_pairs
                + math.log(spectrum[first])
                + math.log(spectrum[second])
                + float(best[first, second])
            )
            terms.append(value)
            rows.append((value, first, second))
    aggregate = float(logsumexp(terms))
    dominant, first, second = max(rows)
    return {
        "margin_bits": -aggregate / LOG2,
        "dominant_first_weight": first,
        "dominant_second_weight": second,
        "dominant_pointwise_margin_bits": -dominant / LOG2,
    }


def randomstep_case(
    *, key: str, exponent: int, memory_bits: int, bad_weight: int,
    tilts: list[float]
) -> dict[str, object]:
    constituent = CONSTITUENTS[key]
    spectrum = load_spectrum(constituent)
    outer_rows = (1 << exponent) // constituent.dimension
    best = np.full((constituent.block_bits + 1, constituent.block_bits + 1), math.inf)
    for log_surprisal in tilts:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = randomstep.step_matrices(z, memory_bits)
        region_matrices = uniform_coefficients(
            randomstep.log_entries(zero),
            randomstep.log_entries(active),
            outer_rows,
            2,
        )
        coefficients = pair_support_coefficients(
            region_matrices, constituent.block_bits
        )
        moments = np.logaddexp(coefficients[..., 0, 0], coefficients[..., 0, 1])
        best = np.minimum(best, np.minimum(0.0, moments + bad_weight * surprisal))
    result = aggregate_exact_spectrum(
        best=best, spectrum=spectrum, outer_rows=outer_rows
    )
    result.update(
        {
            "inner": "RandomStepConv",
            "configuration": f"m{memory_bits}",
            "state_bits": memory_bits,
            "step_bits": 1,
            "persistence_exponent": float(memory_bits),
            "epochs_per_region": outer_rows,
        }
    )
    return result


def rm2sub_case(
    *, key: str, exponent: int, config: dict[str, object], bad_weight: int,
    tilts: list[float]
) -> dict[str, object]:
    constituent = CONSTITUENTS[key]
    spectrum = load_spectrum(constituent)
    outer_rows = (1 << exponent) // constituent.dimension
    step_bits = int(config["step_bits"])
    state_bits = int(config["state_bits"])
    if outer_rows % step_bits:
        raise ValueError("RM2Sub epoch length must divide each transposed region")
    epochs_per_region = outer_rows // step_bits
    nonactivation = load_uniform_nonactivation(Path(config["b"]), step_bits)
    live_spectrum = load_nonzero_spectrum(
        Path(config["a"]), step_bits, state_bits
    )
    best = np.full((constituent.block_bits + 1, constituent.block_bits + 1), math.inf)
    witnesses = np.full_like(best, math.nan)
    for log_surprisal in tilts:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        impulses = splitstate_impulse_matrices(
            z=z,
            step_bits=step_bits,
            state_bits=state_bits,
            constituent_distance=int(config["a_distance"]),
            live_moment_order=3,
            live_model="factorial-moment",
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        with np.errstate(divide="ignore"):
            impulse_logs = np.log(impulses)
        region_matrices = regular_region_log_matrices(
            impulse_logs, step_bits, epochs_per_region, 2
        )
        coefficients = pair_support_coefficients(
            region_matrices, constituent.block_bits
        )
        moments = np.logaddexp(coefficients[..., 0, 0], coefficients[..., 0, 1])
        values = np.minimum(0.0, moments + bad_weight * surprisal)
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = log_surprisal
    result = aggregate_exact_spectrum(
        best=best, spectrum=spectrum, outer_rows=outer_rows
    )
    result["dominant_log_surprisal"] = float(
        witnesses[
            int(result["dominant_first_weight"]),
            int(result["dominant_second_weight"]),
        ]
    )
    result.update(
        {
            "inner": "RM2Sub",
            "configuration": config["name"],
            "state_bits": state_bits,
            "step_bits": step_bits,
            "persistence_exponent": state_bits + math.log2(step_bits),
            "epochs_per_region": epochs_per_region,
            "a_minimum_distance": config["a_distance"],
            "kernel_minimum_distance": config["kernel_distance"],
            "kernel_weight_four": config["kernel_weight_four"],
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponent", type=int, default=13)
    parser.add_argument("--constituents", nargs="+", default=["ebch32", "ebch128"])
    parser.add_argument("--state-bits", type=int, nargs="+", default=[16, 19])
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
    cases: list[dict[str, object]] = []

    for key in args.constituents:
        constituent = CONSTITUENTS[key]
        spectrum = load_spectrum(constituent)
        outer_spectrum_path(key, spectrum)
        message_bits = 1 << args.message_exponent
        outer_rows = message_bits // constituent.dimension
        output_bits = outer_rows * constituent.block_bits
        bad_weight = math.floor(distance * output_bits)
        common = {
            "constituent": constituent.name,
            "block_bits": constituent.block_bits,
            "dimension": constituent.dimension,
            "message_exponent": args.message_exponent,
            "message_bits": message_bits,
            "outer_rows": outer_rows,
            "output_bits": output_bits,
            "bad_weight": bad_weight,
        }
        for state_bits in args.state_bits:
            config = t64_configuration(state_bits)
            memory_bits = state_bits + 6
            print(f"case,{key},RSC,m{memory_bits},Q2", flush=True)
            row = randomstep_case(
                key=key,
                exponent=args.message_exponent,
                memory_bits=memory_bits,
                bad_weight=bad_weight,
                tilts=tilts,
            )
            row.update(
                {
                    "comparison": "same-live-state-reset-rate",
                    "reference_rm2sub_state_bits": state_bits,
                    **common,
                }
            )
            cases.append(row)
            print(f"case,{key},RM2Sub,t64_s{state_bits},Q2", flush=True)
            row = rm2sub_case(
                key=key,
                exponent=args.message_exponent,
                config=config,
                bad_weight=bad_weight,
                tilts=tilts,
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
        "schema": "rm2sub-q2-calibration-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "quantity": "exact-spectrum occupation-two union-bound margin in bits",
        "integer_event": (
            "output weight at most floor((distance_numerator/"
            "distance_denominator)*N)"
        ),
        "probability_space": {
            "outer": "one fixed exact-spectrum BCH constituent reused in every row",
            "active_rows": "an unordered pair of distinct outer rows",
            "routing": (
                "independent uniform coordinate permutation per row and one "
                "uniform position permutation per transposed region"
            ),
            "randomstep": "independent random linear transitions at every bit",
            "rm2sub": (
                "one fixed audited A/B pair; independent nonzero state "
                "multipliers per epoch"
            ),
        },
        "parameters": {
            "message_exponent": args.message_exponent,
            "constituents": args.constituents,
            "state_bits": args.state_bits,
            "step_bits": 64,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "persistence_matching_rule": "RandomStepConv M=s+6",
            "tilt_count": len(tilts),
        },
        "cases": cases,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "Only occupation two is covered.",
            "The finite tilt grid is not proved optimal.",
            "The selected RM2Sub A/B pairs are audited samples, not proven optimal.",
        ],
    }
    args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    fields = [
        "constituent", "block_bits", "dimension", "message_exponent",
        "message_bits", "outer_rows", "output_bits", "bad_weight", "inner",
        "configuration", "comparison", "reference_rm2sub_state_bits",
        "step_bits", "state_bits", "persistence_exponent", "epochs_per_region",
        "a_minimum_distance", "kernel_minimum_distance", "kernel_weight_four",
        "margin_bits", "dominant_first_weight", "dominant_second_weight",
        "dominant_pointwise_margin_bits",
    ]
    with args.csv_output.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cases)
    print(f"json={args.json_output}")
    print(f"csv={args.csv_output}")


if __name__ == "__main__":
    main()
