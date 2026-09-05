#!/usr/bin/env python3
"""Probe pure RM spectrum bands with the finite dense RM2Sub transfer."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_rm2sub_dense_occupation import (  # noqa: E402
    DenseEnvelope,
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402
from small_k_replay.evaluate_rm2sub_fixed_rm_dense import (  # noqa: E402
    LOG2,
    log_choose,
    log_row_moment_power,
)


BANDS = (
    (1, 47),
    (48, 63),
    (64, 95),
    (96, 159),
    (160, 191),
    (192, 320),
    (321, 352),
    (353, 416),
    (417, 512),
)


def band_envelope(
    spectrum: dict[int, int], block_bits: int, lower: int, upper: int
) -> tuple[float, float, int]:
    weights = [
        weight
        for weight, count in spectrum.items()
        if count and lower <= weight <= upper
    ]

    def maximum_log_density(probability: float) -> tuple[float, int]:
        rows = [
            (
                math.log(spectrum[weight])
                - log_choose(block_bits, weight)
                - weight * math.log(probability)
                - (block_bits - weight) * math.log1p(-probability),
                weight,
            )
            for weight in weights
        ]
        return max(rows)

    result = minimize_scalar(
        lambda probability: maximum_log_density(float(probability))[0],
        bounds=(1e-8, 1.0 - 1e-8),
        method="bounded",
        options={"xatol": 1e-13},
    )
    probability = float(result.x)
    envelope, maximizing_weight = maximum_log_density(probability)
    return probability, envelope, maximizing_weight


def band_envelope_at_probability(
    spectrum: dict[int, int],
    block_bits: int,
    lower: int,
    upper: int,
    probability: float,
) -> tuple[float, int]:
    return max(
        (
            math.log(count)
            - log_choose(block_bits, weight)
            - weight * math.log(probability)
            - (block_bits - weight) * math.log1p(-probability),
            weight,
        )
        for weight, count in spectrum.items()
        if count and lower <= weight <= upper
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, default=64)
    parser.add_argument("--state-bits", type=int, default=14)
    parser.add_argument("--occupations", type=int, nargs="+", default=[29, 30, 32])
    parser.add_argument("--band-indexes", type=int, nargs="+")
    parser.add_argument("--output", type=Path, default=HERE / "rm2sub_fixed_rm_dense_band_probe_d100.json")
    args = parser.parse_args()

    step_bits = args.step_bits
    state_bits = args.state_bits
    message_bits = 1 << 16
    output_bits = 2 * message_bits
    outer_rows = 256
    outer_bits = 512
    bad_weight = output_bits // 10
    config = configuration(step_bits, state_bits + int(math.log2(step_bits)))
    selection = json.loads(Path(config["selection"]).read_text(encoding="utf-8"))
    selected = selection["selected"]
    generators = [int(value, 16) for value in selected["A_generator_words_hex"]]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generators, columns=columns, output_bits=step_bits
    )
    weights = state_code_weights(generators, step_bits)
    classes, counts, association = association_matrix(weights)
    dense = DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=state_bits,
        step_bits=step_bits,
        delta=0.10,
    )
    spectrum = load_spectrum(CONSTITUENTS["rm49"])
    selected_band_ranges = (
        list(BANDS)
        if args.band_indexes is None
        else [BANDS[index] for index in args.band_indexes]
    )
    bands = [
        (lower, upper, *band_envelope(spectrum, outer_bits, lower, upper))
        for lower, upper in selected_band_ranges
    ]
    rows = []
    epsilon = 1e-10
    epochs = output_bits // step_bits
    for occupation in args.occupations:
        alpha = occupation / outer_rows
        for lower, upper, bit_probability, envelope_log, maximizing_weight in bands:

            def objective(point: np.ndarray) -> float:
                mark_probability = float(point[0])
                surprisal = float(point[1])
                row_bit_probability = float(point[2])
                selected_envelope, _weight = band_envelope_at_probability(
                    spectrum,
                    outer_bits,
                    lower,
                    upper,
                    row_bit_probability,
                )
                active_bit_probability = mark_probability * row_bit_probability
                transfer, _details = dense.transfer(
                    candidate_probability=2.0 * active_bit_probability,
                    surprisal=surprisal,
                )
                log_binomial_mass = (
                    log_choose(outer_rows, occupation)
                    + occupation * math.log(mark_probability)
                    + (outer_rows - occupation) * math.log1p(-mark_probability)
                )
                raw = (
                    -outer_bits * log_binomial_mass
                    + bad_weight * surprisal
                    + log_row_moment_power(transfer, epochs)
                )
                return (
                    log_choose(outer_rows, occupation)
                    + occupation * selected_envelope
                    + min(0.0, raw)
                )

            starts = [
                (
                    min(1.0 - epsilon, factor * alpha),
                    max(epsilon, factor * alpha),
                    min(1.0 - epsilon, max(epsilon, row_probability)),
                )
                for factor in (1.0, 1.6)
                for row_probability in (
                    bit_probability,
                    (lower + upper) / (2.0 * outer_bits),
                    0.5,
                )
            ]
            best = None
            for start in starts:
                result = minimize(
                    objective,
                    x0=np.asarray(start),
                    bounds=(
                        (epsilon, 1.0 - epsilon),
                        (epsilon, 8.0),
                        (epsilon, 1.0 - epsilon),
                    ),
                    method="Nelder-Mead",
                    options={"maxiter": 700, "xatol": 1e-10, "fatol": 1e-9},
                )
                if best is None or float(result.fun) < float(best.fun):
                    best = result
            assert best is not None
            selected_probability = float(best.x[2])
            selected_envelope, selected_weight = band_envelope_at_probability(
                spectrum,
                outer_bits,
                lower,
                upper,
                selected_probability,
            )
            row = {
                "occupation": occupation,
                "band_lower": lower,
                "band_upper": upper,
                "minimum_envelope_bit_probability": bit_probability,
                "minimum_envelope_log2": envelope_log / LOG2,
                "selected_band_bit_probability": selected_probability,
                "selected_band_envelope_log2": selected_envelope / LOG2,
                "selected_band_envelope_maximizing_weight": selected_weight,
                "mark_probability": float(best.x[0]),
                "surprisal": float(best.x[1]),
                "margin_bits_diagnostic": -float(best.fun) / LOG2,
                "optimizer_success": bool(best.success),
            }
            rows.append(row)
            print(
                f"Q,{occupation},band,{lower}-{upper},"
                f"margin,{row['margin_bits_diagnostic']:.6f}",
                flush=True,
            )
    payload = {
        "schema": "rm2sub-fixed-rm-dense-pure-band-probe-v1",
        "status": "BINARY64_PARTIAL_DIAGNOSTIC",
        "parameters": {
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "outer_bits": outer_bits,
            "step_bits": step_bits,
            "state_bits": state_bits,
            "bad_weight": bad_weight,
            "occupations": args.occupations,
        },
        "bands": [
            {
                "lower": lower,
                "upper": upper,
                "bit_probability": probability,
                "envelope_log2": envelope / LOG2,
                "maximizing_weight": weight,
            }
            for lower, upper, probability, envelope, weight in bands
        ],
        "rows": rows,
        "limitations": [
            "Only pure band faces are evaluated; mixed band compositions remain open.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "Continuous optimizers are diagnostic and not globally certified.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
