#!/usr/bin/env python3
"""Diagnose a finite dense bridge for one fixed RM outer and RM2Sub."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


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


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "rm2sub_fixed_rm_dense_d100.json"


def log_choose(total: int, selected: int) -> float:
    if not 0 <= selected <= total:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def holder_grid() -> np.ndarray:
    return 1.0 + np.exp2(np.linspace(-20.0, 20.0, 1281))


def density_statistics(
    spectrum: dict[int, int], block_bits: int, orders: np.ndarray
) -> tuple[np.ndarray, float, int]:
    """Moments of the full nonzero RM counting density relative to uniform."""
    reference_logs = []
    density_logs = []
    weights = []
    for weight, count in spectrum.items():
        if weight == 0 or count == 0:
            continue
        shell = log_choose(block_bits, weight)
        reference_logs.append(shell - block_bits * LOG2)
        density_logs.append(math.log(count) - shell + block_bits * LOG2)
        weights.append(weight)
    reference = np.asarray(reference_logs)
    density = np.asarray(density_logs)
    moments = np.asarray(
        [float(logsumexp(reference + order * density)) for order in orders]
    )
    index = int(np.argmax(density))
    return moments, float(density[index]), weights[index]


def log_row_moment_power(transfer: np.ndarray, power: int) -> float:
    """Return log(e_0^T transfer^power 1) with positive rescaling."""
    row = np.zeros(transfer.shape[0], dtype=np.float64)
    row[0] = 1.0
    log_scale = 0.0
    for _ in range(power):
        row = row @ transfer
        scale = float(np.max(row))
        if scale == 0.0:
            return -math.inf
        row /= scale
        log_scale += math.log(scale)
    return log_scale + math.log(float(np.sum(row)))


def holder_correction(
    *,
    occupation: int,
    raw_event_log: float,
    orders: np.ndarray,
    density_log_moments: np.ndarray,
    maximum_density_log: float,
) -> tuple[float, float | str]:
    clipped = min(0.0, raw_event_log)
    values = (
        occupation * density_log_moments / orders
        + (1.0 - 1.0 / orders) * clipped
    )
    index = int(np.argmin(values))
    best = float(values[index])
    witness: float | str = float(orders[index])
    infinity = occupation * maximum_density_log + clipped
    if infinity < best:
        best = infinity
        witness = "infinity"
    return best, witness


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, default=64)
    parser.add_argument("--state-bits", type=int, default=14)
    parser.add_argument("--message-exponent", type=int, default=16)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=[29, 30, 32, 40, 48, 64, 96, 128, 160, 192, 224, 256],
    )
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    constituent = CONSTITUENTS["rm49"]
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("outer dimension does not divide the message length")
    outer_rows = message_bits // constituent.dimension
    output_bits = 2 * message_bits
    if output_bits % args.step_bits:
        parser.error("epoch length does not divide the output length")
    if any(not 1 <= value <= outer_rows for value in args.occupations):
        parser.error("occupations must lie in [1,L]")
    expected_persistence = args.state_bits + int(math.log2(args.step_bits))
    config = configuration(args.step_bits, expected_persistence)
    if int(config["state_bits"]) != args.state_bits:
        raise AssertionError("configuration lookup returned the wrong state size")

    selection = json.loads(Path(config["selection"]).read_text(encoding="utf-8"))
    selected = selection["selected"]
    generators = [int(value, 16) for value in selected["A_generator_words_hex"]]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generators,
        columns=columns,
        output_bits=args.step_bits,
    )
    weights = state_code_weights(generators, args.step_bits)
    classes, counts, association = association_matrix(weights)
    dense = DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=args.state_bits,
        step_bits=args.step_bits,
        delta=args.distance_numerator / args.distance_denominator,
    )

    orders = holder_grid()
    spectrum = load_spectrum(constituent)
    density_log_moments, maximum_density_log, maximum_density_weight = (
        density_statistics(spectrum, constituent.block_bits, orders)
    )
    bad_weight = math.floor(
        args.distance_numerator * output_bits / args.distance_denominator
    )
    epochs = output_bits // args.step_bits
    epsilon = 1e-10

    rows = []
    for occupation in args.occupations:
        alpha = occupation / outer_rows

        def objective(point: np.ndarray, envelope: str) -> float:
            probability = float(point[0])
            surprisal = float(point[1])
            if envelope == "three_state":
                transfer, _details = dense.transfer(
                    candidate_probability=probability,
                    surprisal=surprisal,
                )
            else:
                transfer, _details = dense.transfer_small_density(
                    candidate_probability=probability,
                    surprisal=surprisal,
                )
            log_binomial_mass = (
                log_choose(outer_rows, occupation)
                + occupation * math.log(probability)
                + (outer_rows - occupation) * math.log1p(-probability)
            )
            raw_event = (
                -constituent.block_bits * log_binomial_mass
                + bad_weight * surprisal
                + log_row_moment_power(transfer, epochs)
            )
            correction, _holder = holder_correction(
                occupation=occupation,
                raw_event_log=raw_event,
                orders=orders,
                density_log_moments=density_log_moments,
                maximum_density_log=maximum_density_log,
            )
            return log_choose(outer_rows, occupation) + correction

        starts = [
            (min(1.0 - epsilon, max(epsilon, factor * alpha)), max(epsilon, factor * alpha))
            for factor in (1.0, 1.3, 1.6, 2.0)
        ]
        starts.extend(
            [
                (min(1.0 - epsilon, max(epsilon, alpha)), value)
                for value in (0.05, 0.2, 0.8, 1.5)
            ]
        )
        best = None
        best_envelope = ""
        for envelope in ("three_state", "four_state"):
            for start in starts:
                result = minimize(
                    lambda point: objective(point, envelope),
                    x0=np.asarray(start),
                    bounds=((epsilon, 1.0 - epsilon), (epsilon, 8.0)),
                    method="Nelder-Mead",
                    options={"maxiter": 700, "xatol": 1e-10, "fatol": 1e-9},
                )
                if best is None or float(result.fun) < float(best.fun):
                    best = result
                    best_envelope = envelope
        assert best is not None
        probability = float(best.x[0])
        surprisal = float(best.x[1])
        transfer = (
            dense.transfer(
                candidate_probability=probability, surprisal=surprisal
            )[0]
            if best_envelope == "three_state"
            else dense.transfer_small_density(
                candidate_probability=probability, surprisal=surprisal
            )[0]
        )
        log_binomial_mass = (
            log_choose(outer_rows, occupation)
            + occupation * math.log(probability)
            + (outer_rows - occupation) * math.log1p(-probability)
        )
        raw_event = (
            -constituent.block_bits * log_binomial_mass
            + bad_weight * surprisal
            + log_row_moment_power(transfer, epochs)
        )
        correction, holder_order = holder_correction(
            occupation=occupation,
            raw_event_log=raw_event,
            orders=orders,
            density_log_moments=density_log_moments,
            maximum_density_log=maximum_density_log,
        )
        total = log_choose(outer_rows, occupation) + correction
        row = {
            "occupation": occupation,
            "alpha": alpha,
            "margin_bits_diagnostic": -total / LOG2,
            "log2_upper_diagnostic": total / LOG2,
            "envelope": best_envelope,
            "candidate_probability": probability,
            "surprisal": surprisal,
            "holder_order": holder_order,
            "uniform_reference_event_log2_upper": min(0.0, raw_event) / LOG2,
            "outer_location_log2": log_choose(outer_rows, occupation) / LOG2,
            "optimizer_success": bool(best.success),
        }
        rows.append(row)
        print(
            f"Q,{occupation},margin,{row['margin_bits_diagnostic']:.6f},"
            f"envelope,{best_envelope},holder,{holder_order}",
            flush=True,
        )

    total_log = float(
        logsumexp(np.asarray([row["log2_upper_diagnostic"] * LOG2 for row in rows]))
    )
    payload = {
        "schema": "rm2sub-fixed-rm-dense-holder-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": (
            "one fixed RM(4,9) constituent repeated in every row, uniform "
            "routing, and one fixed audited RM2Sub A/B pair"
        ),
        "probability_space": {
            "outer": "one fixed authenticated RM(4,9) constituent",
            "routing": "independent uniform row-coordinate and region permutations",
            "inner": "independent nonzero field scalar in every RM2Sub epoch",
        },
        "parameters": {
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "outer_block_bits": constituent.block_bits,
            "outer_dimension": constituent.dimension,
            "bad_weight": bad_weight,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "persistence_exponent": expected_persistence,
            "epochs": epochs,
            "occupations": args.occupations,
        },
        "exact_inner_association": {
            "weight_classes": [int(value) for value in classes],
            "class_counts": [int(value) for value in counts],
            "association_matrix": [
                [int(value) for value in row] for row in association
            ],
        },
        "holder_change_of_measure": {
            "reference": "uniform measure on F_2^512 for every active row",
            "density": "complete nonzero RM(4,9) counting measure after a uniform coordinate permutation",
            "order_grid_size": int(len(orders)),
            "maximum_log2_density": maximum_density_log / LOG2,
            "maximum_density_weight": maximum_density_weight,
        },
        "occupation_rows": rows,
        "partial_margin_bits": -total_log / LOG2,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "Continuous optimization is diagnostic and not proved globally optimal.",
            "Only the listed occupations are covered.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"partial_margin={payload['partial_margin_bits']:.6f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
