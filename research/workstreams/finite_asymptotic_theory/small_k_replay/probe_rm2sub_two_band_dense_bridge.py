#!/usr/bin/env python3
"""Combine sparse and dense bounds for two-band RM2Sub compositions.

The companion sparse probe is strongest when few input bits survive.  This
probe supplies a dense reference for each fixed two-band composition.  In
each transposed region, sample a categorical mark at every position and
condition on the required counts of the two row types.  After independent
Bernoulli thinning within each type, the input bits are independent with a
known common marginal.  The dense RM2Sub epoch envelope then applies.

The evaluator takes the pointwise minimum of the sparse and dense upper
bounds before summing compositions.  Nearest-binary64 arithmetic and finite
witness grids make every reported margin diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
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
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402
from small_k_replay.probe_rm2sub_two_band_bridge import (  # noqa: E402
    LOG2,
    log_multinomial_rows,
)


DEFAULT_SPARSE = HERE / "rm2sub_two_band_bridge_all_partners_probe_d100.json"
DEFAULT_OUTPUT = HERE / "rm2sub_two_band_dense_bridge_probe_d100.json"


def log_power_moment(transfer: np.ndarray, power: int) -> float:
    """Return log(e_0^T transfer^power 1) by scaled binary powering."""
    matrix_scale = float(np.max(transfer))
    if matrix_scale == 0.0:
        return -math.inf
    matrix = transfer / matrix_scale
    matrix_log_scale = math.log(matrix_scale)
    row = np.zeros(transfer.shape[0], dtype=np.float64)
    row[0] = 1.0
    row_log_scale = 0.0
    exponent = power
    while exponent:
        if exponent & 1:
            row = row @ matrix
            scale = float(np.max(row))
            if scale == 0.0:
                return -math.inf
            row /= scale
            row_log_scale += matrix_log_scale + math.log(scale)
        exponent >>= 1
        if exponent:
            matrix = matrix @ matrix
            scale = float(np.max(matrix))
            if scale == 0.0:
                return -math.inf
            matrix /= scale
            matrix_log_scale = 2.0 * matrix_log_scale + math.log(scale)
    return row_log_scale + math.log(float(np.sum(row)))


def categorical_count_log_mass(
    total: int, first: int, second: int, first_mark: float, second_mark: float
) -> float:
    zero_mark = 1.0 - first_mark - second_mark
    if min(first_mark, second_mark, zero_mark) <= 0.0:
        return -math.inf
    return (
        log_multinomial_rows(total, first, second)
        + first * math.log(first_mark)
        + second * math.log(second_mark)
        + (total - first - second) * math.log(zero_mark)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sparse-receipt", type=Path, default=DEFAULT_SPARSE)
    parser.add_argument("--partner-band-indexes", type=int, nargs="+")
    parser.add_argument(
        "--mark-scales",
        type=float,
        nargs="+",
        default=[0.75, 1.0, 1.25, 1.5, 2.0],
    )
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[value / 4.0 for value in range(-32, 1)],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    sparse = json.loads(args.sparse_receipt.read_text(encoding="utf-8"))
    if sparse["schema"] != "rm2sub-fixed-rm-two-band-sparse-bridge-v1":
        parser.error("sparse receipt has the wrong schema")
    parameters = sparse["parameters"]
    step_bits = int(parameters["step_bits"])
    state_bits = int(parameters["state_bits"])
    outer_rows = int(parameters["outer_rows"])
    outer_bits = int(parameters["outer_block_bits"])
    output_bits = int(parameters["output_bits"])
    bad_weight = int(parameters["bad_weight"])
    minimum_occupation, maximum_occupation = map(
        int, parameters["occupation_interval"]
    )
    epochs = output_bits // step_bits

    persistence = state_bits + int(math.log2(step_bits))
    config = configuration(step_bits, persistence)
    selection = json.loads(Path(config["selection"]).read_text(encoding="utf-8"))
    selected = selection["selected"]
    generators = [int(value, 16) for value in selected["A_generator_words_hex"]]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generators,
        columns=columns,
        output_bits=step_bits,
    )
    weights = state_code_weights(generators, step_bits)
    classes, counts, association = association_matrix(weights)
    dense = DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=state_bits,
        step_bits=step_bits,
        delta=int(parameters["distance_numerator"])
        / int(parameters["distance_denominator"]),
    )

    anchor = sparse["anchor_band"]
    available = {
        int(row["partner_band_index"]): row for row in sparse["pair_rows"]
    }
    partner_indexes = (
        sorted(available)
        if args.partner_band_indexes is None
        else args.partner_band_indexes
    )
    if any(index not in available for index in partner_indexes):
        parser.error("requested partner is absent from the sparse receipt")

    pair_rows = []
    displayed_pair_values = []
    for partner_index in partner_indexes:
        sparse_pair = available[partner_index]
        partner = sparse_pair["partner_band"]
        sparse_occupations = {
            int(row["occupation"]): row
            for row in sparse_pair["occupation_rows"]
        }
        occupation_rows = []
        for occupation in range(minimum_occupation, maximum_occupation + 1):
            sparse_splits = {
                (int(row["anchor_rows"]), int(row["partner_rows"])): row
                for row in sparse_occupations[occupation]["split_rows"]
            }
            split_rows = []
            mixed_values = []
            for first_count in range(1, occupation):
                second_count = occupation - first_count
                location_log = log_multinomial_rows(
                    outer_rows, first_count, second_count
                )
                outer_log = (
                    location_log
                    + first_count * float(anchor["envelope_log2"]) * LOG2
                    + second_count * float(partner["envelope_log2"]) * LOG2
                )
                dense_best = math.inf
                dense_witness: dict[str, float | str] = {}
                for mark_scale in args.mark_scales:
                    first_mark = mark_scale * first_count / outer_rows
                    second_mark = mark_scale * second_count / outer_rows
                    count_log_mass = categorical_count_log_mass(
                        outer_rows,
                        first_count,
                        second_count,
                        first_mark,
                        second_mark,
                    )
                    if not math.isfinite(count_log_mass):
                        continue
                    live_probability = (
                        first_mark * float(anchor["probability"])
                        + second_mark * float(partner["probability"])
                    )
                    candidate_probability = 2.0 * live_probability
                    if not 0.0 < candidate_probability < 1.0:
                        continue
                    for log_surprisal in args.log_surprisals:
                        surprisal = math.exp(log_surprisal)
                        for envelope in ("three_state", "four_state"):
                            transfer = (
                                dense.transfer(
                                    candidate_probability=candidate_probability,
                                    surprisal=surprisal,
                                )[0]
                                if envelope == "three_state"
                                else dense.transfer_small_density(
                                    candidate_probability=candidate_probability,
                                    surprisal=surprisal,
                                )[0]
                            )
                            raw_event = (
                                -outer_bits * count_log_mass
                                + bad_weight * surprisal
                                + log_power_moment(transfer, epochs)
                            )
                            candidate = outer_log + min(0.0, raw_event)
                            if candidate < dense_best:
                                dense_best = candidate
                                dense_witness = {
                                    "mark_scale": mark_scale,
                                    "first_mark_probability": first_mark,
                                    "second_mark_probability": second_mark,
                                    "live_bit_probability": live_probability,
                                    "log_surprisal": log_surprisal,
                                    "envelope": envelope,
                                }
                sparse_row = sparse_splits[(first_count, second_count)]
                sparse_value = float(sparse_row["log2_upper_diagnostic"]) * LOG2
                combined = min(sparse_value, dense_best)
                selected_method = "sparse" if sparse_value <= dense_best else "dense"
                mixed_values.append(combined)
                split_rows.append(
                    {
                        "anchor_rows": first_count,
                        "partner_rows": second_count,
                        "sparse_margin_bits_diagnostic": -sparse_value / LOG2,
                        "dense_margin_bits_diagnostic": -dense_best / LOG2,
                        "selected_method": selected_method,
                        "combined_margin_bits_diagnostic": -combined / LOG2,
                        "combined_gap_to_40_bits": -combined / LOG2 - 40.0,
                        "dense_witness": dense_witness,
                    }
                )
            aggregate = float(logsumexp(np.asarray(mixed_values)))
            displayed_pair_values.append(aggregate)
            dominant = min(
                split_rows,
                key=lambda row: float(row["combined_margin_bits_diagnostic"]),
            )
            occupation_rows.append(
                {
                    "occupation": occupation,
                    "combined_mixed_log2_upper_diagnostic": aggregate / LOG2,
                    "combined_mixed_margin_bits_diagnostic": -aggregate / LOG2,
                    "combined_mixed_gap_to_40_bits": -aggregate / LOG2 - 40.0,
                    "dominant_mixed_split": {
                        "anchor_rows": dominant["anchor_rows"],
                        "partner_rows": dominant["partner_rows"],
                        "selected_method": dominant["selected_method"],
                        "combined_margin_bits_diagnostic": dominant[
                            "combined_margin_bits_diagnostic"
                        ],
                    },
                    "split_rows": split_rows,
                }
            )
        pair_aggregate = float(
            logsumexp(
                np.asarray(
                    [
                        row["combined_mixed_log2_upper_diagnostic"] * LOG2
                        for row in occupation_rows
                    ]
                )
            )
        )
        pair_rows.append(
            {
                "partner_band_index": partner_index,
                "partner_band": partner,
                "combined_interval_log2_upper_diagnostic": pair_aggregate / LOG2,
                "combined_interval_margin_bits_diagnostic": -pair_aggregate / LOG2,
                "combined_interval_gap_to_40_bits": -pair_aggregate / LOG2 - 40.0,
                "occupation_rows": occupation_rows,
            }
        )
        print(
            f"pair,{int(anchor['index'])},{partner_index},"
            f"combined_margin,{(-pair_aggregate / LOG2):.6f},"
            f"gap40,{(-pair_aggregate / LOG2 - 40.0):.6f}",
            flush=True,
        )

    displayed_aggregate = float(logsumexp(np.asarray(displayed_pair_values)))
    payload = {
        "schema": "rm2sub-fixed-rm-two-band-sparse-dense-bridge-v1",
        "status": "BINARY64_PARTIAL_DIAGNOSTIC",
        "source_sparse_receipt": str(args.sparse_receipt),
        "construction": sparse["construction"],
        "probability_space": sparse["probability_space"],
        "parameters": {
            **parameters,
            "partner_band_indexes": partner_indexes,
            "dense_mark_scales": args.mark_scales,
            "dense_log_surprisals": args.log_surprisals,
        },
        "anchor_band": anchor,
        "pair_rows": pair_rows,
        "union_of_displayed_mixed_pairs": {
            "log2_upper_diagnostic": displayed_aggregate / LOG2,
            "margin_bits_diagnostic": -displayed_aggregate / LOG2,
            "gap_to_40_bits": -displayed_aggregate / LOG2 - 40.0,
        },
        "proof_reduction": [
            *sparse["proof_reduction"],
            "For the dense alternative, sample independent categorical position marks in each transposed region and condition on the fixed two-band row counts.",
            "Under the unconditioned categorical law, Bernoulli thinning gives independent input bits with the displayed common marginal.",
            "Apply the dense RM2Sub epoch envelope and charge the exact categorical conditioning probability in every region.",
            "Take the smaller valid upper bound for each composition before summing positive terms.",
        ],
        "limitations": [
            "Only compositions supported on the anchor band and one displayed partner band are included.",
            "Compositions using three or more bands remain open.",
            "Pure-band faces are excluded from the displayed unions.",
            "The dense mark ratios are fixed by the target composition and only the displayed common scales are tested.",
            "The band Bernoulli probabilities are inherited from the sparse receipt and are not jointly optimized with the dense witness.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite witness grids are not asserted optimal.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        "displayed_pair_union_margin,"
        f"{payload['union_of_displayed_mixed_pairs']['margin_bits_diagnostic']:.6f}",
        flush=True,
    )
    print(f"output={args.output}", flush=True)


if __name__ == "__main__":
    main()
