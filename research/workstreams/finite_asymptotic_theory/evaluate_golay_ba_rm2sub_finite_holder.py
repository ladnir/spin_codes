#!/usr/bin/env python3
"""Finite full-spectrum Holder diagnostic for independent conditioned BA rows.

The regular-bulk receipt uses a Bernoulli-product envelope for the outer
counting measure.  At candidate probability 1/2, its stored inner term is

    Q * log(envelope mass) + min(0, log P_U[bad]),

where U is the uniform distribution on F_2^B.  This program removes that
envelope and applies Holder directly to an upper bound on the conditioned BA
counting measure.  The calculation is binary64 and is not an outward-rounded
certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


LOG2 = math.log(2.0)
WORKSTREAM = Path(__file__).resolve().parent


def log_choose(total: int, count: int) -> float:
    if count < 0 or count > total:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(count + 1)
        - math.lgamma(total - count + 1)
    )


def holder_grid() -> np.ndarray:
    """Return a grid dense near p=1 and extending to large p."""
    return 1.0 + np.exp2(np.linspace(-20.0, 20.0, 1281))


def load_spectrum(path: Path, outer_bits: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = np.full(outer_bits + 1, -math.inf, dtype=np.float64)
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        value = row["log2_expected_multiplicity"]
        if value is not None:
            result[weight] = float(value) * LOG2
    result[0] = -math.inf
    return result


def envelope_at_half(payload: dict[str, object]) -> float:
    entries = payload["envelope"]["bernoulli_envelopes"]
    for row in entries:
        if abs(float(row["probability"]) - 0.5) < 1e-15:
            return float(row["log2_mass"]) * LOG2
    raise ValueError("the bulk receipt has no probability-1/2 envelope")


def evaluate(
    *,
    bulk_path: Path,
    spectrum_path: Path,
    one_active_path: Path,
) -> dict[str, object]:
    bulk = json.loads(bulk_path.read_text(encoding="utf-8"))
    one_active = json.loads(one_active_path.read_text(encoding="utf-8"))
    outer_bits = int(bulk["parameters"]["outer_bits"])
    outer_blocks = int(bulk["parameters"]["outer_blocks"])
    spectrum = load_spectrum(spectrum_path, outer_bits)
    old_envelope = envelope_at_half(bulk)

    # U is uniform on F_2^B.  For a support S of weight w, the conditioned
    # expected BA counting measure is at most A_w / C(B,w).  Thus its density
    # M=d mu/dU has logarithm log A_w-log C(B,w)+B log 2.
    log_reference = []
    log_counting_density = []
    density_weights = []
    for weight in range(1, outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(outer_bits, weight)
        log_reference.append(shell - outer_bits * LOG2)
        log_counting_density.append(multiplicity - shell + outer_bits * LOG2)
        density_weights.append(weight)
    log_reference_array = np.asarray(log_reference, dtype=np.float64)
    log_density_array = np.asarray(log_counting_density, dtype=np.float64)

    p_values = holder_grid()
    log_moments = np.asarray(
        [
            float(
                logsumexp(
                    log_reference_array + float(p) * log_density_array
                )
            )
            for p in p_values
        ],
        dtype=np.float64,
    )
    max_index = int(np.argmax(log_density_array))
    max_log_density = float(log_density_array[max_index])

    rows = []
    contribution_logs = []
    for source in bulk["occupation_rows"]:
        occupation = int(source["active_regular_outer_blocks"])
        probability = float(source["best_candidate_probability"])
        if abs(probability - 0.5) >= 1e-15:
            raise ValueError(
                "raw-inner recovery requires probability 1/2 for every row"
            )
        stored_inner = float(source["inner_log2_upper"]) * LOG2
        raw_inner = stored_inner - occupation * old_envelope
        # The bulk evaluator clips the raw Chernoff exponent at zero before it
        # adds the envelope.  Allow a small floating-point tolerance only.
        if raw_inner > 1e-7:
            raise ArithmeticError("recovered inner log probability is positive")
        raw_inner = min(0.0, raw_inner)

        corrections = (
            occupation * log_moments / p_values
            + (1.0 - 1.0 / p_values) * raw_inner
        )
        best_index = int(np.argmin(corrections))
        best_correction = float(corrections[best_index])
        best_p: float | str = float(p_values[best_index])

        infinity_correction = occupation * max_log_density + raw_inner
        if infinity_correction < best_correction:
            best_correction = infinity_correction
            best_p = "infinity"

        location_log = log_choose(outer_blocks, occupation)
        contribution = location_log + best_correction
        contribution_logs.append(contribution)
        rows.append(
            {
                "active_outer_blocks": occupation,
                "reference_inner_log2_upper": raw_inner / LOG2,
                "best_holder_p": best_p,
                "holder_counting_measure_log2_upper": (
                    best_correction / LOG2
                ),
                "outer_location_log2": location_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "pointwise_margin_bits": -contribution / LOG2,
            }
        )

    multi_log = float(logsumexp(np.asarray(contribution_logs)))

    selection = one_active["outer_selection"]
    good_probability_lower = float(selection["good_event_probability_lower"])
    lower, upper = map(int, selection["asymptotic_weight_window"])
    central_terms = []
    for row in one_active["one_active"]["weight_rows"]:
        weight = int(row["outer_weight"])
        if lower <= weight <= upper:
            central_terms.append(float(row["pointwise_log2_upper"]) * LOG2)
    conditional_one_log = (
        float(logsumexp(np.asarray(central_terms)))
        - math.log(good_probability_lower)
    )
    combined_log = float(np.logaddexp(conditional_one_log, multi_log))
    dominant = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]

    return {
        "schema": "golay-ba3-rm2sub-finite-holder-v1",
        "status": "BINARY64_FINITE_DIAGNOSTIC",
        "inputs": {
            "bulk_uniform_reference": str(bulk_path),
            "conditioned_spectrum_upper": str(spectrum_path),
            "one_active": str(one_active_path),
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_blocks": outer_blocks,
            "active_block_values": [
                int(row["active_outer_blocks"]) for row in rows
            ],
            "holder_grid_size": int(len(p_values)),
            "holder_p_min": float(p_values[0]),
            "holder_p_max": float(p_values[-1]),
        },
        "change_of_measure": {
            "reference": "uniform distribution on F_2^B",
            "measure": (
                "upper bound on the expected nonzero counting measure of one "
                "BA draw conditioned on the recorded tail-free event"
            ),
            "log2_total_upper_measure": (
                float(logsumexp(log_reference_array + log_density_array))
                / LOG2
            ),
            "log2_max_counting_density": max_log_density / LOG2,
            "maximizing_weight": density_weights[max_index],
            "removed_old_envelope_log2_per_row": old_envelope / LOG2,
        },
        "conditional_one_active_log2_upper": conditional_one_log / LOG2,
        "conditional_one_active_margin_bits": -conditional_one_log / LOG2,
        "sampled_multi_active_log2_upper": multi_log / LOG2,
        "sampled_multi_active_margin_bits": -multi_log / LOG2,
        "sampled_combined_log2_upper": combined_log / LOG2,
        "sampled_combined_margin_bits": -combined_log / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": [
            "The outer BA draws are independent across rows and are independently conditioned on the tail-free event.",
            "The listed multi-active occupations are sampled diagnostics, not a union over every occupation.",
            "The inner receipt and this Holder calculation use nearest-binary64 arithmetic, not outward rounding.",
            "The conditioned spectrum is a pointwise expectation upper bound obtained from Markov conditioning.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--conditioned-spectrum", type=Path, required=True)
    parser.add_argument("--one-active", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            WORKSTREAM
            / "golay_ba3_rm2sub_finite_holder_corrected_q768_9728_k20_d11.json"
        ),
    )
    args = parser.parse_args()
    result = evaluate(
        bulk_path=args.bulk,
        spectrum_path=args.conditioned_spectrum,
        one_active_path=args.one_active,
    )
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "one_active_margin_bits": result[
                    "conditional_one_active_margin_bits"
                ],
                "sampled_multi_active_margin_bits": result[
                    "sampled_multi_active_margin_bits"
                ],
                "sampled_combined_margin_bits": result[
                    "sampled_combined_margin_bits"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
