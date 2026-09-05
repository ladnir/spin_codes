#!/usr/bin/env python3
"""Finite Renyi diagnostic for fixed BCH250 with row-local fanout.

The fixed [250,125,>=38] BCH map is reused in every outer row.  Each row
uses an independent 128-layer ParityFanout wrapper.  The spectrum input is a
pointwise upper bound on the expected counting measure of one wrapped row.

The reduction is exact in form.  Optimization and arithmetic are nearest
binary64, except for the cancellation-safe final transfer evaluation.
Consequently, this program emits a diagnostic and not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_binomial_point_mass,
    log_choose,
    log_matrix_power_row_sum,
    logistic,
    logit,
)


B = 250
L = 8448
N = B * L
EPOCHS = N // 128
DISTANCE = 232_320
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l128_packing_expected_envelope.json"
)


def load_spectrum_logs(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = np.full(B + 1, -math.inf, dtype=np.float64)
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        value = row.get("log2_expected_multiplicity")
        if value is not None:
            result[weight] = float(value) * math.log(2.0)
    return result


def outer_renyi_log_moment(
    spectrum_logs: np.ndarray, value_probability: float, order: float
) -> float:
    terms = []
    for weight in range(1, B + 1):
        log_a = float(spectrum_logs[weight])
        if not math.isfinite(log_a):
            continue
        log_reference_shell = (
            log_choose(B, weight)
            + weight * math.log(value_probability)
            + (B - weight) * math.log1p(-value_probability)
        )
        terms.append(order * log_a + (1.0 - order) * log_reference_shell)
    return float(logsumexp(np.asarray(terms, dtype=np.float64)))


def reference_bad_log_upper(
    *,
    envelope,
    occupation: int,
    candidate_probability: float,
    value_probability: float,
    surprisal: float,
) -> tuple[float, dict[str, float]]:
    actual_bit_probability = candidate_probability * value_probability
    if not 0.0 < actual_bit_probability < 1.0:
        return math.inf, {}
    transfer, details = envelope.transfer(
        candidate_probability=2.0 * actual_bit_probability,
        surprisal=surprisal,
    )
    if actual_bit_probability <= 0.5:
        log_moment = log_matrix_power_row_sum(transfer, EPOCHS)
        arithmetic_method = "binary64"
    else:
        log_moment = high_precision_transfer_log_moment(
            envelope, actual_bit_probability, surprisal, epochs=EPOCHS
        )
        arithmetic_method = "mpmath-100-digit"
    log_conditioning_probability = log_binomial_point_mass(
        L, occupation, candidate_probability
    )
    raw = (
        log_moment
        + DISTANCE * surprisal
        - B * log_conditioning_probability
    )
    return min(0.0, raw), {
        "actual_bit_probability": actual_bit_probability,
        "log_inner_moment": log_moment,
        "log_region_conditioning_probability": log_conditioning_probability,
        "raw_reference_bad_log_upper": raw,
        "transfer_arithmetic_method": arithmetic_method,
        **details,
    }


def objective(
    point: np.ndarray,
    *,
    envelope,
    spectrum_logs: np.ndarray,
    occupation: int,
) -> tuple[float, dict[str, float]]:
    candidate_probability = logistic(float(point[0]))
    value_probability = logistic(float(point[1]))
    surprisal = math.exp(min(4.0, max(-16.0, float(point[2]))))
    order = 1.0 + math.exp(min(14.0, max(-14.0, float(point[3]))))
    outer_moment = outer_renyi_log_moment(
        spectrum_logs, value_probability, order
    )
    reference_bad, details = reference_bad_log_upper(
        envelope=envelope,
        occupation=occupation,
        candidate_probability=candidate_probability,
        value_probability=value_probability,
        surprisal=surprisal,
    )
    if not math.isfinite(reference_bad):
        return math.inf, {}
    conjugate_weight = (order - 1.0) / order
    value = (
        log_choose(L, occupation)
        + occupation * outer_moment / order
        + conjugate_weight * reference_bad
    )
    return value, {
        "candidate_probability": candidate_probability,
        "value_probability": value_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "outer_renyi_log_moment": outer_moment,
        "conjugate_weight": conjugate_weight,
        **details,
    }


def optimize_occupation(
    envelope, spectrum_logs: np.ndarray, occupation: int
) -> dict[str, object]:
    alpha = occupation / L
    starts = []
    for order in (1.01, 1.1, 2.0, 8.0, 32.0, 128.0):
        for value_probability in (0.4, 0.5, 0.6):
            candidate_probability = min(0.999, max(0.001, alpha))
            surprisal = max(1e-5, 1.6 * alpha * value_probability)
            starts.append(
                np.asarray(
                    [
                        logit(candidate_probability),
                        logit(value_probability),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ],
                    dtype=np.float64,
                )
            )
    best = None
    for start in starts:
        result = minimize(
            lambda point: objective(
                point,
                envelope=envelope,
                spectrum_logs=spectrum_logs,
                occupation=occupation,
            )[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1200, "xatol": 1e-8, "fatol": 1e-8},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = objective(
        best.x,
        envelope=envelope,
        spectrum_logs=spectrum_logs,
        occupation=occupation,
    )
    return {
        "occupation": occupation,
        "occupation_density": alpha,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(101, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048,
                 3072, 4096, 6144, 8192, L),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if any(not 1 <= value <= L for value in args.occupations):
        raise ValueError(f"occupations must lie in [1,{L}]")

    envelope = load_envelope(args.selection, DISTANCE / N)
    spectrum_logs = load_spectrum_logs(args.spectrum)
    rows = []
    for occupation in args.occupations:
        row = optimize_occupation(envelope, spectrum_logs, occupation)
        rows.append(row)
        print(
            f"occupation={occupation},margin_bits={row['margin_bits']:.6f},"
            f"order={row['holder_order']:.6f}",
            flush=True,
        )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-rowlocal-fanout-renyi-dense-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_dimension": 125,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "inner_epochs": EPOCHS,
            "fanout_layers": 128,
            "occupations": list(args.occupations),
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
        },
        "method": {
            "outer": (
                "pointwise expected spectrum upper for one fixed BCH250 map "
                "with independent row-local fanout wrappers"
            ),
            "change_of_measure": "Renyi/Holder against Bernoulli-y row values",
            "region_relaxation": (
                "iid Bernoulli-p candidate positions divided by the exact "
                "Binomial(8448,p) mass at Q in each of 250 regions"
            ),
            "inner": (
                "finite three-state RM2Sub transfer powered through 16500 epochs"
            ),
        },
        "worst_sampled": worst,
        "rows": rows,
        "all_sampled_below_2^-40": all(
            float(row["margin_bits"]) > 40.0 for row in rows
        ),
        "limitations": [
            "The occupation set is sampled and is not an interval cover.",
            "Optimization and most arithmetic use nearest binary64.",
            "The final transfer evaluation is cancellation-safe but is not outward rounded.",
            "The spectrum bounds maximize each output shell separately; the resulting common majorant is valid but need not be jointly attainable.",
            "The probability space resamples the fanout wrapper independently in each row; it does not reuse one fanout composition.",
        ],
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"output={args.output}")
    print(f"worst_margin_bits={worst['margin_bits']:.9f}")


if __name__ == "__main__":
    main()
