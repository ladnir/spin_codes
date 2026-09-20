#!/usr/bin/env python3
"""Optimize the exact fixed-occupation comparison on the pure bulk face."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_bch250_parityfanout_three_band_qL import (  # noqa: E402
    B,
    DISTANCE,
    EPOCHS,
    L,
    N,
    band_majorant,
    load_spectrum,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_matrix_power_row_sum,
)


DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l128_packing_expected_envelope.json"
)
BULK = (17, 233)


def logistic(value: float) -> float:
    if value >= 0.0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def logit(value: float) -> float:
    return math.log(value) - math.log1p(-value)


def fixed_occupation_log_ratio(occupation: int, reference: float) -> tuple[float, int]:
    weights = np.arange(occupation + 1, dtype=np.float64)
    log_source = (
        gammaln(occupation + 1.0)
        - gammaln(weights + 1.0)
        - gammaln(occupation - weights + 1.0)
        - occupation * math.log(2.0)
    )
    log_target = (
        gammaln(L + 1.0)
        - gammaln(weights + 1.0)
        - gammaln(L - weights + 1.0)
        + weights * math.log(reference)
        + (L - weights) * math.log1p(-reference)
    )
    ratios = log_source - log_target
    index = int(np.argmax(ratios))
    return float(ratios[index]), index


def fast_inner_log(envelope, reference: float, surprisal: float) -> float:
    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * reference,
        surprisal=surprisal,
    )
    if np.min(transfer) < -1e-12:
        return math.inf
    return log_matrix_power_row_sum(np.maximum(transfer, 0.0), EPOCHS)


def optimize_occupation(envelope, bulk_log_majorant: float, occupation: int) -> dict[str, object]:
    log_support_count = (
        math.lgamma(L + 1)
        - math.lgamma(occupation + 1)
        - math.lgamma(L - occupation + 1)
    )

    def decode(point: np.ndarray) -> tuple[float, float]:
        reference = min(0.499999999, max(1e-10, logistic(float(point[0]))))
        surprisal = math.exp(min(3.0, max(-6.0, float(point[1]))))
        return reference, surprisal

    def objective(point: np.ndarray) -> float:
        reference, surprisal = decode(point)
        log_ratio, _ = fixed_occupation_log_ratio(occupation, reference)
        inner = fast_inner_log(envelope, reference, surprisal)
        if not math.isfinite(inner):
            return math.inf
        return (
            log_support_count
            + occupation * bulk_log_majorant
            + B * log_ratio
            + min(0.0, inner + DISTANCE * surprisal)
        )

    density = occupation / (2.0 * L)
    starts = []
    for multiplier in (0.7, 1.0, 1.4):
        reference = min(0.49, max(1e-6, density * multiplier))
        for surprisal in (1.0, 2.0, 3.0, 5.0):
            starts.append(np.asarray([logit(reference), math.log(surprisal)]))
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"maxiter": 800, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    reference, surprisal = decode(best.x)
    log_ratio, maximizing_weight = fixed_occupation_log_ratio(occupation, reference)
    inner = high_precision_transfer_log_moment(
        envelope, reference, surprisal, epochs=EPOCHS
    )
    value = (
        log_support_count
        + occupation * bulk_log_majorant
        + B * log_ratio
        + min(0.0, inner + DISTANCE * surprisal)
    )
    return {
        "occupation": occupation,
        "margin_bits": -value / math.log(2.0),
        "log2_upper": value / math.log(2.0),
        "reference_bit_probability": reference,
        "surprisal": surprisal,
        "fixed_occupation_log2_ratio_per_region": log_ratio / math.log(2.0),
        "ratio_maximizing_weight": maximizing_weight,
        "optimizer_success": bool(best.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    bulk = band_majorant(spectrum, *BULK)
    envelope = load_envelope(args.selection, DISTANCE / N)
    occupations = (
        101, 128, 192, 256, 384, 512, 768, 1024,
        1536, 2048, 3072, 4096, 6144, 8192, L,
    )
    rows = []
    for occupation in occupations:
        row = optimize_occupation(
            envelope, float(bulk["log_majorant"]), occupation
        )
        rows.append(row)
        print(
            f"occupation,{occupation},margin,{row['margin_bits']:.9f},"
            f"reference,{row['reference_bit_probability']:.12f}",
            flush=True,
        )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-rowlocal-fanout-pure-occupation-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "bulk_band": list(BULK),
            "bulk_log2_majorant": float(bulk["log_majorant"]) / math.log(2.0),
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
        },
        "worst_sampled": worst,
        "rows": rows,
        "limitations": [
            "Only the pure bulk face and fifteen sampled occupations are evaluated.",
            "The spectrum envelope is an expectation for independently sampled row-local fanout wrappers.",
            "Optimization and coefficient ratios use nearest binary64 arithmetic.",
            "The final RM2Sub transfer uses 100 decimal digits but is not outward rounded.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")
    print(f"worst_margin_bits,{worst['margin_bits']:.9f}")


if __name__ == "__main__":
    main()
