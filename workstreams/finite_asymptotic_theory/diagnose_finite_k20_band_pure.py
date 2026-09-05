#!/usr/bin/env python3
"""Binary64 probe of pure spectrum bands in the finite B=240 dense range.

This program tests the band widths needed by the positive-mixture reduction.
It is a diagnostic, not a certificate: it covers only compositions in which
all active rows use one band reference.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize, minimize_scalar


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    DISTANCE,
    L,
    log_choose,
    load_spectrum_logs,
    logit,
    logistic,
    reference_bad_log_upper,
)


DEFAULT_BANDS = ((23, 54), (55, 94), (95, 145), (146, 185), (186, 217))
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_pure_diagnostic.json"


def band_majorant(spectrum_logs: np.ndarray, lower: int, upper: int) -> dict[str, float | int]:
    def log_mass(value_probability: float) -> float:
        return max(
            float(spectrum_logs[weight])
            - log_choose(B, weight)
            - weight * math.log(value_probability)
            - (B - weight) * math.log1p(-value_probability)
            for weight in range(lower, upper + 1)
            if math.isfinite(float(spectrum_logs[weight]))
        )

    result = minimize_scalar(
        log_mass,
        bounds=(max(1e-6, lower / B - 0.2), min(1.0 - 1e-6, upper / B + 0.2)),
        method="bounded",
        options={"xatol": 1e-13},
    )
    value_probability = float(result.x)
    return {
        "lower_weight": lower,
        "upper_weight": upper,
        "value_probability": value_probability,
        "log_majorant": log_mass(value_probability),
    }


def optimize_pure(
    envelope,
    occupation: int,
    band: dict[str, float | int],
    distance: int,
) -> dict[str, object]:
    alpha = occupation / L
    spectrum_logs = load_spectrum_logs()
    lower = int(band["lower_weight"])
    upper = int(band["upper_weight"])

    def log_majorant_at(value_probability: float) -> float:
        return max(
            float(spectrum_logs[weight])
            - log_choose(B, weight)
            - weight * math.log(value_probability)
            - (B - weight) * math.log1p(-value_probability)
            for weight in range(lower, upper + 1)
            if math.isfinite(float(spectrum_logs[weight]))
        )

    def objective(point: np.ndarray) -> tuple[float, dict[str, float]]:
        candidate_probability = logistic(float(point[0]))
        value_probability = logistic(float(point[1]))
        surprisal = math.exp(min(4.0, max(-16.0, float(point[2]))))
        log_majorant = log_majorant_at(value_probability)
        reference_bad, details = reference_bad_log_upper(
            envelope=envelope,
            occupation=occupation,
            candidate_probability=candidate_probability,
            value_probability=value_probability,
            surprisal=surprisal,
            distance=distance,
        )
        value = log_choose(L, occupation) + occupation * log_majorant + reference_bad
        return value, {
            "candidate_probability": candidate_probability,
            "value_probability": value_probability,
            "log2_majorant": log_majorant / math.log(2.0),
            "surprisal": surprisal,
            **details,
        }

    starts = []
    for scale in (0.75, 1.0, 1.25):
        candidate = min(1.0 - 1e-9, max(1e-6, alpha * scale))
        for value_probability in (
            max(1e-4, lower / B),
            (lower + upper) / (2.0 * B),
            min(1.0 - 1e-4, upper / B),
            0.5,
        ):
            for surprisal in (0.2, 0.7, 1.5, 2.5):
                starts.append(
                    np.asarray(
                        [logit(candidate), logit(value_probability), math.log(surprisal)]
                    )
                )
    best = None
    for start in starts:
        result = minimize(
            lambda point: objective(point)[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 800, "xatol": 1e-9, "fatol": 1e-8},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = objective(best.x)
    return {
        "occupation": occupation,
        "band": [int(band["lower_weight"]), int(band["upper_weight"])],
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(65, 512, 2048, 3072, 4096, 6144, 8192, 8704, 8832),
    )
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--distance", type=int, default=DISTANCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    spectrum_logs = load_spectrum_logs()
    bands = [band_majorant(spectrum_logs, *limits) for limits in DEFAULT_BANDS]
    envelope = load_envelope(args.selection, args.distance / (B * L))
    rows = []
    for occupation in args.occupations:
        for band in bands:
            row = optimize_pure(envelope, occupation, band, args.distance)
            rows.append(row)
            print(
                f"Q={occupation},band={row['band']},margin={row['margin_bits']:.6f}",
                flush=True,
            )
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-band-pure-diagnostic-v1",
        "status": "BINARY64_PURE_COMPOSITION_DIAGNOSTIC",
        "parameters": {"outer_bits": B, "outer_rows": L, "distance": args.distance},
        "bands": [
            {
                **band,
                "log2_majorant": float(band["log_majorant"]) / math.log(2.0),
            }
            for band in bands
        ],
        "rows": rows,
        "all_sampled_below_2^-40": all(float(row["margin_bits"]) > 40.0 for row in rows),
        "limitations": [
            "Only pure band compositions are tested; mixed compositions remain open.",
            "The occupation set is sampled.",
            "Optimization and arithmetic use nearest binary64.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
