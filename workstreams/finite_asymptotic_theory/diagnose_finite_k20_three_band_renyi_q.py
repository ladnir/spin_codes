#!/usr/bin/env python3
"""Sample the three-band Renyi bound at arbitrary occupations Q."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import (  # noqa: E402
    DISTANCE,
    fast_inner_log,
    log_multinomial,
    softmax,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    L,
    N,
    high_precision_transfer_log_moment,
    load_spectrum_logs,
    logistic,
    logit,
)
from diagnose_finite_k20_three_band_renyi_qL import (  # noqa: E402
    BANDS,
    band_log_moment,
)


OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_three_band_renyi_q_d109.json"


def integer_composition(total: int, fractions: np.ndarray) -> tuple[int, int, int]:
    raw = total * fractions
    counts = np.floor(raw).astype(int)
    remainder = total - int(np.sum(counts))
    order = np.argsort(-(raw - counts))
    counts[order[:remainder]] += 1
    return tuple(int(value) for value in counts)


def optimize_composition(
    envelope,
    spectrum: np.ndarray,
    occupation: float,
    active_counts: tuple[float, float, float],
    *,
    cover_mode: bool = False,
) -> dict[str, object]:
    if not math.isclose(sum(active_counts), occupation, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError("active counts must sum to occupation")
    counts = (L - occupation, *active_counts)
    frequencies = (np.asarray(counts, dtype=np.float64) + 0.5) / (L + 2.0)
    log_type_count = log_multinomial(counts)

    def decode(point: np.ndarray):
        probabilities = softmax(point[:4])
        values = np.asarray(
            [
                min(1.0 - 1e-12, max(1e-12, logistic(float(value))))
                for value in point[4:7]
            ]
        )
        surprisal = math.exp(min(4.0, max(-12.0, float(point[7]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[8]))))
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray, high_precision: bool = False):
        probabilities, values, surprisal, order = decode(point)
        if np.any(probabilities <= 0.0):
            return math.inf, None
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(BANDS, values)
            ]
        )
        bit_probability = float(np.dot(probabilities[1:], values))
        inner = (
            high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
            if high_precision
            else fast_inner_log(envelope, bit_probability, surprisal)
        )
        if not math.isfinite(inner):
            return math.inf, None
        log_conditioning = log_type_count + float(
            np.dot(np.asarray(counts), np.log(probabilities))
        )
        reference_bad = min(
            0.0,
            inner + DISTANCE * surprisal - B * log_conditioning,
        )
        conjugate = (order - 1.0) / order
        value = (
            log_type_count
            + float(np.dot(np.asarray(active_counts), moments)) / order
            + conjugate * reference_bad
        )
        details = {
            "counts": list(counts),
            "active_counts": list(active_counts),
            "active_fractions": [count / occupation for count in active_counts],
            "reference_type_probabilities": probabilities.tolist(),
            "reference_value_probabilities": values.tolist(),
            "reference_bit_probability": bit_probability,
            "surprisal": surprisal,
            "holder_order": order,
            "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
            "log2_region_conditioning_probability": log_conditioning / math.log(2.0),
            "raw_reference_bad_log2_upper": (
                inner + DISTANCE * surprisal - B * log_conditioning
            )
            / math.log(2.0),
        }
        return value, details

    initial_values = (0.32, 0.5, 0.68)
    alpha = occupation / L
    starts = []
    surprisals = (
        (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha))
        if cover_mode
        else (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha), 2.1)
    )
    orders = (32.0,) if cover_mode else (2.0, 8.0, 32.0)
    for order in orders:
        for surprisal in surprisals:
            starts.append(
                (
                    order,
                    surprisal,
                    np.asarray(
                    [
                        *np.log(frequencies),
                        *(logit(value) for value in initial_values),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ]
                    ),
                )
            )
    best = None
    best_start = None
    for start_order, start_surprisal, start in starts:
        result = minimize(
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={
                "maxiter": 600 if cover_mode else 900,
                "xatol": 1e-8,
                "fatol": 1e-7,
            },
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
            best_start = (start_order, start_surprisal)
    assert best is not None
    assert best_start is not None
    value, details = evaluate(best.x, high_precision=True)
    assert details is not None
    return {
        "occupation": occupation,
        "occupation_density": occupation / L,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        "optimizer_start_order": best_start[0],
        "optimizer_start_surprisal": best_start[1],
        **details,
    }


def sampled_compositions(occupation: int, grid: int) -> list[tuple[int, int, int]]:
    result = set()
    for first in range(grid + 1):
        for second in range(grid + 1 - first):
            fractions = np.asarray(
                [first / grid, second / grid, 1.0 - (first + second) / grid]
            )
            result.add(integer_composition(occupation, fractions))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(65, 128, 256, 512, 1024, 2048, 4096, 6144, 8192, 8704, 8832),
    )
    parser.add_argument("--grid", type=int, default=2)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if any(not 65 <= occupation <= L for occupation in args.occupations):
        raise ValueError(f"occupations must lie in [65,{L}]")
    spectrum = load_spectrum_logs()
    envelope = load_envelope(args.selection, DISTANCE / N)
    rows = []
    for occupation in args.occupations:
        for active_counts in sampled_compositions(occupation, args.grid):
            row = optimize_composition(envelope, spectrum, occupation, active_counts)
            rows.append(row)
            print(
                f"Q={occupation},counts={active_counts},margin={row['margin_bits']:.6f}",
                flush=True,
            )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-three-band-renyi-q-d109-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "occupations": list(args.occupations),
            "composition_grid": args.grid,
        },
        "worst_sampled": worst,
        "rows": rows,
        "limitations": [
            "The occupation and composition sets are sampled, not complete covers.",
            "Optimization uses binary64 and final transfer evaluation uses 100-digit non-interval arithmetic.",
            "The inactive-category Renyi reduction and density-above-half transfer require written proof audits.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst_sampled": worst, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
