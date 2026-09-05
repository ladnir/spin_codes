#!/usr/bin/env python3
"""Three-band Renyi composition diagnostic at Q=L and distance 10.9%."""

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
from diagnose_finite_k20_band_compositions_qL import (  # noqa: E402
    DISTANCE,
    fast_inner_log,
    integer_composition,
    log_multinomial,
    softmax,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    L,
    N,
    high_precision_transfer_log_moment,
    load_spectrum_logs,
    log_choose,
    logistic,
    logit,
)


BANDS = ((23, 94), (95, 145), (146, 217))
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_three_band_renyi_qL_d109.json"


def band_log_moment(
    spectrum: np.ndarray,
    band: tuple[int, int],
    value_probability: float,
    order: float,
) -> float:
    terms = []
    for weight in range(band[0], band[1] + 1):
        log_a = float(spectrum[weight])
        if not math.isfinite(log_a):
            continue
        log_reference_shell = (
            log_choose(B, weight)
            + weight * math.log(value_probability)
            + (B - weight) * math.log1p(-value_probability)
        )
        terms.append(order * log_a + (1.0 - order) * log_reference_shell)
    return float(logsumexp(np.asarray(terms)))


def optimize_composition(envelope, spectrum: np.ndarray, counts: tuple[int, int, int]) -> dict[str, object]:
    frequencies = (np.asarray(counts, dtype=np.float64) + 0.5) / (L + 1.5)
    log_type_count = log_multinomial(counts)

    def decode(point: np.ndarray):
        probabilities = softmax(point[:3])
        values = np.asarray(
            [
                min(1.0 - 1e-12, max(1e-12, logistic(float(value))))
                for value in point[3:6]
            ]
        )
        surprisal = math.exp(min(4.0, max(-8.0, float(point[6]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[7]))))
        return probabilities, values, surprisal, order

    def objective(point: np.ndarray) -> float:
        probabilities, values, surprisal, order = decode(point)
        if np.any(probabilities <= 0.0):
            return math.inf
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(BANDS, values)
            ]
        )
        bit_probability = float(np.dot(probabilities, values))
        inner = fast_inner_log(envelope, bit_probability, surprisal)
        if not math.isfinite(inner):
            return math.inf
        log_conditioning = log_type_count + float(
            np.dot(np.asarray(counts), np.log(probabilities))
        )
        reference_bad = min(
            0.0,
            inner + DISTANCE * surprisal - B * log_conditioning,
        )
        conjugate = (order - 1.0) / order
        return (
            log_type_count
            + float(np.dot(np.asarray(counts), moments)) / order
            + conjugate * reference_bad
        )

    initial_values = (0.32, 0.5, 0.68)
    starts = []
    for order in (2.0, 8.0, 32.0):
        for surprisal in (2.1,):
            starts.append(
                np.asarray(
                    [
                        *np.log(frequencies),
                        *(logit(value) for value in initial_values),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ]
                )
            )
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"maxiter": 700, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    probabilities, values, surprisal, order = decode(best.x)
    moments = np.asarray(
        [
            band_log_moment(spectrum, band, value, order)
            for band, value in zip(BANDS, values)
        ]
    )
    bit_probability = float(np.dot(probabilities, values))
    inner = high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
    log_conditioning = log_type_count + float(
        np.dot(np.asarray(counts), np.log(probabilities))
    )
    reference_bad = min(0.0, inner + DISTANCE * surprisal - B * log_conditioning)
    conjugate = (order - 1.0) / order
    value = (
        log_type_count
        + float(np.dot(np.asarray(counts), moments)) / order
        + conjugate * reference_bad
    )
    return {
        "counts": list(counts),
        "fractions": [count / L for count in counts],
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "reference_type_probabilities": probabilities.tolist(),
        "reference_value_probabilities": values.tolist(),
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
        "log2_region_conditioning_probability": log_conditioning / math.log(2.0),
        "optimizer_success": bool(best.success),
    }


def compositions(grid: int) -> list[tuple[int, int, int]]:
    result = set()
    for first in range(grid + 1):
        for second in range(grid + 1 - first):
            fractions = np.asarray(
                [first / grid, second / grid, 1.0 - (first + second) / grid]
            )
            result.add(integer_composition(fractions))
    generator = np.random.default_rng(0x3BA240)
    return sorted(result)


def add_random_compositions(
    result: set[tuple[int, int, int]], random_per_concentration: int
) -> list[tuple[int, int, int]]:
    generator = np.random.default_rng(0x3BA240)
    for concentration in (0.2, 0.5, 1.0, 2.0, 5.0):
        for _ in range(random_per_concentration):
            result.add(integer_composition(generator.dirichlet(np.full(3, concentration))))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--grid", type=int, default=8)
    parser.add_argument("--random-per-concentration", type=int, default=20)
    args = parser.parse_args()
    spectrum = load_spectrum_logs()
    envelope = load_envelope(args.selection, DISTANCE / N)
    rows = []
    points = add_random_compositions(
        set(compositions(args.grid)), args.random_per_concentration
    )
    for index, counts in enumerate(points, start=1):
        row = optimize_composition(envelope, spectrum, counts)
        rows.append(row)
        print(
            f"composition={index}/{len(points)},margin={row['margin_bits']:.6f},counts={counts}",
            flush=True,
        )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-three-band-renyi-qL-d109-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "sampled_compositions": len(points),
        },
        "worst_sampled": worst,
        "rows": rows,
        "limitations": [
            "The composition set is sampled, not a complete cover.",
            "Optimization uses binary64 and final transfer evaluation uses 100-digit non-interval arithmetic.",
            "The multiband Renyi reduction and density-above-half transfer require written proof audits.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst_sampled": worst, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
