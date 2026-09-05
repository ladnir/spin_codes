#!/usr/bin/env python3
"""Probe arbitrary-Q three-band EBCH32--BA-3 compositions at 10.9%."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import log_multinomial, softmax  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    logistic,
    logit,
)
from diagnose_ebch32_ba_k20_three_band_qL import (  # noqa: E402
    B,
    BANDS,
    EPOCHS,
    L,
    N,
    band_log_moment,
    conditioned_spectrum,
    fast_inner_log,
    log_type_probability,
)


RELATIVE_DISTANCE = float(os.environ.get("SPIN_EBCH_RELATIVE_DISTANCE", "0.109"))
DISTANCE = math.floor(RELATIVE_DISTANCE * N)
OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_Q_PROBE_OUTPUT",
    "ebch32_ba3_finite_B256_three_band_q_d109_probe.json",
)
FIVE_BANDS = (
    (BANDS[0][0], 58),
    (59, 100),
    (101, 155),
    (156, 197),
    (198, BANDS[-1][1]),
)
FIVE_OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_FIVE_PROBE_OUTPUT",
    "ebch32_ba3_finite_B256_five_band_q_d109_probe.json",
)


def optimize(envelope, spectrum, active_counts: tuple[int, ...], bands) -> dict[str, object]:
    occupation = sum(active_counts)
    counts = (L - occupation, *active_counts)
    frequencies = (np.asarray(counts, dtype=np.float64) + 0.5) / (L + 2.0)
    log_type_count = log_multinomial(counts)

    def decode(point: np.ndarray):
        type_count = len(counts)
        band_count = len(active_counts)
        probabilities = softmax(point[:type_count])
        values = np.asarray(
            [
                min(1.0 - 1e-12, max(1e-12, logistic(float(x))))
                for x in point[type_count : type_count + band_count]
            ]
        )
        surprisal = math.exp(
            min(4.0, max(-12.0, float(point[type_count + band_count])))
        )
        order = 1.0 + math.exp(
            min(10.0, max(-8.0, float(point[type_count + band_count + 1])))
        )
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray, high_precision: bool = False):
        probabilities, values, surprisal, order = decode(point)
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(bands, values)
            ]
        )
        bit_probability = float(np.dot(probabilities[1:], values))
        inner = (
            high_precision_transfer_log_moment(
                envelope, bit_probability, surprisal, epochs=EPOCHS
            )
            if high_precision
            else fast_inner_log(envelope, bit_probability, surprisal)
        )
        log_conditioning = log_type_count + log_type_probability(
            counts, probabilities
        )
        raw = inner + DISTANCE * surprisal - B * log_conditioning
        conjugate = (order - 1.0) / order
        value = (
            log_type_count
            + float(np.dot(np.asarray(active_counts), moments)) / order
            + conjugate * min(0.0, raw)
        )
        return value, probabilities, values, surprisal, order, raw

    alpha = occupation / L
    initial_values = tuple(
        float(value)
        for value in np.linspace(0.32, 0.68, len(active_counts))
    )
    starts = []
    for order in (2.0, 8.0, 32.0):
        for surprisal in (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha), 2.1):
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
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 900, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, probabilities, values, surprisal, order, raw = evaluate(
        best.x, high_precision=True
    )
    return {
        "counts": list(counts),
        "active_counts": list(active_counts),
        "occupation": occupation,
        "margin_bits": -value / math.log(2.0),
        "reference_type_probabilities": probabilities.tolist(),
        "reference_value_probabilities": values.tolist(),
        "surprisal": surprisal,
        "holder_order": order,
        "raw_reference_bad_log2_upper": raw / math.log(2.0),
        "optimizer_success": bool(best.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--five-band", action="store_true")
    args = parser.parse_args()
    spectrum, good = conditioned_spectrum()
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    if args.five_band:
        bands = FIVE_BANDS
        probes = [
            tuple(3065 if index == band else 0 for index in range(5))
            for band in range(5)
        ]
        probes.extend(
            [
                (1532, 1533, 0, 0, 0),
                (2299, 766, 0, 0, 0),
                (766, 2299, 0, 0, 0),
            ]
        )
        output = FIVE_OUTPUT
    else:
        bands = BANDS
        probes = [
            (65, 0, 0),
            (0, 65, 0),
            (0, 0, 65),
            (3029, 3, 33),
            (3018, 4, 32),
            (4096, 0, 0),
            (0, 4096, 0),
            (0, 0, 4096),
        ]
        output = OUTPUT
    rows = []
    for counts in probes:
        row = optimize(envelope, spectrum, counts, bands)
        rows.append(row)
        print(f"counts={counts},margin={row['margin_bits']:.6f}", flush=True)
    payload = {
        "schema": "ebch32-ba3-finite-b256-three-band-q-d109-probe-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "distance": DISTANCE,
            "relative_distance": RELATIVE_DISTANCE,
            "good_event_probability_lower": good,
            "bands": [list(band) for band in bands],
        },
        "worst_sampled": min(rows, key=lambda row: float(row["margin_bits"])),
        "rows": rows,
        "limitations": [
            "The occupation/composition list is sampled, not a cover.",
            "Optimization uses binary64; final transfers use non-interval 100-digit arithmetic.",
        ],
    }
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst_sampled": payload["worst_sampled"], "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
