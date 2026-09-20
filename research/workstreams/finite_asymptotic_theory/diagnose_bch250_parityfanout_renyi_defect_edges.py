#!/usr/bin/env python3
"""Probe thin BCH250 fanout faces with a bandwise Renyi transfer.

The input is a pointwise expected-spectrum envelope.  The Renyi calculation
is valid for a fixed spectrum below that envelope, but no claim is made that
one sampled and reused fanout composition satisfies the envelope.  Numerical
optimization uses binary64.  The final RM2Sub transfer uses 100 decimal
digits without outward rounding.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_bch250_parityfanout_three_band_qL import (  # noqa: E402
    B,
    DISTANCE,
    EPOCHS,
    L,
    N,
    load_spectrum,
    log_choose,
    log_multinomial,
    softmax,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_matrix_power_row_sum,
)


BANDS = ((1, 16), (17, 233), (234, 249))
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l128_packing_expected_envelope.json"
)


def logistic(value: float) -> float:
    if value >= 0.0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def logit(value: float) -> float:
    return math.log(value) - math.log1p(-value)


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
    if not terms:
        return -math.inf
    return float(logsumexp(np.asarray(terms)))


def fast_inner_log(envelope, bit_probability: float, surprisal: float) -> float:
    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * bit_probability,
        surprisal=surprisal,
    )
    if np.min(transfer) < -1e-12:
        return math.inf
    return log_matrix_power_row_sum(np.maximum(transfer, 0.0), EPOCHS)


def optimize_composition(
    envelope,
    spectrum: np.ndarray,
    counts: tuple[int, int, int],
) -> dict[str, object]:
    active = [index for index, count in enumerate(counts) if count]
    active_counts = np.asarray([counts[index] for index in active], dtype=np.float64)
    frequencies = active_counts / L
    log_type_count = log_multinomial(counts)
    initial_values = np.asarray((0.064, 0.5, 0.936), dtype=np.float64)[active]

    def decode(point: np.ndarray):
        count = len(active)
        probabilities = softmax(point[:count])
        values = np.asarray(
            [
                min(1.0 - 1e-12, max(1e-12, logistic(float(value))))
                for value in point[count : 2 * count]
            ]
        )
        surprisal = math.exp(min(4.0, max(-8.0, float(point[-2]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[-1]))))
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray, high_precision: bool) -> float:
        probabilities, values, surprisal, order = decode(point)
        moments = np.asarray(
            [
                band_log_moment(spectrum, BANDS[index], value, order)
                for index, value in zip(active, values)
            ]
        )
        bit_probability = float(np.dot(probabilities, values))
        if high_precision:
            inner = high_precision_transfer_log_moment(
                envelope,
                bit_probability,
                surprisal,
                epochs=EPOCHS,
            )
        else:
            inner = fast_inner_log(envelope, bit_probability, surprisal)
        if not math.isfinite(inner):
            return math.inf
        log_conditioning = log_type_count + float(
            np.dot(active_counts, np.log(probabilities))
        )
        reference_bad = min(
            0.0,
            inner + DISTANCE * surprisal - B * log_conditioning,
        )
        conjugate = (order - 1.0) / order
        return (
            log_type_count
            + float(np.dot(active_counts, moments)) / order
            + conjugate * reference_bad
        )

    starts = []
    for order in (1.25, 2.0, 8.0, 32.0):
        for surprisal in (1.3, 2.1, 2.7):
            starts.append(
                np.concatenate(
                    (
                        np.log(np.maximum(frequencies, 1e-12)),
                        np.asarray([logit(value) for value in initial_values]),
                        [math.log(surprisal), math.log(order - 1.0)],
                    )
                )
            )
    best = None
    for start in starts:
        result = minimize(
            lambda point: evaluate(point, False),
            start,
            method="Nelder-Mead",
            options={"maxiter": 1200, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value = evaluate(best.x, True)
    probabilities, values, surprisal, order = decode(best.x)
    full_probabilities = [0.0] * len(BANDS)
    full_values = [0.0] * len(BANDS)
    for index, probability, shell_value in zip(active, probabilities, values):
        full_probabilities[index] = float(probability)
        full_values[index] = float(shell_value)
    return {
        "counts": list(counts),
        "margin_bits": -value / math.log(2.0),
        "log2_upper": value / math.log(2.0),
        "reference_type_probabilities": full_probabilities,
        "reference_value_probabilities": full_values,
        "surprisal": surprisal,
        "holder_order": order,
        "optimizer_success": bool(best.success),
    }


def probe_compositions() -> list[tuple[int, int, int]]:
    result = [(0, L, 0)]
    for defect_count in (1, 2, 4, 8, 16):
        result.append((defect_count, L - defect_count, 0))
        result.append((0, L - defect_count, defect_count))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    envelope = load_envelope(args.selection, DISTANCE / N)
    rows = []
    for counts in probe_compositions():
        row = optimize_composition(envelope, spectrum, counts)
        rows.append(row)
        print(f"counts,{counts},margin,{row['margin_bits']:.9f}", flush=True)
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-parityfanout-renyi-defect-edges-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
        },
        "worst_sampled": worst,
        "rows": rows,
        "limitations": [
            "The pointwise input is an expected-spectrum envelope, not a fixed-code guarantee.",
            "Only eleven thin-face probes are evaluated.",
            "Optimization and Renyi moments use binary64; the final RM2Sub transfer is not outward rounded.",
            "The bandwise Renyi reduction still requires a written proof audit for this probability space.",
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
