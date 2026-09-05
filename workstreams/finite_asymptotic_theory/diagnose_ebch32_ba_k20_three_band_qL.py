#!/usr/bin/env python3
"""Three-band all-active diagnostic for EBCH32--BA-3 at k=2^20."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from compare_ebch32_ba_k20 import expected_ba_log_spectrum  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import (  # noqa: E402
    log_multinomial,
    softmax,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_choose,
    logistic,
    logit,
    log_matrix_power_row_sum,
)


B = 256
L = 8192
N = B * L
EPOCHS = N // 128
LOWER_WEIGHT = int(os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22"))
if not 1 <= LOWER_WEIGHT < B // 2:
    raise ValueError("SPIN_EBCH_LOWER_WEIGHT must lie in [1,127]")
WINDOW = (LOWER_WEIGHT, B - LOWER_WEIGHT)
BANDS = ((LOWER_WEIGHT, 100), (101, 155), (156, B - LOWER_WEIGHT))
LOCAL_SPECTRUM = {0: 1, 8: 620, 12: 13888, 16: 36518, 20: 13888, 24: 620, 32: 1}
PARITY_FANOUT = (
    (31, 33) if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") == "1" else None
)
OUTPUT = WORKSTREAM / "ebch32_ba3_finite_B256_three_band_qL_diagnostic.json"


def logadd(values) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return -math.inf
    return float(logsumexp(np.asarray(finite, dtype=np.float64)))


def conditioned_spectrum() -> tuple[np.ndarray, float]:
    spectrum = expected_ba_log_spectrum(
        B, 32, 16, LOCAL_SPECTRUM, parity_fanout=PARITY_FANOUT
    )
    tail = logadd(list(spectrum[1 : WINDOW[0]]) + list(spectrum[WINDOW[1] + 1 :]))
    good = 1.0 - math.exp(tail)
    if good <= 0.0:
        raise ArithmeticError("conditioning lower bound vanished")
    result = np.full(B + 1, -math.inf, dtype=np.float64)
    for weight in range(WINDOW[0], WINDOW[1] + 1):
        if math.isfinite(float(spectrum[weight])):
            result[weight] = float(spectrum[weight]) - math.log(good)
    return result, good


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
        reference = (
            log_choose(B, weight)
            + weight * math.log(value_probability)
            + (B - weight) * math.log1p(-value_probability)
        )
        terms.append(order * log_a + (1.0 - order) * reference)
    return float(logsumexp(np.asarray(terms)))


def fast_inner_log(envelope, bit_probability: float, surprisal: float) -> float:
    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * bit_probability,
        surprisal=surprisal,
    )
    if np.min(transfer) < -1e-12:
        return math.inf
    return log_matrix_power_row_sum(np.maximum(transfer, 0.0), EPOCHS)


def log_type_probability(counts, probabilities: np.ndarray) -> float:
    total = 0.0
    for count, probability in zip(counts, probabilities):
        if not count:
            continue
        if probability <= 0.0:
            return -math.inf
        total += float(count) * math.log(float(probability))
    return total


def optimize(
    envelope,
    spectrum: np.ndarray,
    counts: tuple[int, int, int],
    distance: int,
) -> dict[str, object]:
    frequencies = (np.asarray(counts, dtype=np.float64) + 0.5) / (L + 1.5)
    log_type_count = log_multinomial(counts)

    def decode(point: np.ndarray):
        probabilities = softmax(point[:3])
        values = np.asarray(
            [min(1.0 - 1e-12, max(1e-12, logistic(float(x)))) for x in point[3:6]]
        )
        surprisal = math.exp(min(4.0, max(-8.0, float(point[6]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[7]))))
        return probabilities, values, surprisal, order

    def objective(point: np.ndarray) -> float:
        probabilities, values, surprisal, order = decode(point)
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
        log_conditioning = log_type_count + log_type_probability(
            counts, probabilities
        )
        reference_bad = min(
            0.0, inner + distance * surprisal - B * log_conditioning
        )
        conjugate = (order - 1.0) / order
        return (
            log_type_count
            + float(np.dot(np.asarray(counts), moments)) / order
            + conjugate * reference_bad
        )

    starts = []
    for order in (2.0, 8.0, 32.0):
        starts.append(
            np.asarray(
                [
                    *np.log(frequencies),
                    *(logit(value) for value in (0.32, 0.5, 0.68)),
                    math.log(2.1),
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
    inner = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal, epochs=EPOCHS
    )
    log_conditioning = log_type_count + log_type_probability(counts, probabilities)
    raw = inner + distance * surprisal - B * log_conditioning
    conjugate = (order - 1.0) / order
    value = (
        log_type_count
        + float(np.dot(np.asarray(counts), moments)) / order
        + conjugate * min(0.0, raw)
    )
    return {
        "counts": list(counts),
        "margin_bits": -value / math.log(2.0),
        "reference_type_probabilities": probabilities.tolist(),
        "reference_value_probabilities": values.tolist(),
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "raw_reference_bad_log2_upper": raw / math.log(2.0),
        "optimizer_success": bool(best.success),
    }


def compositions(grid: int) -> list[tuple[int, int, int]]:
    result = set()
    for first in range(grid + 1):
        for second in range(grid + 1 - first):
            fractions = np.asarray(
                [first / grid, second / grid, 1.0 - (first + second) / grid]
            )
            raw = fractions * L
            counts = np.floor(raw).astype(int)
            counts[np.argsort(-(raw - counts))[: L - int(np.sum(counts))]] += 1
            result.add(tuple(int(value) for value in counts))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=2)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    args = parser.parse_args()
    spectrum, good = conditioned_spectrum()
    results = []
    for relative_distance in (0.11, 0.109):
        distance = math.floor(relative_distance * N)
        envelope = load_envelope(args.selection, distance / N)
        rows = []
        for counts in compositions(args.grid):
            row = optimize(envelope, spectrum, counts, distance)
            rows.append(row)
            print(
                f"delta={relative_distance},counts={counts},margin={row['margin_bits']:.6f}",
                flush=True,
            )
        results.append(
            {
                "relative_distance": relative_distance,
                "distance": distance,
                "worst_sampled": min(rows, key=lambda row: float(row["margin_bits"])),
                "rows": rows,
            }
        )
    payload = {
        "schema": "ebch32-ba3-finite-b256-three-band-qL-diagnostic-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "window": list(WINDOW),
            "bands": [list(band) for band in BANDS],
            "good_event_probability_lower": good,
            "parity_fanout": list(PARITY_FANOUT) if PARITY_FANOUT else None,
            "conditioning_cost_all_rows_bits": -L * math.log2(good),
            "composition_grid": args.grid,
        },
        "results": results,
        "limitations": [
            "The composition grid is sampled, not a complete cover.",
            "The BA spectrum and optimizer use nearest binary64 arithmetic.",
            "Final transfers use 100-digit non-interval arithmetic.",
            "The three-band Renyi reduction still requires the same proof audit as the Golay calculation.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"results": [{"relative_distance": row["relative_distance"], "worst": row["worst_sampled"]} for row in results], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
