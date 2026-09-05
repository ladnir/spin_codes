#!/usr/bin/env python3
"""Sample finite five-band compositions at Q=L and relative distance 10.9%."""

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
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    EPOCHS,
    L,
    N,
    high_precision_transfer_log_moment,
    log_matrix_power_row_sum,
)


PURE_RECEIPT = (
    WORKSTREAM
    / "golay_ba3_rm2sub_finite_B240_band_pure_w23_217_d109_qL_diagnostic.json"
)
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_compositions_qL_d109_diagnostic.json"
DISTANCE = 231_045


def log_multinomial(counts: tuple[int, ...]) -> float:
    return math.lgamma(sum(counts) + 1.0) - sum(
        math.lgamma(count + 1.0) for count in counts
    )


def softmax(point: np.ndarray) -> np.ndarray:
    shifted = point - np.max(point)
    values = np.exp(shifted)
    return values / np.sum(values)


def load_band_witnesses() -> list[dict[str, object]]:
    payload = json.loads(PURE_RECEIPT.read_text(encoding="utf-8"))
    rows = payload["rows"]
    if any(int(row["occupation"]) != L for row in rows):
        raise ValueError("pure receipt must contain only Q=L rows")
    return [
        {
            "band": [int(value) for value in row["band"]],
            "value_probability": float(row["value_probability"]),
            "log_majorant": float(row["log2_majorant"]) * math.log(2.0),
        }
        for row in rows
    ]


def integer_composition(fractions: np.ndarray) -> tuple[int, ...]:
    raw = fractions * L
    counts = np.floor(raw).astype(int)
    remainder = L - int(np.sum(counts))
    order = np.argsort(-(raw - counts))
    counts[order[:remainder]] += 1
    return tuple(int(value) for value in counts)


def sampled_compositions(
    band_count: int, random_per_concentration: int = 12
) -> list[tuple[int, ...]]:
    result: set[tuple[int, ...]] = set()
    for index in range(band_count):
        fractions = np.zeros(band_count)
        fractions[index] = 1.0
        result.add(integer_composition(fractions))
    for left in range(band_count):
        for right in range(left + 1, band_count):
            for share in (0.1, 0.25, 0.5, 0.75, 0.9):
                fractions = np.zeros(band_count)
                fractions[left] = share
                fractions[right] = 1.0 - share
                result.add(integer_composition(fractions))
    result.add(integer_composition(np.full(band_count, 1.0 / band_count)))
    generator = np.random.default_rng(0xBA240)
    for concentration in (0.2, 0.5, 1.0, 2.0, 5.0):
        for _ in range(random_per_concentration):
            result.add(
                integer_composition(
                    generator.dirichlet(np.full(band_count, concentration))
                )
            )
    return sorted(result)


def fast_inner_log(envelope, bit_probability: float, surprisal: float) -> float:
    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * bit_probability,
        surprisal=surprisal,
    )
    if np.min(transfer) < -1e-12:
        return math.inf
    transfer = np.maximum(transfer, 0.0)
    value = log_matrix_power_row_sum(transfer, EPOCHS)
    return value if math.isfinite(value) else math.inf


def optimize_composition(
    envelope,
    bands: list[dict[str, object]],
    counts: tuple[int, ...],
) -> dict[str, object]:
    active = [index for index, count in enumerate(counts) if count]
    active_counts = tuple(counts[index] for index in active)
    frequencies = np.asarray(active_counts, dtype=np.float64) / L
    values = np.asarray(
        [float(bands[index]["value_probability"]) for index in active]
    )
    outer_log = sum(
        count * float(bands[index]["log_majorant"])
        for index, count in enumerate(counts)
    )
    log_type_count = log_multinomial(counts)

    def objective(point: np.ndarray) -> float:
        probabilities = softmax(point[:-1])
        surprisal = math.exp(min(4.0, max(-8.0, float(point[-1]))))
        bit_probability = float(np.dot(probabilities, values))
        inner_log = fast_inner_log(envelope, bit_probability, surprisal)
        if not math.isfinite(inner_log):
            return math.inf
        log_conditioning = log_type_count + sum(
            count * math.log(probability)
            for count, probability in zip(active_counts, probabilities)
        )
        reference_bad = min(
            0.0,
            inner_log + DISTANCE * surprisal - B * log_conditioning,
        )
        return log_type_count + outer_log + reference_bad

    start_logits = np.log(np.maximum(frequencies, 1e-12))
    starts = [
        np.concatenate((start_logits, [math.log(surprisal)]))
        for surprisal in (0.7, 1.3, 2.0, 2.7)
    ]
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"maxiter": 1000, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    probabilities = softmax(best.x[:-1])
    surprisal = math.exp(min(4.0, max(-8.0, float(best.x[-1]))))
    bit_probability = float(np.dot(probabilities, values))
    inner_log = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal
    )
    log_conditioning = log_type_count + sum(
        count * math.log(probability)
        for count, probability in zip(active_counts, probabilities)
    )
    reference_bad = min(
        0.0,
        inner_log + DISTANCE * surprisal - B * log_conditioning,
    )
    value = log_type_count + outer_log + reference_bad
    full_probabilities = [0.0] * len(counts)
    for index, probability in zip(active, probabilities):
        full_probabilities[index] = float(probability)
    return {
        "counts": list(counts),
        "fractions": [count / L for count in counts],
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "reference_type_probabilities": full_probabilities,
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "log2_type_count": log_type_count / math.log(2.0),
        "log2_region_conditioning_probability": log_conditioning / math.log(2.0),
        "optimizer_success": bool(best.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--random-per-concentration", type=int, default=12)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    bands = load_band_witnesses()
    envelope = load_envelope(args.selection, DISTANCE / N)
    compositions = sampled_compositions(
        len(bands), args.random_per_concentration
    )
    rows = []
    for index, counts in enumerate(compositions, start=1):
        row = optimize_composition(envelope, bands, counts)
        rows.append(row)
        if not args.quiet:
            print(
                f"composition={index}/{len(compositions)},margin={row['margin_bits']:.6f},"
                f"counts={counts}",
                flush=True,
            )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-band-compositions-qL-d109-v1",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "sampled_compositions": len(compositions),
        },
        "bands": bands,
        "worst_sampled": worst,
        "all_sampled_below_2^-40": all(
            float(row["margin_bits"]) > 40.0 for row in rows
        ),
        "rows": rows,
        "limitations": [
            "The composition set is sampled and does not cover every integer composition.",
            "Band Bernoulli witnesses are fixed from pure-composition optimization.",
            "Optimization uses binary64; each final transfer evaluation uses 100-digit arithmetic.",
            "The Bernoulli transfer extension above density one half requires a written proof and outward implementation.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst_sampled": worst, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
