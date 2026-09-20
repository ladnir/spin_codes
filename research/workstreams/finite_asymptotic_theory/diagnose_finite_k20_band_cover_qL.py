#!/usr/bin/env python3
"""Adaptive convex-box cover of all five-band compositions at Q=L, d=10.9%."""

from __future__ import annotations

import argparse
import itertools
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
    PURE_RECEIPT,
    fast_inner_log,
    load_band_witnesses,
    softmax,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    L,
    N,
    high_precision_transfer_log_moment,
)


OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_cover_qL_d109.json"
COMPOSITION_COUNT = math.comb(L + 4, 4)
REQUIRED_POINTWISE_MARGIN = 40.0 + math.log2(COMPOSITION_COUNT)


def log_multinomial_real(counts: np.ndarray) -> float:
    return math.lgamma(L + 1.0) - sum(math.lgamma(float(value) + 1.0) for value in counts)


def polytope_vertices(lower: np.ndarray, upper: np.ndarray) -> list[np.ndarray]:
    vertices: list[np.ndarray] = []
    for selectors in itertools.product((0, 1), repeat=4):
        point = np.asarray(
            [upper[i] if selectors[i] else lower[i] for i in range(4)],
            dtype=np.float64,
        )
        if float(np.sum(point)) <= L + 1e-8:
            vertices.append(point)
    for free in range(4):
        fixed = [index for index in range(4) if index != free]
        for selectors in itertools.product((0, 1), repeat=3):
            point = np.zeros(4, dtype=np.float64)
            for index, selector in zip(fixed, selectors):
                point[index] = upper[index] if selector else lower[index]
            point[free] = L - float(np.sum(point))
            if lower[free] - 1e-8 <= point[free] <= upper[free] + 1e-8:
                point[free] = min(upper[free], max(lower[free], point[free]))
                vertices.append(point)
    unique: dict[tuple[float, ...], np.ndarray] = {}
    for point in vertices:
        key = tuple(round(float(value), 9) for value in point)
        unique[key] = point
    return list(unique.values())


def full_counts(point: np.ndarray) -> np.ndarray:
    final = L - float(np.sum(point))
    if final < -1e-7:
        raise ArithmeticError("point lies outside the composition simplex")
    return np.concatenate((point, [max(0.0, final)]))


def optimize_witness(envelope, bands: list[dict[str, object]], counts: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray([float(band["value_probability"]) for band in bands])
    outer = np.asarray([float(band["log_majorant"]) for band in bands])
    log_type_count = log_multinomial_real(counts)

    def objective(point: np.ndarray) -> float:
        probabilities = softmax(point[:-1])
        if np.any(probabilities <= 0.0):
            return math.inf
        surprisal = math.exp(min(4.0, max(-8.0, float(point[-1]))))
        bit_probability = float(np.dot(probabilities, values))
        inner = fast_inner_log(envelope, bit_probability, surprisal)
        if not math.isfinite(inner):
            return math.inf
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = inner + DISTANCE * surprisal - B * log_conditioning
        return log_type_count + float(np.dot(counts, outer)) + min(0.0, raw)

    frequencies = (counts + 0.5) / (L + 2.5)
    logits = np.log(frequencies)
    starts = [
        np.concatenate((logits, [math.log(surprisal)]))
        for surprisal in (0.7, 1.3, 2.0, 2.7)
    ]
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"maxiter": 600, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    return softmax(best.x[:-1]), math.exp(float(best.x[-1]))


def certify_box_diagnostic(
    envelope,
    bands: list[dict[str, object]],
    lower: np.ndarray,
    upper: np.ndarray,
) -> dict[str, object] | None:
    vertices4 = polytope_vertices(lower, upper)
    if not vertices4:
        return None
    vertices = [full_counts(point) for point in vertices4]
    center = np.mean(np.asarray(vertices), axis=0)
    probabilities, surprisal = optimize_witness(envelope, bands, center)
    values = np.asarray([float(band["value_probability"]) for band in bands])
    outer = np.asarray([float(band["log_majorant"]) for band in bands])
    bit_probability = float(np.dot(probabilities, values))
    inner = high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
    constant = inner + DISTANCE * surprisal
    vertex_rows = []
    for counts in vertices:
        log_type_count = log_multinomial_real(counts)
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = constant - B * log_conditioning
        value = log_type_count + float(np.dot(counts, outer)) + min(0.0, raw)
        vertex_rows.append(
            {
                "counts": [float(item) for item in counts],
                "raw_reference_log_upper": raw,
                "margin_bits": -value / math.log(2.0),
            }
        )
    if any(float(row["raw_reference_log_upper"]) >= 0.0 for row in vertex_rows):
        return None
    minimum_margin = min(float(row["margin_bits"]) for row in vertex_rows)
    if minimum_margin <= REQUIRED_POINTWISE_MARGIN:
        return None
    return {
        "lower_first_four_counts": [float(value) for value in lower],
        "upper_first_four_counts": [float(value) for value in upper],
        "reference_type_probabilities": [float(value) for value in probabilities],
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "vertex_count": len(vertex_rows),
        "minimum_vertex_margin_bits": minimum_margin,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--max-boxes", type=int, default=20_000)
    parser.add_argument("--max-depth", type=int, default=18)
    args = parser.parse_args()
    bands = load_band_witnesses()
    envelope = load_envelope(args.selection, DISTANCE / N)
    pending = [(np.zeros(4), np.full(4, float(L)), 0)]
    accepted = []
    rejected = []
    processed = 0
    while pending:
        lower, upper, depth = pending.pop()
        if not polytope_vertices(lower, upper):
            continue
        processed += 1
        if processed > args.max_boxes:
            raise RuntimeError("box budget exhausted")
        receipt = certify_box_diagnostic(envelope, bands, lower, upper)
        if receipt is not None:
            receipt["depth"] = depth
            accepted.append(receipt)
            continue
        widths = upper - lower
        split_dimension = int(np.argmax(widths))
        if depth >= args.max_depth or widths[split_dimension] < 1.0:
            rejected.append(
                {
                    "lower": lower.tolist(),
                    "upper": upper.tolist(),
                    "depth": depth,
                }
            )
            continue
        middle = (lower[split_dimension] + upper[split_dimension]) / 2.0
        left_upper = upper.copy()
        left_upper[split_dimension] = middle
        right_lower = lower.copy()
        right_lower[split_dimension] = middle
        pending.append((right_lower, upper.copy(), depth + 1))
        pending.append((lower.copy(), left_upper, depth + 1))
        if processed % 10 == 0:
            print(
                f"processed={processed},accepted={len(accepted)},pending={len(pending)},"
                f"rejected={len(rejected)}",
                flush=True,
            )
    minimum = min(
        float(row["minimum_vertex_margin_bits"]) for row in accepted
    ) if accepted else -math.inf
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-band-cover-qL-d109-v1",
        "status": "BINARY64_DIAGNOSTIC_CONVEX_BOX_COVER" if not rejected else "INCOMPLETE_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "composition_count": COMPOSITION_COUNT,
            "composition_count_log2": math.log2(COMPOSITION_COUNT),
            "required_pointwise_margin_bits": REQUIRED_POINTWISE_MARGIN,
        },
        "source_pure_receipt": str(PURE_RECEIPT),
        "processed_boxes": processed,
        "accepted_boxes": len(accepted),
        "rejected_boxes": rejected,
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": minimum - math.log2(COMPOSITION_COUNT),
        "boxes": accepted,
        "limitations": [
            "Optimization uses nearest binary64 and accepted witnesses are not outward rounded.",
            "Each final inner transfer uses 100-digit arithmetic but not interval arithmetic.",
            "The convexity and Bernoulli-density extension require a written proof audit.",
            "This cover applies only to Q=L.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: payload[key] for key in (
        "status", "processed_boxes", "accepted_boxes", "minimum_accepted_vertex_margin_bits", "aggregate_margin_from_max_term_bits"
    )}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
