#!/usr/bin/env python3
"""Adaptive triangle cover for the three-band Renyi Q=L bound at 10.9%."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import DISTANCE  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    L,
    N,
    high_precision_transfer_log_moment,
    load_spectrum_logs,
)
from diagnose_finite_k20_three_band_renyi_qL import (  # noqa: E402
    BANDS,
    optimize_composition,
)


OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_three_band_renyi_cover_qL_d109.json"
COMPOSITION_COUNT = math.comb(L + 2, 2)
REQUIRED_POINTWISE_MARGIN = 40.0 + math.log2(COMPOSITION_COUNT)


def log_multinomial_real(counts: np.ndarray) -> float:
    return math.lgamma(L + 1.0) - sum(
        math.lgamma(float(value) + 1.0) for value in counts
    )


def certify_triangle(envelope, spectrum, vertices: np.ndarray) -> dict[str, object] | None:
    center = np.mean(vertices, axis=0)
    witness = optimize_composition(envelope, spectrum, tuple(center))
    probabilities = np.asarray(witness["reference_type_probabilities"])
    values = np.asarray(witness["reference_value_probabilities"])
    moments = np.asarray(witness["band_log2_renyi_moments"]) * math.log(2.0)
    surprisal = float(witness["surprisal"])
    order = float(witness["holder_order"])
    conjugate = (order - 1.0) / order
    log_multinomial_coefficient = 1.0 - B * conjugate
    if log_multinomial_coefficient >= 0.0 or np.any(probabilities <= 0.0):
        return None
    bit_probability = float(np.dot(probabilities, values))
    inner = high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
    constant = inner + DISTANCE * surprisal
    rows = []
    for counts in vertices:
        log_type_count = log_multinomial_real(counts)
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = constant - B * log_conditioning
        value = (
            log_type_count
            + float(np.dot(counts, moments)) / order
            + conjugate * min(0.0, raw)
        )
        rows.append((raw, -value / math.log(2.0)))
    if any(raw >= 0.0 for raw, _ in rows):
        return None
    margin = min(margin for _, margin in rows)
    if margin <= REQUIRED_POINTWISE_MARGIN:
        return None
    return {
        "vertices": vertices.tolist(),
        "reference_type_probabilities": probabilities.tolist(),
        "reference_value_probabilities": values.tolist(),
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
        "minimum_vertex_margin_bits": margin,
    }


def subdivide(vertices: np.ndarray) -> list[np.ndarray]:
    edges = [
        (float(np.sum((vertices[left] - vertices[right]) ** 2)), left, right)
        for left in range(3)
        for right in range(left + 1, 3)
    ]
    _, left, right = max(edges)
    midpoint = (vertices[left] + vertices[right]) / 2.0
    first = vertices.copy()
    second = vertices.copy()
    first[left] = midpoint
    second[right] = midpoint
    return [first, second]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--max-triangles", type=int, default=5000)
    args = parser.parse_args()
    spectrum = load_spectrum_logs()
    envelope = load_envelope(args.selection, DISTANCE / N)
    initial = np.asarray(
        [[float(L), 0.0, 0.0], [0.0, float(L), 0.0], [0.0, 0.0, float(L)]]
    )
    pending = [(initial, 0)]
    accepted = []
    rejected = []
    processed = 0
    while pending:
        vertices, depth = pending.pop()
        processed += 1
        if processed > args.max_triangles:
            rejected.append({"reason": "triangle budget exhausted", "pending": len(pending) + 1})
            break
        receipt = certify_triangle(envelope, spectrum, vertices)
        if receipt is not None:
            receipt["depth"] = depth
            accepted.append(receipt)
        elif depth < args.max_depth:
            pending.extend((child, depth + 1) for child in subdivide(vertices))
        else:
            rejected.append({"reason": "depth exhausted", "depth": depth, "vertices": vertices.tolist()})
        if processed % 25 == 0:
            print(
                f"processed={processed},accepted={len(accepted)},pending={len(pending)},"
                f"rejected={len(rejected)}",
                flush=True,
            )
    minimum = min(
        float(row["minimum_vertex_margin_bits"]) for row in accepted
    ) if accepted else -math.inf
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-three-band-renyi-cover-qL-d109-v1",
        "status": "BINARY64_DIAGNOSTIC_CONVEX_TRIANGLE_COVER" if not rejected and not pending else "INCOMPLETE_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "composition_count": COMPOSITION_COUNT,
            "composition_count_log2": math.log2(COMPOSITION_COUNT),
            "required_pointwise_margin_bits": REQUIRED_POINTWISE_MARGIN,
        },
        "processed_triangles": processed,
        "accepted_triangles": len(accepted),
        "pending_triangles": len(pending),
        "rejected_triangles": rejected,
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": minimum - math.log2(COMPOSITION_COUNT),
        "triangles": accepted,
        "limitations": [
            "Optimization and witnesses use nearest binary64 and are not outward rounded.",
            "Final transfer evaluations use 100-digit arithmetic but not intervals.",
            "The multiband Renyi reduction, convexity argument, and density-above-half extension require written audits.",
            "This cover applies only to Q=L.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "processed_triangles": processed,
        "accepted_triangles": len(accepted),
        "pending_triangles": len(pending),
        "rejected_triangles": len(rejected),
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": payload["aggregate_margin_from_max_term_bits"],
    }, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
