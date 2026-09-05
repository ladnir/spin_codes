#!/usr/bin/env python3
"""Barycentric simplex cover of all five-band compositions at Q=L, d=10.9%."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import Delaunay


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import (  # noqa: E402
    DISTANCE,
    PURE_RECEIPT,
    load_band_witnesses,
)
from diagnose_finite_k20_band_cover_qL import (  # noqa: E402
    COMPOSITION_COUNT,
    REQUIRED_POINTWISE_MARGIN,
    log_multinomial_real,
    optimize_witness,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    L,
    N,
    high_precision_transfer_log_moment,
)


OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_simplex_cover_qL_d109.json"


def initial_simplices(grid: int) -> list[np.ndarray]:
    points = []
    for first in range(grid + 1):
        for second in range(grid + 1 - first):
            for third in range(grid + 1 - first - second):
                for fourth in range(grid + 1 - first - second - third):
                    points.append((first, second, third, fourth))
    coordinates = np.asarray(points, dtype=np.float64)
    triangulation = Delaunay(coordinates, qhull_options="QJ")
    result = []
    scale = L / grid
    for indices in triangulation.simplices:
        first_four = coordinates[indices] * scale
        fifth = L - np.sum(first_four, axis=1)
        vertices = np.column_stack((first_four, fifth))
        if np.min(vertices) >= -1e-7:
            result.append(np.maximum(vertices, 0.0))
    return result


def certify_simplex(envelope, bands, vertices: np.ndarray) -> dict[str, object] | None:
    center = np.mean(vertices, axis=0)
    probabilities, surprisal = optimize_witness(envelope, bands, center)
    values = np.asarray([float(band["value_probability"]) for band in bands])
    outer = np.asarray([float(band["log_majorant"]) for band in bands])
    bit_probability = float(np.dot(probabilities, values))
    inner = high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
    constant = inner + DISTANCE * surprisal
    rows = []
    for counts in vertices:
        log_type_count = log_multinomial_real(counts)
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = constant - B * log_conditioning
        value = log_type_count + float(np.dot(counts, outer)) + min(0.0, raw)
        rows.append((raw, -value / math.log(2.0)))
    if any(raw >= 0.0 for raw, _ in rows):
        return None
    margin = min(margin for _, margin in rows)
    if margin <= REQUIRED_POINTWISE_MARGIN:
        return None
    return {
        "vertices": [[float(value) for value in row] for row in vertices],
        "reference_type_probabilities": [float(value) for value in probabilities],
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "minimum_vertex_margin_bits": margin,
    }


def subdivide(vertices: np.ndarray) -> list[np.ndarray]:
    longest = None
    for left in range(len(vertices)):
        for right in range(left + 1, len(vertices)):
            length = float(np.sum((vertices[left] - vertices[right]) ** 2))
            if longest is None or length > longest[0]:
                longest = (length, left, right)
    assert longest is not None
    _, left, right = longest
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
    parser.add_argument("--grid", type=int, default=4)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--max-simplices", type=int, default=20_000)
    args = parser.parse_args()
    bands = load_band_witnesses()
    envelope = load_envelope(args.selection, DISTANCE / N)
    initial = initial_simplices(args.grid)
    pending = [(vertices, 0) for vertices in initial]
    accepted = []
    rejected = []
    processed = 0
    while pending:
        vertices, depth = pending.pop()
        processed += 1
        if processed > args.max_simplices:
            raise RuntimeError("simplex budget exhausted")
        receipt = certify_simplex(envelope, bands, vertices)
        if receipt is not None:
            receipt["depth"] = depth
            accepted.append(receipt)
        elif depth < args.max_depth:
            pending.extend((child, depth + 1) for child in subdivide(vertices))
        else:
            rejected.append(
                {
                    "vertices": vertices.tolist(),
                    "depth": depth,
                }
            )
        if processed % 100 == 0:
            print(
                f"processed={processed},accepted={len(accepted)},pending={len(pending)},"
                f"rejected={len(rejected)}",
                flush=True,
            )
    minimum = min(
        float(row["minimum_vertex_margin_bits"]) for row in accepted
    ) if accepted else -math.inf
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-band-simplex-cover-qL-d109-v1",
        "status": "BINARY64_DIAGNOSTIC_CONVEX_SIMPLEX_COVER" if not rejected else "INCOMPLETE_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "initial_grid": args.grid,
            "initial_simplices": len(initial),
            "composition_count": COMPOSITION_COUNT,
            "composition_count_log2": math.log2(COMPOSITION_COUNT),
            "required_pointwise_margin_bits": REQUIRED_POINTWISE_MARGIN,
        },
        "source_pure_receipt": str(PURE_RECEIPT),
        "processed_simplices": processed,
        "accepted_simplices": len(accepted),
        "rejected_simplices": rejected,
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": minimum - math.log2(COMPOSITION_COUNT),
        "simplices": accepted,
        "limitations": [
            "The Delaunay triangulation topology is generated in binary64 and requires exact rational reconstruction.",
            "Optimization and accepted witnesses are not outward rounded.",
            "Final inner evaluations use 100-digit arithmetic but not intervals.",
            "The convexity and Bernoulli-density extension require a written proof audit.",
            "This cover applies only to Q=L.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "processed_simplices": processed,
        "accepted_simplices": len(accepted),
        "rejected_simplices": len(rejected),
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": payload["aggregate_margin_from_max_term_bits"],
    }, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
