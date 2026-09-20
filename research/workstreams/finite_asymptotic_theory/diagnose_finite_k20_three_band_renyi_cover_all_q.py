#!/usr/bin/env python3
"""Adaptive tetrahedral cover of all Q >= 65 three-band compositions."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import Delaunay


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
from diagnose_finite_k20_three_band_renyi_q import optimize_composition  # noqa: E402
from diagnose_finite_k20_three_band_renyi_qL import BANDS  # noqa: E402


MIN_OCCUPATION = 65
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_three_band_renyi_cover_all_q_d109.json"
COMPOSITION_COUNT = math.comb(L + 3, 3) - math.comb(MIN_OCCUPATION + 2, 3)
REQUIRED_POINTWISE_MARGIN = 40.0 + math.log2(COMPOSITION_COUNT)


def log_multinomial_real(counts: np.ndarray) -> float:
    return math.lgamma(L + 1.0) - sum(
        math.lgamma(float(value) + 1.0) for value in counts
    )


def certify_with_witness(envelope, vertices: np.ndarray, witness) -> dict[str, object] | None:
    probabilities = np.asarray(witness["reference_type_probabilities"])
    values = np.asarray(witness["reference_value_probabilities"])
    moments = np.asarray(witness["band_log2_renyi_moments"]) * math.log(2.0)
    surprisal = float(witness["surprisal"])
    order = float(witness["holder_order"])
    conjugate = (order - 1.0) / order
    log_multinomial_coefficient = 1.0 - B * conjugate
    if log_multinomial_coefficient >= 0.0 or np.any(probabilities <= 0.0):
        return None
    bit_probability = float(np.dot(probabilities[1:], values))
    inner = high_precision_transfer_log_moment(envelope, bit_probability, surprisal)
    constant = inner + DISTANCE * surprisal
    rows = []
    for active_counts in vertices:
        occupation_at_vertex = float(np.sum(active_counts))
        counts = np.concatenate(([L - occupation_at_vertex], active_counts))
        log_type_count = log_multinomial_real(counts)
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = constant - B * log_conditioning
        value = (
            log_type_count
            + float(np.dot(active_counts, moments)) / order
            + conjugate * min(0.0, raw)
        )
        rows.append((raw, -value / math.log(2.0)))
    # The raw exponent and the final fixed-witness exponent are convex in
    # the four type counts.  Strict negativity and the margin at every
    # vertex therefore extend throughout the tetrahedron.
    if any(raw >= 0.0 for raw, _ in rows):
        return None
    margin = min(margin for _, margin in rows)
    if margin <= REQUIRED_POINTWISE_MARGIN:
        return None
    return {
        "vertices_active_counts": vertices.tolist(),
        "reference_type_probabilities": probabilities.tolist(),
        "reference_value_probabilities": values.tolist(),
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
        "minimum_vertex_margin_bits": margin,
    }


def certify_tetrahedron(
    envelope, spectrum, vertices: np.ndarray, depth: int
) -> dict[str, object] | None:
    center = np.mean(vertices, axis=0)
    occupation = float(np.sum(center))
    active_counts = tuple(float(value) for value in center)
    quick_witness = optimize_composition(
        envelope,
        spectrum,
        occupation,
        active_counts,
        cover_mode=True,
    )
    receipt = certify_with_witness(envelope, vertices, quick_witness)
    if receipt is not None:
        receipt["optimizer_mode"] = "single-start"
        return receipt
    if depth < 24:
        return None
    robust_witness = optimize_composition(
        envelope,
        spectrum,
        occupation,
        active_counts,
        cover_mode=False,
    )
    receipt = certify_with_witness(envelope, vertices, robust_witness)
    if receipt is not None:
        receipt["optimizer_mode"] = "multistart-fallback"
    return receipt


def subdivide(vertices: np.ndarray) -> list[np.ndarray]:
    edges = [
        (float(np.sum((vertices[left] - vertices[right]) ** 2)), left, right)
        for left in range(4)
        for right in range(left + 1, 4)
    ]
    _, left, right = max(edges)
    midpoint = (vertices[left] + vertices[right]) / 2.0
    first = vertices.copy()
    second = vertices.copy()
    first[left] = midpoint
    second[right] = midpoint
    return [first, second]


def initial_tetrahedra() -> list[np.ndarray]:
    points = np.asarray(
        [
            [float(MIN_OCCUPATION), 0.0, 0.0],
            [0.0, float(MIN_OCCUPATION), 0.0],
            [0.0, 0.0, float(MIN_OCCUPATION)],
            [float(L), 0.0, 0.0],
            [0.0, float(L), 0.0],
            [0.0, 0.0, float(L)],
        ]
    )
    triangulation = Delaunay(points)
    return [points[simplex] for simplex in triangulation.simplices]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--max-depth", type=int, default=50)
    parser.add_argument("--max-tetrahedra", type=int, default=5000)
    parser.add_argument(
        "--resume",
        type=Path,
        help="resume an incomplete v2 receipt containing a pending worklist",
    )
    args = parser.parse_args()
    spectrum = load_spectrum_logs()
    envelope = load_envelope(args.selection, DISTANCE / N)
    if args.resume is None:
        pending = [(tetrahedron, 0) for tetrahedron in initial_tetrahedra()]
        accepted = []
        rejected = []
        processed_before = 0
    else:
        prior = json.loads(args.resume.read_text(encoding="utf-8"))
        if prior.get("schema") != "golay-ba3-rm2sub-finite-b240-three-band-renyi-cover-all-q-d109-v2":
            raise ValueError("resume receipt is not a resumable v2 all-Q cover")
        if prior.get("status") != "INCOMPLETE_DIAGNOSTIC":
            raise ValueError("resume receipt is not incomplete")
        pending = [
            (np.asarray(row["vertices"], dtype=np.float64), int(row["depth"]))
            for row in prior["pending_worklist"]
        ]
        accepted = list(prior["tetrahedra"])
        rejected = list(prior["rejected_tetrahedra"])
        processed_before = int(prior["processed_tetrahedra"])
    processed_this_run = 0
    while pending:
        vertices, depth = pending.pop()
        if processed_this_run >= args.max_tetrahedra:
            pending.append((vertices, depth))
            break
        processed_this_run += 1
        receipt = certify_tetrahedron(envelope, spectrum, vertices, depth)
        if receipt is not None:
            receipt["depth"] = depth
            accepted.append(receipt)
        elif depth < args.max_depth:
            pending.extend((child, depth + 1) for child in subdivide(vertices))
        else:
            rejected.append(
                {"reason": "depth exhausted", "depth": depth, "vertices": vertices.tolist()}
            )
        if processed_this_run % 25 == 0:
            print(
                f"processed_this_run={processed_this_run},accepted={len(accepted)},"
                f"pending={len(pending)},"
                f"rejected={len(rejected)}",
                flush=True,
            )
    processed = processed_before + processed_this_run
    minimum = min(
        float(row["minimum_vertex_margin_bits"]) for row in accepted
    ) if accepted else -math.inf
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-three-band-renyi-cover-all-q-d109-v2",
        "status": (
            "BINARY64_DIAGNOSTIC_CONVEX_TETRAHEDRAL_COVER"
            if not rejected and not pending
            else "INCOMPLETE_DIAGNOSTIC"
        ),
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "minimum_occupation": MIN_OCCUPATION,
            "maximum_occupation": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "composition_count": COMPOSITION_COUNT,
            "composition_count_log2": math.log2(COMPOSITION_COUNT),
            "required_pointwise_margin_bits": REQUIRED_POINTWISE_MARGIN,
        },
        "processed_tetrahedra": processed,
        "accepted_tetrahedra": len(accepted),
        "pending_tetrahedra": len(pending),
        "pending_worklist": [
            {"vertices": vertices.tolist(), "depth": depth}
            for vertices, depth in pending
        ],
        "rejected_tetrahedra": rejected,
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": minimum - math.log2(COMPOSITION_COUNT),
        "tetrahedra": accepted,
        "limitations": [
            "Optimization and witnesses use nearest binary64 and are not outward rounded.",
            "Final transfer evaluations use 100-digit arithmetic but not intervals.",
            "The inactive-category Renyi reduction, convexity argument, and density-above-half extension require written audits.",
            "Occupations Q=1,...,64 are covered by separate outward receipts, not this cover.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "processed_tetrahedra": processed,
                "accepted_tetrahedra": len(accepted),
                "pending_tetrahedra": len(pending),
                "rejected_tetrahedra": len(rejected),
                "minimum_accepted_vertex_margin_bits": minimum,
                "aggregate_margin_from_max_term_bits": payload[
                    "aggregate_margin_from_max_term_bits"
                ],
            },
            indent=2,
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
