#!/usr/bin/env python3
"""Adaptive gap-free diagnostic cover of q=65..8192 and all band mixtures."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_renyi_dense import high_precision_transfer_log_moment  # noqa: E402
from diagnose_ebch32_ba_k20_three_band_qL import (  # noqa: E402
    B,
    BANDS,
    EPOCHS,
    L,
    N,
    conditioned_spectrum,
)
from diagnose_ebch32_parityfanout_ba_three_band_q import (  # noqa: E402
    DISTANCE,
    optimize_composition,
)


if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") != "1":
    raise RuntimeError("set SPIN_EBCH_PARITY_FANOUT=1")
if os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22") != "24":
    raise RuntimeError("set SPIN_EBCH_LOWER_WEIGHT=24")

MIN_OCCUPATION = 65
OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_DENSE_COVER_OUTPUT",
    "ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json",
)
COMPOSITION_COUNT = math.comb(L + 3, 3) - math.comb(MIN_OCCUPATION + 2, 3)
REQUIRED_POINTWISE_MARGIN = 40.0 + math.log2(COMPOSITION_COUNT)
SCHEMA = (
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v11"
    if DISTANCE == 230_686
    else f"ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d{DISTANCE}-v1"
)
LEGACY_SCHEMAS = {
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v1",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v2",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v3",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v4",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v5",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v6",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v7",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v8",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v9",
    "ebch32-parityfanout31x33-ba3-b256-three-band-cover-all-q-d11-v10",
}
CELL_CONTRIBUTION_MARGIN = 44.0
MAX_ACCEPTED_CELLS = 8192
LATTICE_ENUMERATION_LIMIT = 4
DELEGATED_EXTERNAL_TYPE = (3344, 1, 21)
DELEGATED_LOG2_CONTRIBUTION_UPPER = -198.0


def log_multinomial_real(counts: np.ndarray) -> float:
    return math.lgamma(L + 1.0) - sum(
        math.lgamma(float(value) + 1.0) for value in counts
    )


def bounding_box_lattice_count(vertices: np.ndarray) -> int:
    """Upper-bound integer active-count triples in a tetrahedron."""
    result = 1
    for coordinate in range(3):
        lower = math.ceil(float(np.min(vertices[:, coordinate])))
        upper = math.floor(float(np.max(vertices[:, coordinate])))
        result *= max(0, upper - lower + 1)
    return result


def bounding_box_points(vertices: np.ndarray):
    ranges = []
    for coordinate in range(3):
        lower = math.ceil(float(np.min(vertices[:, coordinate])))
        upper = math.floor(float(np.max(vertices[:, coordinate])))
        ranges.append(range(lower, upper + 1))
    for first in ranges[0]:
        for second in ranges[1]:
            for third in ranges[2]:
                occupation = first + second + third
                if MIN_OCCUPATION <= occupation <= L:
                    yield (first, second, third)


def delegated_point_receipt(vertices: np.ndarray) -> dict[str, object] | None:
    lower = tuple(
        math.ceil(float(np.min(vertices[:, coordinate])))
        for coordinate in range(3)
    )
    upper = tuple(
        math.floor(float(np.max(vertices[:, coordinate])))
        for coordinate in range(3)
    )
    if lower != DELEGATED_EXTERNAL_TYPE or upper != DELEGATED_EXTERNAL_TYPE:
        return None
    return {
        "vertices_active_counts": vertices.tolist(),
        "optimizer_mode": "delegated-split-low-point",
        "delegated_external_type": list(DELEGATED_EXTERNAL_TYPE),
        "minimum_vertex_margin_bits": -DELEGATED_LOG2_CONTRIBUTION_UPPER,
        "bounding_box_lattice_count": 1,
        "bounding_box_log2_contribution_upper": DELEGATED_LOG2_CONTRIBUTION_UPPER,
        "delegated_outward_receipt": (
            "ebch32_parityfanout31x33_ba3_B256_"
            "split_low_rejected_point_cover_outward.json"
        ),
    }


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
    inner = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal, epochs=EPOCHS
    )
    constant = inner + DISTANCE * surprisal
    rows = []
    for active_counts in vertices:
        occupation = float(np.sum(active_counts))
        counts = np.concatenate(([L - occupation], active_counts))
        log_type_count = log_multinomial_real(counts)
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = constant - B * log_conditioning
        value = (
            log_type_count
            + float(np.dot(active_counts, moments)) / order
            + conjugate * min(0.0, raw)
        )
        rows.append((raw, -value / math.log(2.0)))
    if any(raw >= 0.0 for raw, _ in rows):
        return None
    margin = min(margin for _, margin in rows)
    lattice_count = bounding_box_lattice_count(vertices)
    log2_lattice_count = math.log2(lattice_count) if lattice_count else -math.inf
    # The box can include triples outside the tetrahedron and triples whose
    # occupation lies outside [65,L].  This only enlarges the union bound.
    # The per-cell 2^-44 threshold controls refinement.  The authoritative
    # 40-bit test is the explicit sum of all accepted box contributions.
    required_margin = CELL_CONTRIBUTION_MARGIN + log2_lattice_count
    if lattice_count and margin <= required_margin:
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
        "bounding_box_lattice_count": lattice_count,
        "bounding_box_log2_contribution_upper": (
            log2_lattice_count - margin if lattice_count else None
        ),
        "required_vertex_margin_bits": required_margin,
    }


def certify_tetrahedron(
    envelope,
    spectrum,
    vertices: np.ndarray,
    depth: int,
    seed_witnesses: tuple[dict[str, object], ...],
):
    lattice_count = bounding_box_lattice_count(vertices)
    if lattice_count == 0:
        return {
            "vertices_active_counts": vertices.tolist(),
            "optimizer_mode": "empty-lattice-box",
            "minimum_vertex_margin_bits": math.inf,
            "bounding_box_lattice_count": 0,
            "bounding_box_log2_contribution_upper": None,
            "point_witnesses": [],
        }, ()
    delegated = delegated_point_receipt(vertices)
    if delegated is not None:
        return delegated, ()
    for inherited in seed_witnesses:
        receipt = certify_with_witness(envelope, vertices, inherited)
        if receipt is not None:
            receipt["optimizer_mode"] = "inherited-witness"
            return receipt, ()
    if lattice_count <= LATTICE_ENUMERATION_LIMIT:
        point_receipts = []
        point_logs = []
        for point in bounding_box_points(vertices):
            point_vertices = np.repeat(
                np.asarray(point, dtype=np.float64)[None, :], 4, axis=0
            )
            occupation = float(sum(point))
            witness = optimize_composition(
                envelope,
                spectrum,
                occupation,
                point,
                cover_mode=True,
                seed_witnesses=seed_witnesses,
            )
            point_receipt = certify_with_witness(
                envelope, point_vertices, witness
            )
            if point_receipt is None:
                witness = optimize_composition(
                    envelope,
                    spectrum,
                    occupation,
                    point,
                    cover_mode=False,
                    seed_witnesses=(witness, *seed_witnesses),
                )
                point_receipt = certify_with_witness(
                    envelope, point_vertices, witness
                )
            if point_receipt is None:
                return None, (witness,)
            point_receipts.append(point_receipt)
            point_logs.append(
                float(point_receipt["bounding_box_log2_contribution_upper"])
            )
        if not point_logs:
            return {
                "vertices_active_counts": vertices.tolist(),
                "optimizer_mode": "empty-admissible-lattice-box",
                "minimum_vertex_margin_bits": math.inf,
                "bounding_box_lattice_count": 0,
                "bounding_box_log2_contribution_upper": None,
                "point_witnesses": [],
            }, ()
        largest = max(point_logs)
        aggregate_log2 = largest + math.log2(
            sum(2.0 ** (value - largest) for value in point_logs)
        )
        if aggregate_log2 > -CELL_CONTRIBUTION_MARGIN:
            return None, tuple(
                row for row in point_receipts[:2]
            )
        return {
            "vertices_active_counts": vertices.tolist(),
            "optimizer_mode": "enumerated-lattice-box",
            "minimum_vertex_margin_bits": min(
                float(row["minimum_vertex_margin_bits"])
                for row in point_receipts
            ),
            "bounding_box_lattice_count": len(point_receipts),
            "bounding_box_log2_contribution_upper": aggregate_log2,
            "point_witnesses": point_receipts,
        }, ()
    center = np.mean(vertices, axis=0)
    occupation = float(np.sum(center))
    active_counts = tuple(float(value) for value in center)
    quick = optimize_composition(
        envelope,
        spectrum,
        occupation,
        active_counts,
        cover_mode=True,
        seed_witnesses=seed_witnesses,
    )
    receipt = certify_with_witness(envelope, vertices, quick)
    if receipt is not None:
        receipt["optimizer_mode"] = "single-start"
        return receipt, ()
    return None, (quick,)


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
    # Explicit pulling triangulation of
    # {x >= 0 : MIN_OCCUPATION <= sum(x) <= L}.  Keeping the roots explicit
    # removes a numerical Delaunay implementation from the certificate base.
    low_x = [float(MIN_OCCUPATION), 0.0, 0.0]
    low_y = [0.0, float(MIN_OCCUPATION), 0.0]
    low_z = [0.0, 0.0, float(MIN_OCCUPATION)]
    high_x = [float(L), 0.0, 0.0]
    high_y = [0.0, float(L), 0.0]
    high_z = [0.0, 0.0, float(L)]
    return [
        np.asarray([low_z, high_z, high_y, high_x]),
        np.asarray([low_z, high_y, high_x, low_y]),
        np.asarray([low_z, low_x, high_x, low_y]),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--max-depth", type=int, default=50)
    parser.add_argument("--max-tetrahedra", type=int, default=5000)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    spectrum, good = conditioned_spectrum()
    envelope = load_envelope(args.selection, DISTANCE / N)
    if args.resume is None:
        pending = [(tetrahedron, 0, ()) for tetrahedron in initial_tetrahedra()]
        accepted = []
        rejected = []
        processed_before = 0
    else:
        prior = json.loads(args.resume.read_text(encoding="utf-8"))
        if prior.get("schema") not in ({SCHEMA} | LEGACY_SCHEMAS) or prior.get("status") != "INCOMPLETE_DIAGNOSTIC":
            raise ValueError("resume receipt is not a compatible incomplete cover")
        pending = [
            (
                np.asarray(row["vertices"], dtype=np.float64),
                int(row["depth"]),
                tuple(row.get("seed_witnesses", ())),
            )
            for row in prior["pending_worklist"]
        ]
        # Resume the shallowest saved branch first.  Children created during
        # this run remain depth-first, so one broad branch is completed before
        # the next saved branch begins.
        pending.sort(key=lambda item: item[1], reverse=True)
        accepted = list(prior["tetrahedra"])
        rejected = []
        for row in prior["rejected_tetrahedra"]:
            vertices = np.asarray(row["vertices"], dtype=np.float64)
            delegated = delegated_point_receipt(vertices)
            if delegated is None:
                rejected.append(row)
                continue
            delegated["depth"] = int(row["depth"])
            accepted.append(delegated)
        processed_before = int(prior["processed_tetrahedra"])
    processed_this_run = 0
    while pending:
        vertices, depth, seed_witnesses = pending.pop()
        if processed_this_run >= args.max_tetrahedra:
            pending.append((vertices, depth, seed_witnesses))
            break
        processed_this_run += 1
        receipt, child_seeds = certify_tetrahedron(
            envelope, spectrum, vertices, depth, seed_witnesses
        )
        if receipt is not None:
            receipt["depth"] = depth
            accepted.append(receipt)
        elif depth < args.max_depth:
            pending.extend(
                (child, depth + 1, child_seeds)
                for child in subdivide(vertices)
            )
        else:
            rejected.append(
                {"reason": "depth exhausted", "depth": depth, "vertices": vertices.tolist()}
            )
        if processed_this_run % 25 == 0:
            print(
                f"processed_this_run={processed_this_run},accepted={len(accepted)},"
                f"pending={len(pending)},rejected={len(rejected)}",
                flush=True,
            )
    processed = processed_before + processed_this_run
    minimum = min(
        float(row["minimum_vertex_margin_bits"]) for row in accepted
    ) if accepted else -math.inf
    cell_log2_contributions = []
    for row in accepted:
        vertices = np.asarray(row["vertices_active_counts"], dtype=np.float64)
        lattice_count = bounding_box_lattice_count(vertices)
        if row.get("optimizer_mode") in (
            "enumerated-lattice-box",
            "empty-lattice-box",
            "empty-admissible-lattice-box",
        ):
            lattice_count = int(row["bounding_box_lattice_count"])
        else:
            row["bounding_box_lattice_count"] = lattice_count
        if lattice_count:
            log2_contribution = (
                float(row["bounding_box_log2_contribution_upper"])
                if row.get("optimizer_mode") == "enumerated-lattice-box"
                else math.log2(lattice_count)
                - float(row["minimum_vertex_margin_bits"])
            )
            row["bounding_box_log2_contribution_upper"] = log2_contribution
            cell_log2_contributions.append(log2_contribution)
        else:
            row["bounding_box_log2_contribution_upper"] = None
    if cell_log2_contributions:
        largest_cell = max(cell_log2_contributions)
        aggregate_log2 = largest_cell + math.log2(
            sum(2.0 ** (value - largest_cell) for value in cell_log2_contributions)
        )
    else:
        aggregate_log2 = -math.inf
    nonempty_minimum = min(
        float(row["minimum_vertex_margin_bits"])
        for row in accepted
        if int(row["bounding_box_lattice_count"]) > 0
    ) if cell_log2_contributions else math.inf
    diagnostic_success = (
        not rejected
        and not pending
        and len(accepted) <= MAX_ACCEPTED_CELLS
        and aggregate_log2 <= -40.0
    )
    payload = {
        "schema": SCHEMA,
        "status": (
            "BINARY64_DIAGNOSTIC_CONVEX_TETRAHEDRAL_COVER"
            if diagnostic_success
            else "COMPLETE_DIAGNOSTIC_BOUND_FAILED"
            if not rejected and not pending
            else "INCOMPLETE_DIAGNOSTIC"
        ),
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "conditioning_good_probability_diagnostic": good,
            "bands": [list(band) for band in BANDS],
            "composition_count": COMPOSITION_COUNT,
            "composition_count_log2": math.log2(COMPOSITION_COUNT),
            "required_pointwise_margin_bits": REQUIRED_POINTWISE_MARGIN,
            "cell_contribution_margin_bits": CELL_CONTRIBUTION_MARGIN,
            "maximum_accepted_cells_resource_guard": MAX_ACCEPTED_CELLS,
            "lattice_enumeration_box_limit": LATTICE_ENUMERATION_LIMIT,
            "delegated_external_type": list(DELEGATED_EXTERNAL_TYPE),
            "delegated_log2_contribution_upper": DELEGATED_LOG2_CONTRIBUTION_UPPER,
            "cell_accounting": (
                "Each cell is charged the number of integer triples in its "
                "axis-aligned coordinate bounding box. Boxes may overlap and "
                "may include inadmissible triples; both effects overcount."
            ),
        },
        "processed_tetrahedra": processed,
        "accepted_tetrahedra": len(accepted),
        "pending_tetrahedra": len(pending),
        "pending_worklist": [
            {
                "vertices": vertices.tolist(),
                "depth": depth,
                "seed_witnesses": list(seed_witnesses),
            }
            for vertices, depth, seed_witnesses in pending
        ],
        "rejected_tetrahedra": rejected,
        "minimum_accepted_vertex_margin_bits": minimum,
        "minimum_nonempty_box_vertex_margin_bits": nonempty_minimum,
        "aggregate_margin_from_max_term_bits": minimum - math.log2(COMPOSITION_COUNT),
        "aggregate_log2_upper_from_cell_boxes": aggregate_log2,
        "aggregate_margin_from_cell_boxes_bits": -aggregate_log2,
        "tetrahedra": accepted,
        "limitations": [
            "Optimization and witnesses use nearest binary64 and are not outward rounded.",
            "Final transfer evaluations use 100-digit arithmetic but not intervals.",
            "A completed diagnostic cover still requires an outward verifier.",
            "Cell-box counts are exact integers, but the displayed aggregate uses nearest binary64.",
            "The inactive-category Renyi reduction and convexity argument require a written audit.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "processed_tetrahedra": processed,
        "accepted_tetrahedra": len(accepted),
        "pending_tetrahedra": len(pending),
        "rejected_tetrahedra": len(rejected),
        "minimum_accepted_vertex_margin_bits": minimum,
        "aggregate_margin_from_max_term_bits": payload["aggregate_margin_from_max_term_bits"],
        "aggregate_margin_from_cell_boxes_bits": payload["aggregate_margin_from_cell_boxes_bits"],
    }, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
