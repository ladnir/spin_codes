#!/usr/bin/env python3
"""Build and exactly certify a regular triangulation of g=4 profile anchors.

Qhull sees normalized binary64 lifted points and proposes lower facets.  It is
not trusted: each proposed facet is reconstructed over ``Fraction``, must be a
strict supporting lower hyperplane with no extra coplanar anchor, and the
projected cells must pass the independent exact facet/boundary/volume audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.spatial import ConvexHull, QhullError

import certify_packet_group_g4_anchor_mesh as verifier


SCHEMA = "packet-group-g4-exact-regular-root-mesh-v1"
PERTURB_PRIME = 1_000_003
PERTURB_SCALE = 1 << 20


_WORKER_ANCHORS: tuple[tuple[Fraction, ...], ...] | None = None
_WORKER_HEIGHTS: tuple[Fraction, ...] | None = None


def lifting_height(profile: tuple[Fraction, ...], rank: int) -> Fraction:
    if not 0 <= rank < PERTURB_PRIME - 1:
        raise ValueError("regular-mesh rank exceeds deterministic perturbation range")
    squared_norm = sum((profile[index] / verifier.M) ** 2 for index in range(1, 5))
    # Since gcd(5, PERTURB_PRIME-1)=1, x -> x^5 permutes nonzero residues.
    # This gives distinct nonlinear deterministic perturbations below 2^-20.
    # Exact strict-facet and topology audits decide whether the finite
    # perturbation was sufficient; Qhull never gets that authority.
    residue = pow(rank + 1, 5, PERTURB_PRIME)
    return squared_norm + Fraction(residue, PERTURB_PRIME * PERTURB_SCALE)


def solve_linear(
    matrix: Iterable[Iterable[Fraction]], right: Iterable[Fraction]
) -> tuple[Fraction, ...]:
    rows = [list(map(Fraction, row)) + [Fraction(value)] for row, value in zip(matrix, right)]
    size = len(rows)
    if size == 0 or any(len(row) != size + 1 for row in rows):
        raise ValueError("linear system must be square")
    for column in range(size):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        if pivot is None:
            raise ValueError("singular lower-facet interpolation system")
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
        pivot_value = rows[column][column]
        rows[column] = [value / pivot_value for value in rows[column]]
        for row in range(size):
            if row == column or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(rows[row], rows[column])
            ]
    return tuple(row[-1] for row in rows)


def lower_plane(
    facet: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
    heights: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    matrix = [[Fraction(1), *anchors[index][1:]] for index in facet]
    return solve_linear(matrix, [heights[index] for index in facet])


def plane_value(coefficients: tuple[Fraction, ...], profile: tuple[Fraction, ...]) -> Fraction:
    return coefficients[0] + sum(
        coefficients[index] * profile[index] for index in range(1, 5)
    )


def exact_lower_facet(
    facet: tuple[int, ...],
    anchors: tuple[tuple[Fraction, ...], ...],
    heights: tuple[Fraction, ...],
) -> tuple[Fraction, int]:
    determinant = verifier.simplex_determinant(facet, anchors)
    if not determinant:
        raise ValueError(f"Qhull proposed exactly degenerate projected facet {facet}")
    coefficients = lower_plane(facet, anchors, heights)
    minimum_slack = None
    for index, (profile, height) in enumerate(zip(anchors, heights)):
        slack = height - plane_value(coefficients, profile)
        if index in facet:
            if slack:
                raise ValueError(f"facet interpolation failed exactly for {facet}")
            continue
        if slack <= 0:
            relation = "coplanar" if slack == 0 else "below"
            raise ValueError(
                f"proposed lower facet {facet} has nonvertex anchor {index} {relation} its plane"
            )
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
    assert minimum_slack is not None
    return minimum_slack, 1 if determinant > 0 else -1


def initialize_audit_worker(
    anchors: tuple[tuple[Fraction, ...], ...],
    heights: tuple[Fraction, ...],
) -> None:
    """Install immutable audit data once per process, not once per facet."""

    global _WORKER_ANCHORS, _WORKER_HEIGHTS
    _WORKER_ANCHORS = anchors
    _WORKER_HEIGHTS = heights


def audit_facet_chunk(
    task: tuple[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, list[tuple[Fraction, int]]]:
    chunk_index, facets = task
    if _WORKER_ANCHORS is None or _WORKER_HEIGHTS is None:
        raise RuntimeError("regular-mesh audit worker was not initialized")
    return chunk_index, [
        exact_lower_facet(facet, _WORKER_ANCHORS, _WORKER_HEIGHTS)
        for facet in facets
    ]


def audit_lower_facets(
    facets: list[tuple[int, ...]],
    anchors: tuple[tuple[Fraction, ...], ...],
    heights: tuple[Fraction, ...],
    *,
    workers: int,
    chunk_size: int,
    progress_every: int,
) -> list[tuple[Fraction, int]]:
    """Audit independent facets in parallel while retaining facet order."""

    if workers == 1:
        result = []
        for completed, facet in enumerate(facets, 1):
            result.append(exact_lower_facet(facet, anchors, heights))
            if progress_every and completed % progress_every == 0:
                print(f"exact_facets_audited={completed}/{len(facets)}", flush=True)
        return result

    chunks = [
        (index, tuple(facets[start : start + chunk_size]))
        for index, start in enumerate(range(0, len(facets), chunk_size))
    ]
    ordered: list[list[tuple[Fraction, int]] | None] = [None] * len(chunks)
    completed = 0
    next_progress = progress_every
    with ProcessPoolExecutor(
        max_workers=min(workers, len(chunks)),
        initializer=initialize_audit_worker,
        initargs=(anchors, heights),
    ) as executor:
        futures = [executor.submit(audit_facet_chunk, task) for task in chunks]
        for future in as_completed(futures):
            chunk_index, rows = future.result()
            ordered[chunk_index] = rows
            completed += len(rows)
            if progress_every and completed >= next_progress:
                print(f"exact_facets_audited={completed}/{len(facets)}", flush=True)
                next_progress = ((completed // progress_every) + 1) * progress_every
    if any(rows is None for rows in ordered):
        raise RuntimeError("regular-mesh parallel audit lost a facet chunk")
    return [row for rows in ordered for row in rows or ()]


def propose_lower_facets(
    anchors: tuple[tuple[Fraction, ...], ...], heights: tuple[Fraction, ...]
) -> list[tuple[int, ...]]:
    points = np.asarray(
        [
            [*(float(value / verifier.M) for value in profile[1:]), float(height)]
            for profile, height in zip(anchors, heights)
        ],
        dtype=np.float64,
    )
    try:
        hull = ConvexHull(points, qhull_options="Qx Qt Qc")
    except QhullError as error:
        raise ValueError(f"Qhull could not propose a lifted convex hull: {error}") from error
    facets = {
        tuple(sorted(int(index) for index in simplex))
        for simplex, equation in zip(hull.simplices, hull.equations)
        if float(equation[-2]) < 0.0
    }
    if not facets:
        raise ValueError("Qhull proposed no lower facets")
    return sorted(facets)


def compact_hash(value) -> str:
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def discard_projected_degenerate_facets(
    facets: list[tuple[int, ...]],
    anchors: tuple[tuple[Fraction, ...], ...],
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Remove only proposals with an exactly zero projected determinant.

    These are artifacts of Qhull's lifted triangulation, not projected
    four-simplices.  Every retained proposal still undergoes the strict exact
    lower-support audit; no other failure is downgraded or discarded.
    """

    retained = []
    discarded = []
    for facet in facets:
        if verifier.simplex_determinant(facet, anchors) == 0:
            discarded.append(facet)
        else:
            retained.append(facet)
    if not retained:
        raise ValueError("all proposed lower facets are projected-degenerate")
    return retained, discarded


def build(
    ledger: dict,
    *,
    workers: int = 1,
    chunk_size: int = 16,
    progress_every: int = 0,
) -> dict:
    if int(ledger.get("group_bits", verifier.GROUP_BITS)) != verifier.GROUP_BITS:
        raise ValueError("regular mesh requires group_bits=4")
    anchors = tuple(sorted(verifier.load_anchors(ledger)))
    heights = tuple(lifting_height(profile, rank) for rank, profile in enumerate(anchors))
    proposed_facets = propose_lower_facets(anchors, heights)
    facets, degenerate_proposals = discard_projected_degenerate_facets(
        proposed_facets, anchors
    )
    audited = audit_lower_facets(
        facets,
        anchors,
        heights,
        workers=workers,
        chunk_size=chunk_size,
        progress_every=progress_every,
    )
    rows = []
    minimum_slack = None
    orientations = {1: 0, -1: 0}
    for cell_id, (facet, (slack, orientation)) in enumerate(zip(facets, audited)):
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
        orientations[orientation] += 1
        rows.append({"cell_id": f"r{cell_id:06d}", "anchor_indices": list(facet)})
    geometry = verifier.verify_mesh(
        {"complete_cover": True, "simplices": len(rows)}, anchors, rows
    )
    height_strings = [str(value) for value in heights]
    return {
        "schema": SCHEMA,
        "status": "EXACT_CERTIFIED_G4_REGULAR_ROOT_TRIANGULATION",
        "group_bits": verifier.GROUP_BITS,
        "complete_exact_regular_mesh": True,
        "anchors": len(anchors),
        "simplices": len(rows),
        "exact_degenerate_simplices": 0,
        "qhull_lower_facet_proposals": len(proposed_facets),
        "exact_projected_degenerate_proposals_discarded": len(
            degenerate_proposals
        ),
        "exact_projected_degenerate_proposal_list_sha256": compact_hash(
            [list(facet) for facet in degenerate_proposals]
        ),
        "anchor_profiles": [[str(value) for value in profile] for profile in anchors],
        "root_cell_ledger": rows,
        "lifting": {
            "base": "sum_{i=1}^4 (a_i/M)^2",
            "perturbation": "((rank+1)^5 mod 1000003)/(1000003*2^20)",
            "rank_order": "lexicographic exact anchor profile order",
            "heights_sha256": compact_hash(height_strings),
            "minimum_exact_nonvertex_lower_slack": str(minimum_slack),
            "positive_projected_orientations": orientations[1],
            "negative_projected_orientations": orientations[-1],
        },
        "exact_geometry": {"mode": "exact_regular_triangulation", **geometry},
        "source_anchor_profiles_sha256": compact_hash(ledger.get("anchor_profiles")),
    }


def run_self_test() -> None:
    # Eight exact clipped-hull vertices exercise Qhull proposal, exact support
    # audit, and the independent topology/volume verifier cheaply.
    anchors = []
    for weight in range(1, 5):
        boundary = [Fraction(0)] * 5
        boundary[0] = verifier.M - Fraction(verifier.MINIMUM_PHYSICAL_WEIGHT, weight)
        boundary[weight] = Fraction(verifier.MINIMUM_PHYSICAL_WEIGHT, weight)
        anchors.append(tuple(boundary))
        pure = [Fraction(0)] * 5
        pure[weight] = verifier.M
        anchors.append(tuple(pure))
    ledger = {
        "group_bits": verifier.GROUP_BITS,
        "anchor_profiles": [[str(value) for value in row] for row in anchors],
    }
    serial = build(ledger, workers=1, chunk_size=2)
    parallel = build(ledger, workers=2, chunk_size=2)
    if serial != parallel:
        raise SystemExit("regular-mesh self-test: serial/parallel reports differ")
    synthetic_anchors = tuple(sorted(verifier.load_anchors(ledger)))
    midpoint = tuple(
        (synthetic_anchors[0][coordinate] + synthetic_anchors[1][coordinate]) / 2
        for coordinate in range(5)
    )
    augmented = (*synthetic_anchors, midpoint)
    degenerate = tuple(sorted((0, 1, 2, 3, len(augmented) - 1)))
    nondegenerate = tuple(serial["root_cell_ledger"][0]["anchor_indices"])
    retained, discarded = discard_projected_degenerate_facets(
        sorted((degenerate, nondegenerate)), augmented
    )
    if retained != [nondegenerate] or discarded != [degenerate]:
        raise SystemExit("regular-mesh self-test: exact degeneracy filter mismatch")
    print(f"self_test_anchors={serial['anchors']}")
    print(f"self_test_simplices={serial['simplices']}")
    print("status=PASS_G4_EXACT_REGULAR_MESH_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--facet-chunk-size", type=int, default=16)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.input is None or args.output is None:
        raise SystemExit("regular mesh: --input and --output are required")
    if args.workers <= 0 or args.facet_chunk_size <= 0 or args.progress_every < 0:
        raise SystemExit("regular mesh: invalid parallel audit parameter")
    ledger = json.loads(args.input.read_text(encoding="utf-8"))
    report = build(
        ledger,
        workers=args.workers,
        chunk_size=args.facet_chunk_size,
        progress_every=args.progress_every,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"anchors={report['anchors']} simplices={report['simplices']}")
    print("exact_degenerate_simplices=0")
    print(
        "exact_projected_degenerate_proposals_discarded="
        f"{report['exact_projected_degenerate_proposals_discarded']}"
    )
    print("exact_lower_facets=PASS exact_mesh_audit=PASS")


if __name__ == "__main__":
    main()
