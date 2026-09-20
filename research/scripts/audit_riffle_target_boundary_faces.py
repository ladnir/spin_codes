#!/usr/bin/env python3
"""Audit every packet-support face that can occur at target boundary weight."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

from audit_riffle_boundary_uniform_local_candidate import (
    active_states,
    primitivity_exponent,
)
from probe_riffle_small_boundary_saddle import lattice_index


EXACT_TWO_DIMENSIONAL_FACES = {(0, 1, 3), (0, 2, 4)}
GOAL25_OPEN_FACES = {
    (0, 1, 2, 3),
    (0, 1, 3, 4),
    (0, 1, 2, 3, 4),
}


def status(face: tuple[int, ...]) -> str:
    if len(face) == 2 or face in EXACT_TWO_DIMENSIONAL_FACES:
        return "EXACT_FORMULA_GOALS_23_24"
    if face in GOAL25_OPEN_FACES:
        return "AUDITED_AS_OPEN_FACE_IN_GOAL_25"
    return "NEW_TARGET_OPEN_FACE"


def main() -> None:
    rows = []
    for nonzero_count in range(1, 5):
        for subset in itertools.combinations(range(1, 5), nonzero_count):
            face = (0,) + subset
            states = active_states(face)
            rows.append(
                {
                    "active_packet_weights": list(face),
                    "dimension": len(face) - 1,
                    "active_state_weights": list(states),
                    "primitivity_exponent": primitivity_exponent(face, states),
                    "lattice_indices_lengths_6_through_20": [
                        lattice_index(face, length) for length in range(6, 21)
                    ],
                    "prior_status": status(face),
                }
            )

    payload = {
        "schema": "riffle-target-boundary-face-audit-v1",
        "evidence_label": "EXACT_GRAPH_AND_LATTICE_AUDIT",
        "face_count": len(rows),
        "exact_formula_faces": sum(
            row["prior_status"] == "EXACT_FORMULA_GOALS_23_24" for row in rows
        ),
        "goal25_open_faces": sum(
            row["prior_status"] == "AUDITED_AS_OPEN_FACE_IN_GOAL_25" for row in rows
        ),
        "new_target_open_faces": sum(
            row["prior_status"] == "NEW_TARGET_OPEN_FACE" for row in rows
        ),
        "all_faces_primitive": all(row["primitivity_exponent"] for row in rows),
        "rows": rows,
        "scope": (
            "The target boundary permits every nonempty subset of packet "
            "weights one through four. Graph reachability, primitivity, and "
            "finite-length lattice indices use exact integer arithmetic."
        ),
    }
    output = Path(
        "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/"
        "receipts/goal26_target_boundary_faces.json"
    )
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(output),
                "face_count": payload["face_count"],
                "exact_formula_faces": payload["exact_formula_faces"],
                "goal25_open_faces": payload["goal25_open_faces"],
                "new_target_open_faces": payload["new_target_open_faces"],
                "all_faces_primitive": payload["all_faces_primitive"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
