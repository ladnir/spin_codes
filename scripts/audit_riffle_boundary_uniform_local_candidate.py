#!/usr/bin/env python3
"""Audit a factor-four local bound on all remaining exact boundary data."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_riffle_small_boundary_saddle import (
    boundary_transition,
    lattice_index,
    log_mass_moments,
)
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS


OPEN_FACES = (
    (0, 1, 2, 3),
    (0, 1, 3, 4),
    (0, 1, 2, 3, 4),
)


def is_closed_face(active: list[int]) -> bool:
    weights = set(active)
    return weights.issubset({0, 2, 4}) or weights.issubset({0, 1, 3})


def logsumexp2(values: list[float]) -> float:
    natural = np.asarray(values, dtype=np.float64) * math.log(2.0)
    return float(np.logaddexp.reduce(natural)) / math.log(2.0)


def active_states(face: tuple[int, ...]) -> tuple[int, ...]:
    forward = {0}
    changed = True
    while changed:
        changed = False
        for old in tuple(forward):
            for packet in face:
                for new in range(5):
                    if boundary_transition(old, new, packet) and new not in forward:
                        forward.add(new)
                        changed = True

    reverse = {0}
    changed = True
    while changed:
        changed = False
        for new in tuple(reverse):
            for packet in face:
                for old in range(5):
                    if boundary_transition(old, new, packet) and old not in reverse:
                        reverse.add(old)
                        changed = True
    return tuple(sorted(forward & reverse))


def primitivity_exponent(face: tuple[int, ...], states: tuple[int, ...]) -> int:
    adjacency = np.zeros((len(states), len(states)), dtype=bool)
    index = {state: i for i, state in enumerate(states)}
    for old in states:
        for new in states:
            adjacency[index[old], index[new]] = any(
                boundary_transition(old, new, packet) for packet in face
            )
    power = adjacency.copy()
    for exponent in range(1, 65):
        if bool(np.all(power)):
            return exponent
        power = (power.astype(np.uint8) @ adjacency.astype(np.uint8)) > 0
    raise RuntimeError("the active boundary graph is not primitive")


def face_audit(face: tuple[int, ...]) -> dict[str, object]:
    states = active_states(face)
    dimension = len(face) - 1
    _log_mass, _gradient, covariance = log_mass_moments(
        face, np.zeros(dimension, dtype=np.float64), 64
    )
    eigenvalues = np.linalg.eigvalsh(covariance)
    lattice_indices = {lattice_index(face, length) for length in range(6, 21)}
    return {
        "active_packet_weights": list(face),
        "dimension": dimension,
        "active_state_weights": list(states),
        "primitivity_exponent": primitivity_exponent(face, states),
        "lattice_indices_for_lengths_6_through_20": sorted(lattice_indices),
        "zero_tilt_covariance_minimum_eigenvalue_at_length_64": float(
            eigenvalues[0]
        ),
        "zero_tilt_covariance_maximum_eigenvalue_at_length_64": float(
            eigenvalues[-1]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor", type=float, default=4.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.factor >= 1.0:
        raise ValueError("the candidate factor must be at least one")
    safety_bits = math.log2(args.factor)

    comparisons = []
    shell_rows = []
    for data_blocks in (4, 8, 12, 15):
        path = DEFAULT_RECEIPTS / (
            f"goal22_boundary_saddle_b{data_blocks}_p2_h16_d8.json"
        )
        receipt = json.loads(path.read_text(encoding="utf-8"))
        candidate_terms = []
        for row in receipt["type_rows"]:
            if is_closed_face(row["active_packet_weights"]):
                candidate_terms.append(row["exact_expected_contribution_log2"])
                continue
            error = row["saddle_approximation_error_bits"]
            comparisons.append(
                {
                    "source": path.name,
                    "histogram_h0_through_h4": row[
                        "histogram_h0_through_h4"
                    ],
                    "dimension": row["saddle_dimension"],
                    "gaussian_error_bits": error,
                    "candidate_slack_bits": error + safety_bits,
                }
            )
            candidate_terms.append(
                row["saddle_expected_contribution_log2"] + safety_bits
            )
        candidate_log2 = logsumexp2(candidate_terms)
        exact_log2 = receipt["exact_shell_expected_count_log2"]
        shell_rows.append(
            {
                "data_blocks": data_blocks,
                "binary_output_length": receipt["binary_output_length"],
                "exact_shell_log2": exact_log2,
                "candidate_shell_bound_log2": candidate_log2,
                "candidate_shell_loss_bits": candidate_log2 - exact_log2,
            }
        )

    scaling_path = DEFAULT_RECEIPTS / "goal22_boundary_saddle_scaling_s4.json"
    scaling = json.loads(scaling_path.read_text(encoding="utf-8"))
    for row in scaling["rows"]:
        if row["saddle_dimension"] < 3:
            continue
        error = row["saddle_approximation_error_bits"]
        comparisons.append(
            {
                "source": scaling_path.name,
                "profile": row["profile"],
                "scale": row["scale"],
                "dimension": row["saddle_dimension"],
                "gaussian_error_bits": error,
                "candidate_slack_bits": error + safety_bits,
            }
        )

    worst = min(comparisons, key=lambda row: row["candidate_slack_bits"])
    payload = {
        "schema": "riffle-boundary-uniform-local-candidate-v1",
        "evidence_label": "EXACT_DATA_AUDIT_OF_UNPROVED_LOCAL_BOUND",
        "candidate_factor": args.factor,
        "candidate_safety_bits": safety_bits,
        "exact_comparisons": len(comparisons),
        "minimum_candidate_slack_bits": worst["candidate_slack_bits"],
        "maximum_required_safety_bits": max(
            0.0,
            max(-row["gaussian_error_bits"] for row in comparisons),
        ),
        "worst_case": worst,
        "candidate_passes_exact_data": all(
            row["candidate_slack_bits"] >= 0.0 for row in comparisons
        ),
        "open_face_audits": [face_audit(face) for face in OPEN_FACES],
        "shell_rows": shell_rows,
        "comparisons": comparisons,
        "scope": (
            "The exact coefficients come from Goal 22. Passing this audit does "
            "not prove the candidate inequality. Graph reachability and lattice "
            "indices are exact; covariance eigenvalues are numerical."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS / "goal25_boundary_uniform_local_candidate.json"
    )
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "candidate_factor": args.factor,
                "exact_comparisons": len(comparisons),
                "passes": payload["candidate_passes_exact_data"],
                "minimum_slack_bits": payload["minimum_candidate_slack_bits"],
                "required_safety_bits": payload["maximum_required_safety_bits"],
                "shell_losses_bits": [
                    row["candidate_shell_loss_bits"] for row in shell_rows
                ],
                "status": (
                    "PASS" if payload["candidate_passes_exact_data"] else "FAIL"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
