#!/usr/bin/env python3
"""Split the small accumulator pair chain into scalar and interaction modes.

For a rank-two pair (x,y), the three nonzero binary characters are x, y,
and x+y.  Functions of their three Hamming weights span the scalar-character
sector.  That sector is invariant under the common-permutation accumulator
chain.  This diagnostic measures the remaining interaction block at small
lengths.  It does not extrapolate the measurements to length 512.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import qr

from analyze_accumulator_pair_chain_small import compositions4, rank_two
from probe_triangle_holder_pair_bound_small import triple_weights
from verify_accumulator_pair_type_kernel import kernel_counts, multinomial


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "accumulator_pair_interaction_B8_probe.json"


def transition_matrix(length: int) -> tuple[list[tuple[int, int, int, int]], np.ndarray, np.ndarray]:
    types = [kind for kind in compositions4(length) if rank_two(kind)]
    index = {kind: position for position, kind in enumerate(types)}
    transition = np.zeros((len(types), len(types)), dtype=np.float64)
    for row, input_type in enumerate(types):
        denominator = multinomial(input_type)
        for output_type, count in kernel_counts(input_type).items():
            transition[row, index[output_type]] = count / denominator
    stationary_counts = np.asarray([multinomial(kind) for kind in types], dtype=float)
    stationary = stationary_counts / stationary_counts.sum()
    return types, transition, stationary


def scalar_basis(
    types: list[tuple[int, int, int, int]], stationary: np.ndarray, length: int
) -> tuple[np.ndarray, int]:
    lifts = np.zeros((len(types), 3 * length), dtype=float)
    for row, kind in enumerate(types):
        weights = triple_weights(kind)
        for character, weight in enumerate(weights):
            lifts[row, character * length + weight - 1] = np.sqrt(stationary[row])
    raw_q, raw_r, _ = qr(lifts, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(raw_r))
    tolerance = max(lifts.shape) * np.finfo(float).eps * diagonal[0]
    rank = int(np.sum(diagonal > tolerance))
    return raw_q[:, :rank], rank


def spectral_norm(matrix: np.ndarray) -> float:
    return float(np.linalg.svd(matrix, compute_uv=False)[0]) if matrix.size else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, choices=range(2, 21), default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    types, transition, stationary = transition_matrix(args.length)
    root = np.sqrt(stationary)
    weighted = root[:, None] * transition / root[None, :]
    basis, scalar_rank = scalar_basis(types, stationary, args.length)
    projector = basis @ basis.T
    identity = np.eye(len(types))
    complement = identity - projector

    scalar_block = basis.T @ weighted @ basis
    invariant_residual = complement @ weighted @ basis
    interaction_to_scalar = basis.T @ weighted @ complement
    interaction_block = complement @ weighted @ complement

    full_values = np.linalg.svd(weighted, compute_uv=False)
    scalar_values = np.linalg.svd(scalar_block, compute_uv=False)
    interaction_values = np.linalg.svd(interaction_block, compute_uv=False)
    payload = {
        "schema": "accumulator-pair-interaction-small-probe-v1",
        "status": "EXACT_KERNEL_WITH_BINARY64_LINEAR_ALGEBRA_DIAGNOSTIC",
        "parameters": {
            "length": args.length,
            "rank_two_pair_types": len(types),
            "scalar_character_span_rank": scalar_rank,
            "expected_scalar_character_span_rank": 3 * args.length - 2,
        },
        "checks": {
            "maximum_stationarity_error": float(
                np.max(np.abs(stationary @ transition - stationary))
            ),
            "scalar_invariance_spectral_residual": spectral_norm(invariant_residual),
        },
        "blocks": {
            "interaction_to_scalar_spectral_norm": spectral_norm(interaction_to_scalar),
            "interaction_block_spectral_norm": float(interaction_values[0]),
            "leading_full_singular_values": full_values[:12].tolist(),
            "leading_scalar_block_singular_values": scalar_values[:12].tolist(),
            "leading_interaction_block_singular_values": interaction_values[:12].tolist(),
        },
        "interpretation": [
            "The scalar-character sector contains every function of wt(x), wt(y), or wt(x+y).",
            "Its invariance is an exact probabilistic identity; only the displayed residual is numerical.",
            "The interaction block retains the shared-overlap information lost by three-marginal bounds.",
        ],
        "limitations": [
            "Kernel entries are exact integer ratios, but all linear algebra uses binary64.",
            "The interaction norm is a worst-function diagnostic, not the target shell-specific variance bound.",
            "Small-length values do not prove a length-512 bound.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "length": args.length,
                "types": len(types),
                "scalar_rank": scalar_rank,
                "invariance_residual": payload["checks"]["scalar_invariance_spectral_residual"],
                "interaction_to_scalar": payload["blocks"]["interaction_to_scalar_spectral_norm"],
                "interaction_norm": payload["blocks"]["interaction_block_spectral_norm"],
            },
            indent=2,
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
