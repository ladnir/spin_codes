#!/usr/bin/env python3
"""Test the boundary saddle correction under proportional histogram scaling."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from probe_riffle_small_boundary_saddle import (
    boundary_transition,
    lattice_index,
    optimize_face,
)
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS


PROFILES = {
    "even_line": (4, 0, 8, 0, 0),
    "three_dimensional": (4, 2, 4, 2, 0),
    "four_dimensional": (5, 1, 4, 1, 1),
}


def exact_boundary_coefficient(histogram: tuple[int, ...]) -> int:
    """Extract [x^h] e0^T L(x)^N e0 with integer arithmetic."""
    active = tuple(packet for packet, count in enumerate(histogram) if count)
    target = tuple(histogram[packet] for packet in active)
    zero = (0,) * len(active)
    current: dict[tuple[tuple[int, ...], int], int] = {(zero, 0): 1}
    for _position in range(sum(histogram)):
        following: defaultdict[tuple[tuple[int, ...], int], int] = defaultdict(int)
        for (counts, old), ways in current.items():
            for coordinate, packet in enumerate(active):
                if counts[coordinate] == target[coordinate]:
                    continue
                next_counts = list(counts)
                next_counts[coordinate] += 1
                next_counts_tuple = tuple(next_counts)
                for new in range(5):
                    multiplicity = boundary_transition(old, new, packet)
                    if multiplicity:
                        following[(next_counts_tuple, new)] += ways * multiplicity
        current = dict(following)
    return current.get((target, 0), 0)


def saddle_row(name: str, histogram: tuple[int, ...], scale: int) -> dict[str, object]:
    packet_positions = sum(histogram)
    exact = exact_boundary_coefficient(histogram)
    if not exact:
        raise RuntimeError("the scaled boundary coefficient is zero")
    bound_log, parameters, covariance, success = optimize_face(
        histogram, packet_positions
    )
    sign, log_determinant = np.linalg.slogdet(covariance)
    if not sign > 0.0:
        raise RuntimeError("the saddle covariance is not positive definite")
    active = tuple(packet for packet, count in enumerate(histogram) if count)
    dimension = len(active) - 1
    index = lattice_index(active, packet_positions)
    gaussian_log_probability = (
        math.log(index)
        - 0.5 * dimension * math.log(2.0 * math.pi)
        - 0.5 * float(log_determinant)
    )
    raw_loss_bits = bound_log / math.log(2.0) - math.log2(exact)
    return {
        "profile": name,
        "scale": scale,
        "histogram_h0_through_h4": list(histogram),
        "packet_positions": packet_positions,
        "input_weight": sum(
            packet * count for packet, count in enumerate(histogram)
        ),
        "output_weight": sum(
            packet * count for packet, count in enumerate(histogram)
        )
        // 2,
        "saddle_dimension": dimension,
        "lattice_index": index,
        "exact_coefficient": str(exact),
        "coefficient_bound_loss_bits": raw_loss_bits,
        "gaussian_log_probability_bits": gaussian_log_probability / math.log(2.0),
        "saddle_approximation_error_bits": (
            raw_loss_bits + gaussian_log_probability / math.log(2.0)
        ),
        "parameters": [float(value) for value in parameters],
        "covariance_log_determinant": float(log_determinant),
        "optimizer_success": success,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-scale", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.maximum_scale < 1:
        raise ValueError("maximum scale must be positive")

    rows = []
    for scale in range(1, args.maximum_scale + 1):
        for name, base in PROFILES.items():
            histogram = tuple(scale * count for count in base)
            rows.append(saddle_row(name, histogram, scale))

    payload = {
        "schema": "riffle-boundary-saddle-scaling-v1",
        "evidence_label": "EXACT_INTEGER_COEFFICIENTS_AND_NUMERICAL_SADDLES",
        "maximum_scale": args.maximum_scale,
        "profiles": {name: list(histogram) for name, histogram in PROFILES.items()},
        "optimizer_failures": sum(int(not row["optimizer_success"]) for row in rows),
        "rows": rows,
        "scope": (
            "Every boundary coefficient is exact. Saddle points, covariance "
            "determinants, and Gaussian local factors are numerical "
            "approximations and are not certified bounds."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / f"goal22_boundary_saddle_scaling_s{args.maximum_scale}.json"
    )
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {}
    for name in PROFILES:
        selected = [row for row in rows if row["profile"] == name]
        summary[name] = [
            {
                "scale": row["scale"],
                "raw_loss_bits": row["coefficient_bound_loss_bits"],
                "saddle_error_bits": row["saddle_approximation_error_bits"],
            }
            for row in selected
        ]
    print(
        json.dumps(
            {
                "output": str(output),
                "optimizer_failures": payload["optimizer_failures"],
                "profiles": summary,
                "status": "PASS" if not payload["optimizer_failures"] else "CHECK",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
