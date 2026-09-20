#!/usr/bin/env python3
"""Fit the narrow finite-SPIN parameter extrapolation used by the paper.

The fit is deliberately diagnostic.  It uses only occupation-one, binary64
results for exact-spectrum BCH constituents under the matched-persistence
RM2Sub schedule.  It does not turn those rows into a distance certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
from contextlib import closing
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_DB = HERE / "spin_landscape.sqlite3"
DEFAULT_OUTPUT = HERE / "parameter_extrapolation.json"


def solve_3x3(matrix: list[list[float]], vector: list[float]) -> list[float]:
    augmented = [row[:] + [rhs] for row, rhs in zip(matrix, vector)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-12:
            raise ArithmeticError("singular normal equations")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(3):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][3] for row in range(3)]


def least_squares(rows: list[tuple[float, float, float]]) -> tuple[list[float], float, float]:
    # margin = intercept + distance_slope * d + exponent_slope * (e - 20)
    normal = [[0.0] * 3 for _ in range(3)]
    rhs = [0.0] * 3
    for distance, exponent, margin in rows:
        features = (1.0, distance, exponent - 20.0)
        for i in range(3):
            rhs[i] += features[i] * margin
            for j in range(3):
                normal[i][j] += features[i] * features[j]
    coefficients = solve_3x3(normal, rhs)
    residuals = []
    for distance, exponent, margin in rows:
        predicted = (
            coefficients[0]
            + coefficients[1] * distance
            + coefficients[2] * (exponent - 20.0)
        )
        residuals.append(margin - predicted)
    rmse = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    max_error = max(abs(value) for value in residuals)
    return coefficients, rmse, max_error


def predict(coefficients: list[float], distance: float, exponent: float) -> float:
    return coefficients[0] + coefficients[1] * distance + coefficients[2] * (exponent - 20.0)


def grouped_cross_validation(
    rows: list[tuple[int, float, float, float]], group_index: int
) -> dict[str, object]:
    predictions = []
    group_value = lambda row: (row[0], row[2])[group_index]
    groups = sorted({group_value(row) for row in rows})
    for group in groups:
        training = [
            (distance, exponent, margin)
            for block, distance, exponent, margin in rows
            if (block, exponent)[group_index] != group
        ]
        coefficients, _, _ = least_squares(training)
        for block, distance, exponent, margin in rows:
            if (block, exponent)[group_index] == group:
                estimate = predict(coefficients, distance, exponent)
                predictions.append(
                    {
                        "held_out_group": group,
                        "block_bits": block,
                        "minimum_distance": distance,
                        "message_exponent": exponent,
                        "observed_margin_bits": margin,
                        "predicted_margin_bits": estimate,
                        "error_bits": estimate - margin,
                    }
                )
    errors = [row["error_bits"] for row in predictions]
    return {
        "groups": groups,
        "rmse_bits": math.sqrt(sum(error * error for error in errors) / len(errors)),
        "maximum_absolute_error_bits": max(abs(error) for error in errors),
        "predictions": predictions,
    }


def fit_distance_constant() -> tuple[float, list[dict[str, float]]]:
    # The B=256 point supplies only its known minimum distance, not a spectrum.
    points = [(32, 8), (64, 12), (128, 22), (256, 38)]
    features = [block / math.log2(block) for block, _ in points]
    coefficient = sum(x * distance for x, (_, distance) in zip(features, points)) / sum(
        x * x for x in features
    )
    return coefficient, [
        {
            "block_bits": block,
            "minimum_distance": distance,
            "normalized_distance_constant": distance * math.log2(block) / block,
        }
        for block, distance in points
    ]


def invert_distance_model(required_distance: float, coefficient: float) -> float:
    low, high = 2.0, 4096.0
    while coefficient * high / math.log2(high) < required_distance:
        high *= 2.0
    for _ in range(100):
        middle = (low + high) / 2.0
        if coefficient * middle / math.log2(middle) < required_distance:
            low = middle
        else:
            high = middle
    return high


def rm2sub_capacity(m: int) -> int:
    return 1 + m + m * (m - 1) // 2


def rm2sub_schedule(message_exponent: int, offset: int) -> dict[str, int]:
    persistence = message_exponent + offset
    # The landscape has audited RM2Sub maps only from t=64 onward.
    for m in range(6, 32):
        state = persistence - m
        if 1 <= state <= rm2sub_capacity(m):
            return {
                "message_exponent": message_exponent,
                "persistence_offset": offset,
                "step_exponent": m,
                "step_bits": 1 << m,
                "state_bits": state,
                "rm2_capacity": rm2sub_capacity(m),
            }
    raise ArithmeticError("no RM2Sub schedule found")


def build(db_path: Path) -> dict[str, object]:
    with closing(sqlite3.connect(db_path)) as db:
        source_rows = db.execute(
            """
            SELECT result_id, block_bits, minimum_distance, message_exponent, margin_bits
            FROM q1_curves
            WHERE study = 'rm2sub_t64_matched_persistence_family_curve'
              AND outer_family = 'bch'
              AND block_bits >= 32
              AND message_exponent BETWEEN 16 AND 20
            ORDER BY block_bits, message_exponent
            """
        ).fetchall()

    fit_rows = [
        (float(distance), float(exponent), float(margin))
        for _, _, distance, exponent, margin in source_rows
    ]
    cross_validation_rows = [
        (int(block), float(distance), float(exponent), float(margin))
        for _, block, distance, exponent, margin in source_rows
    ]
    coefficients, rmse, max_error = least_squares(fit_rows)
    distance_coefficient, distance_points = fit_distance_constant()
    target_margin = 40.0
    projections = []
    for exponent in (16, 20, 24, 28, 32, 36, 40):
        required_distance = (
            target_margin - coefficients[0] - coefficients[2] * (exponent - 20.0)
        ) / coefficients[1]
        continuous_block = invert_distance_model(required_distance, distance_coefficient)
        power_two_block = 1 << math.ceil(math.log2(continuous_block))
        projections.append(
            {
                "message_exponent": exponent,
                "target_margin_bits": target_margin,
                "projected_required_minimum_distance": required_distance,
                "projected_continuous_block_bits": continuous_block,
                "next_power_of_two_block_bits": power_two_block,
            }
        )

    return {
        "schema": "spin-parameter-extrapolation-v1",
        "status": "historical Q1 fit; activation-state transfer pending review; not for current parameter selection",
        "distance_target": 0.10,
        "source_database": {
            "path": db_path.name,
            "sha256": hashlib.sha256(db_path.read_bytes()).hexdigest(),
        },
        "source_results": [
            {
                "result_id": result_id,
                "block_bits": block,
                "minimum_distance": distance,
                "message_exponent": exponent,
                "margin_bits": margin,
            }
            for result_id, block, distance, exponent, margin in source_rows
        ],
        "fit": {
            "formula": "margin = intercept + distance_slope*d + exponent_slope*(log2(k)-20)",
            "scope": (
                "Q=1 binary64 rows; exact-spectrum BCH B=32,64,128; "
                "RM2Sub t=64; p=log2(k)+2; log2(k)=16,...,20"
            ),
            "observation_count": len(fit_rows),
            "intercept": coefficients[0],
            "distance_slope": coefficients[1],
            "exponent_slope": coefficients[2],
            "rmse_bits": rmse,
            "maximum_absolute_residual_bits": max_error,
        },
        "cross_validation": {
            "leave_one_message_exponent_out": grouped_cross_validation(
                cross_validation_rows, 1
            ),
            "leave_one_block_size_out": grouped_cross_validation(cross_validation_rows, 0),
        },
        "bch_distance_model": {
            "formula": "d(B) = coefficient*B/log2(B)",
            "coefficient": distance_coefficient,
            "points": distance_points,
            "caveat": (
                "The B=256 point is the deterministic BCH-derived [256,128,>=38] "
                "subcode, not a pure BCH constituent; the model uses only its known "
                "distance, not its spectrum."
            ),
        },
        "forty_bit_q1_projections": projections,
        "rm2sub_staircases": {
            "aggressive_offset_2": [rm2sub_schedule(exponent, 2) for exponent in range(16, 41, 4)],
            "certificate_anchored_offset_4": [
                rm2sub_schedule(exponent, 4) for exponent in range(16, 41, 4)
            ],
        },
        "limitations": [
            "Source Q1 screens use the historical two-state activation invariant; retain this fit only as a comparison until recomputed.",
            "The regression extrapolates occupation-one diagnostics only.",
            "The fit does not use a BCH [256,128] spectrum or RM2Sub transfer receipt.",
            "The BCH distance law is a size heuristic and is not a spectrum theorem.",
            "Full-distance selection still requires every occupation and outward rounding.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build(args.db)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["fit"], indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
