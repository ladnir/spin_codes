#!/usr/bin/env python3
"""Sum all outer occupations for RandomStepConv Goal 04."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from analyze_riffle_randomstepconv_g4_sigma20_goal02 import (
    TARGET_G,
    TARGET_PACKETS,
    TARGET_SIGMA,
    log_binomial,
    log_common_bound,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal03 import (
    DATA_POSITIONS,
    DEFAULT_INNER_RECEIPT,
    DEFAULT_SPECTRUM,
    LocalMomentModel,
    load_inner_caps,
    load_spectrum,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_randomstepconv_g4_sigma20/receipts/"
    "goal04_all_occupation_closure.json"
)
LOG2 = math.log(2.0)


def log_binomial_tail_from_two(
    log_x: float, positions: int = DATA_POSITIONS
) -> float:
    """Return log(sum_{r=2}^K C(K,r) x^r) stably."""
    expected = positions * math.exp(log_x) if log_x < 700.0 else math.inf
    if expected < 0.1:
        x = math.exp(log_x)
        log_term_two = (
            math.log(positions)
            + math.log(positions - 1)
            - math.log(2.0)
            + 2.0 * log_x
        )
        relative_term = 1.0
        relative_sum = 1.0
        for occupation in range(2, positions):
            relative_term *= (
                (positions - occupation) * x / (occupation + 1)
            )
            relative_sum += relative_term
            if relative_term <= 1e-17 * relative_sum:
                break
        return log_term_two + math.log(relative_sum)

    log_one_plus_x = float(np.logaddexp(0.0, log_x))
    total_log = positions * log_one_plus_x
    removed_log = float(np.logaddexp(0.0, math.log(positions) + log_x))
    gap = total_log - removed_log
    if not gap > 0.0:
        raise ArithmeticError("binomial-tail subtraction lost positivity")
    return total_log + math.log(-math.expm1(-gap))


def occupation_tail_checks() -> list[dict[str, object]]:
    rows = []
    for positions in (3, 5, 8):
        for x in (1e-5, 0.03, 0.5, 2.0):
            exact = sum(
                math.comb(positions, occupation) * x**occupation
                for occupation in range(2, positions + 1)
            )
            observed = math.exp(
                log_binomial_tail_from_two(math.log(x), positions)
            )
            if not math.isclose(observed, exact, rel_tol=3e-13, abs_tol=1e-300):
                raise AssertionError("stable occupation-tail evaluator failed")
            rows.append(
                {
                    "positions": positions,
                    "x": x,
                    "relative_error": abs(observed - exact) / exact,
                    "status": "MATCH",
                }
            )
    return rows


def all_occupation_log_moment(
    model: LocalMomentModel, log_u: float
) -> tuple[float, dict[str, float]]:
    log_m1, log_m2_nonzero, log_m2_all, log_m3 = model.log_moments(log_u)
    log_one_data = math.log(DATA_POSITIONS) + log_m3
    log_many_data = (
        log_m2_nonzero
        + log_m2_all
        - 2.0 * log_m1
        + log_binomial_tail_from_two(log_m1)
    )
    total = float(np.logaddexp(log_one_data, log_many_data))
    if log_m1 > 40.0:
        mode = float(DATA_POSITIONS)
    else:
        x = math.exp(log_m1)
        mode = max(2.0, math.floor((DATA_POSITIONS + 1) * x / (1.0 + x)))
    return total, {
        "log2_one_data_component": log_one_data / LOG2,
        "log2_many_data_component": log_many_data / LOG2,
        "many_data_occupation_mode": mode,
        "log2_M1": log_m1 / LOG2,
    }


def inner_parameters(row: dict[str, object]) -> dict[str, float]:
    z = float(row["left_outward"]["verified_weight_tilt_decimal"])
    radius = float(row["left_outward"]["verified_coefficient_radius_decimal"])
    distance = int(row["distance"])
    q = math.ldexp(1.0, -TARGET_SIGMA)
    b = math.ldexp((1.0 + z) ** TARGET_G, -TARGET_G)
    d = (1.0 - q) * b
    gap = (b - d * radius) / ((1.0 - radius) * (1.0 - d * radius))
    c = radius * gap
    constant = (
        -distance * math.log(z)
        - TARGET_PACKETS * math.log(radius)
        - math.log1p(-radius)
    )
    for support in (int(row["support_left"]), int(row["support_right"])):
        direct = constant + support * math.log(c) - log_binomial(TARGET_PACKETS, support)
        reference = log_common_bound(
            TARGET_PACKETS,
            TARGET_G,
            TARGET_SIGMA,
            support,
            distance,
            z,
            radius,
        )
        if not math.isclose(direct, reference, rel_tol=0.0, abs_tol=2e-7):
            raise AssertionError("analytic inner form failed the Goal 02 gate")
    return {
        "weight_tilt": z,
        "coefficient_radius": radius,
        "log_constant": constant,
        "log_support_base": math.log(c),
    }


def refined_segments(rows: list[dict[str, object]]) -> list[tuple[int, int, dict[str, object]]]:
    boundaries = {18, TARGET_PACKETS + 1, TARGET_PACKETS // 2 + 1}
    for row in rows:
        boundaries.add(int(row["support_left"]))
        boundaries.add(int(row["support_right"]) + 1)
    for boundary in tuple(boundaries):
        reflected = TARGET_PACKETS - boundary + 1
        if 18 <= reflected <= TARGET_PACKETS:
            boundaries.add(reflected)
    ordered = sorted(boundary for boundary in boundaries if 18 <= boundary <= TARGET_PACKETS + 1)
    segments = []
    row_index = 0
    for left, stop in zip(ordered, ordered[1:]):
        right = stop - 1
        while left > int(rows[row_index]["support_right"]):
            row_index += 1
        if right > int(rows[row_index]["support_right"]):
            raise AssertionError("refined segment crosses an inner-certificate boundary")
        segments.append((left, right, rows[row_index]))
    if segments[0][0] != 18 or segments[-1][1] != TARGET_PACKETS:
        raise AssertionError("refined segments do not cover the support range")
    return segments


def optimize_weighted_interval(
    model: LocalMomentModel,
    left: int,
    right: int,
    log_support_base: float,
) -> dict[str, object]:
    unconditioned, unconditioned_meta = all_occupation_log_moment(
        model, log_support_base
    )

    def optimize_side(side: str) -> tuple[float, float, bool, str]:
        def objective(log_surprisal: float) -> float:
            magnitude = math.exp(log_surprisal)
            if side == "lower":
                log_s = -magnitude
                cutoff = right
            else:
                log_s = magnitude
                cutoff = left
            moment, _meta = all_occupation_log_moment(
                model, log_support_base + log_s
            )
            return moment - cutoff * log_s

        grid = np.linspace(-16.0, 3.0, 77)
        values = np.asarray([objective(float(point)) for point in grid])
        index = int(np.argmin(values))
        low = float(grid[max(0, index - 1)])
        high = float(grid[min(len(grid) - 1, index + 1)])
        result = minimize_scalar(
            objective,
            bounds=(low, high),
            method="bounded",
            options={"xatol": 1e-11, "maxiter": 500},
        )
        return float(result.fun), math.exp(float(result.x)), bool(result.success), str(result.message)

    lower, lower_magnitude, lower_success, lower_message = optimize_side("lower")
    upper, upper_magnitude, upper_success, upper_message = optimize_side("upper")
    candidates = (
        ("unconditioned", unconditioned, 0.0, True, "no support tilt"),
        ("lower_tail", lower, -lower_magnitude, lower_success, lower_message),
        ("upper_tail", upper, upper_magnitude, upper_success, upper_message),
    )
    selected = min(candidates, key=lambda item: item[1])
    selected_moment, selected_meta = all_occupation_log_moment(
        model, log_support_base + selected[2]
    )
    return {
        "selected_tail": selected[0],
        "weighted_log_moment_upper": selected[1],
        "support_log_tilt": selected[2],
        "support_tilt": math.exp(selected[2]),
        "optimizer_success": selected[3],
        "optimizer_message": selected[4],
        "moment_at_selected_tilt_log2": selected_moment / LOG2,
        "one_data_component_log2": selected_meta["log2_one_data_component"],
        "many_data_component_log2": selected_meta["log2_many_data_component"],
        "many_data_occupation_mode": selected_meta["many_data_occupation_mode"],
        "unconditioned_many_data_occupation_mode": unconditioned_meta[
            "many_data_occupation_mode"
        ],
    }


def evaluate_distance(
    model: LocalMomentModel, rows: list[dict[str, object]]
) -> dict[str, object]:
    segments = refined_segments(rows)
    results = []
    for left, right, source_row in segments:
        parameters = inner_parameters(source_row)
        negative_log_binomial_left = -log_binomial(TARGET_PACKETS, left)
        negative_log_binomial_right = -log_binomial(TARGET_PACKETS, right)
        if left == right:
            secant_slope = 0.0
            secant_intercept = negative_log_binomial_left
        else:
            secant_slope = (
                negative_log_binomial_right - negative_log_binomial_left
            ) / (right - left)
            secant_intercept = negative_log_binomial_left - secant_slope * left
            midpoint = (left + right) // 2
            if (
                secant_intercept + secant_slope * midpoint + 1e-7
                < -log_binomial(TARGET_PACKETS, midpoint)
            ):
                raise AssertionError("reciprocal-binomial secant failed convexity gate")
        weighted = optimize_weighted_interval(
            model,
            left,
            right,
            parameters["log_support_base"] + secant_slope,
        )
        contribution = (
            parameters["log_constant"]
            + secant_intercept
            + float(weighted["weighted_log_moment_upper"])
        )
        result: dict[str, object] = {
            "support_left": left,
            "support_right": right,
            "source_certificate_left": int(source_row["support_left"]),
            "source_certificate_right": int(source_row["support_right"]),
            "reciprocal_binomial_secant_intercept_log2": secant_intercept / LOG2,
            "reciprocal_binomial_secant_slope_log2": secant_slope / LOG2,
            "log2_contribution_upper": contribution / LOG2,
        }
        result.update(parameters)
        result.update(weighted)
        results.append(result)
        print(
            f"delta,{float(source_row['relative_binary_weight']):.6f},"
            f"segment,{left}:{right},log2,{contribution / LOG2:.6f},"
            f"tail,{weighted['selected_tail']},"
            f"mode,{weighted['many_data_occupation_mode']:.0f}",
            flush=True,
        )
    total = float(logsumexp([float(row["log2_contribution_upper"]) * LOG2 for row in results])) / LOG2
    dominant = max(results, key=lambda row: float(row["log2_contribution_upper"]))
    return {
        "relative_binary_weight": float(rows[0]["relative_binary_weight"]),
        "distance": int(rows[0]["distance"]),
        "refined_segment_count": len(results),
        "total_log2_first_moment_upper": total,
        "dominant_support_interval": [
            int(dominant["support_left"]),
            int(dominant["support_right"]),
        ],
        "dominant_occupation_mode": dominant["many_data_occupation_mode"],
        "segments": results,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--inner-receipt", type=Path, default=DEFAULT_INNER_RECEIPT)
    parser.add_argument("--relative-distances", nargs="+", type=float, default=[0.05, 0.09, 0.12])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spectrum = load_spectrum(args.spectrum)
    model = LocalMomentModel(spectrum)
    tail_checks = occupation_tail_checks()
    grouped = load_inner_caps(args.inner_receipt)
    requested = tuple(args.relative_distances)
    results = []
    for delta in requested:
        if delta not in grouped:
            raise ValueError(f"inner receipt does not contain relative distance {delta}")
        result = evaluate_distance(model, grouped[delta])
        results.append(result)
        print(
            f"result,delta,{delta:.6f},total_log2,"
            f"{result['total_log2_first_moment_upper']:.6f},"
            f"dominant,{result['dominant_support_interval']},"
            f"mode,{result['dominant_occupation_mode']}",
            flush=True,
        )
    payload = {
        "schema": "riffle-randomstepconv-goal04-v1",
        "construction": "Riffle RandomStepConv g=4 sigma=20",
        "evidence": {
            "occupation_sum": "EXACT_ALGEBRAIC_IDENTITY",
            "holder_and_restricted_cauchy": "EXACT",
            "bch_spectrum_and_packet_support_laws": "EXACT_VALIDATED",
            "inner_certificate_interface": "EXACT_FORM_WITH_VALID_GOAL02_TILTS",
            "optimization_and_evaluation": "FLOATING_DIAGNOSTIC",
        },
        "parameters": {
            "data_positions": DATA_POSITIONS,
            "packet_positions": TARGET_PACKETS,
            "relative_distances": list(requested),
        },
        "occupation_tail_checks": tail_checks,
        "results": results,
        "scope": (
            "The algebra covers all nonzero outer words and all occupations. "
            "The reported numerical values are floating-point diagnostics, not "
            "outward-rounded end-to-end certificates."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
