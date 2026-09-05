#!/usr/bin/env python3
"""Evaluate spectrum-only outer moment gates for RandomStepConv Goal 03."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from analyze_riffle_bchblockperm_parallelacc_g4_goal01 import (
    BLOCK_BITS,
    PACKETS_PER_BLOCK,
    packet_support_count,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal02 import (
    TARGET_G,
    TARGET_PACKETS,
    TARGET_SIGMA,
    log_common_bound,
)


FIELD_BITS = 64
FIELD_SIZE = 1 << FIELD_BITS
DATA_POSITIONS = 1 << 14
DEFAULT_OCCUPATIONS = (1, 2, 3, 4, 8, 16)
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
DEFAULT_SPECTRUM = SCRIPT_DIRECTORY / "ebch128_64_spectrum.csv"
DEFAULT_INNER_RECEIPT = Path(
    "constructions/riffle_randomstepconv_g4_sigma20/receipts/"
    "goal02_all_support_envelope.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_randomstepconv_g4_sigma20/receipts/"
    "goal03_bch_moment_gate.json"
)
LOG2 = math.log(2.0)


class LocalMomentModel:
    def __init__(self, spectrum: dict[int, int]) -> None:
        self.weights = np.asarray(sorted(spectrum), dtype=np.int64)
        self.log_multiplicities = np.asarray(
            [math.log(spectrum[int(weight)]) for weight in self.weights],
            dtype=np.float64,
        )
        self.log_support_probabilities = np.full(
            (len(self.weights), PACKETS_PER_BLOCK + 1),
            -math.inf,
            dtype=np.float64,
        )
        for row, weight_value in enumerate(self.weights):
            weight = int(weight_value)
            if weight == 0:
                self.log_support_probabilities[row, 0] = 0.0
                continue
            denominator_log = math.log(math.comb(BLOCK_BITS, weight))
            total = 0
            for support in range(1, PACKETS_PER_BLOCK + 1):
                count = packet_support_count(weight, support)
                if not count:
                    continue
                total += count
                self.log_support_probabilities[row, support] = (
                    math.log(count) - denominator_log
                )
            if total != math.comb(BLOCK_BITS, weight):
                raise AssertionError(f"packet-support law failed at BCH weight {weight}")
        self.nonzero = self.weights != 0
        self.supports = np.arange(PACKETS_PER_BLOCK + 1, dtype=np.float64)

    def log_q(self, log_t: float) -> np.ndarray:
        return logsumexp(
            self.log_support_probabilities + self.supports * log_t,
            axis=1,
        )

    def log_moments(self, log_t: float) -> tuple[float, float, float, float]:
        log_q = self.log_q(log_t)
        log_m1 = float(
            logsumexp(
                self.log_multiplicities[self.nonzero] + log_q[self.nonzero]
            )
        )
        log_m2_nonzero = float(
            logsumexp(
                self.log_multiplicities[self.nonzero]
                + 2.0 * log_q[self.nonzero]
            )
        )
        log_m2_all = float(logsumexp(self.log_multiplicities + 2.0 * log_q))
        log_m3 = float(
            logsumexp(
                self.log_multiplicities[self.nonzero]
                + 3.0 * log_q[self.nonzero]
            )
        )
        return log_m1, log_m2_nonzero, log_m2_all, log_m3


def load_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        spectrum = {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(source)
        }
    if sum(spectrum.values()) != FIELD_SIZE:
        raise AssertionError("BCH spectrum mass is not 2^64")
    if spectrum.get(0) != 1 or spectrum.get(BLOCK_BITS) != 1:
        raise AssertionError("BCH spectrum endpoint check failed")
    if min(weight for weight, count in spectrum.items() if weight and count) != 22:
        raise AssertionError("BCH minimum-weight check failed")
    if any(spectrum.get(weight, 0) != spectrum.get(BLOCK_BITS - weight, 0) for weight in spectrum):
        raise AssertionError("BCH complement symmetry check failed")
    return spectrum


def log_shell_moment(
    model: LocalMomentModel,
    occupation: int,
    log_t: float,
    method: str,
) -> float:
    log_m1, log_m2_nonzero, log_m2_all, log_m3 = model.log_moments(log_t)
    if method == "cauchy":
        if occupation == 1:
            return log_m3
        return (
            (occupation - 2) * log_m1
            + log_m2_nonzero
            + log_m2_all
        )
    if method == "parity_worst":
        # For t<=1, each nonzero BCH block occupies at least six packets.
        guaranteed_parity_packets = 12 if occupation == 1 else 6 if occupation == 2 else 0
        return occupation * log_m1 + guaranteed_parity_packets * log_t
    raise ValueError(f"unknown method {method}")


def optimize_cumulative_count(
    model: LocalMomentModel,
    occupation: int,
    cutoff: int,
    method: str,
) -> dict[str, object]:
    exact_total_log = occupation * math.log(FIELD_SIZE - 1)

    def objective(log_surprisal: float) -> float:
        log_t = -math.exp(log_surprisal)
        return (
            log_shell_moment(model, occupation, log_t, method)
            - cutoff * log_t
        )

    grid = np.linspace(-16.0, 3.0, 77)
    values = np.asarray([objective(float(point)) for point in grid])
    index = int(np.argmin(values))
    low = float(grid[max(0, index - 1)])
    high = float(grid[min(len(grid) - 1, index + 1)])
    optimized = minimize_scalar(
        objective,
        bounds=(low, high),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 500},
    )
    selected_log = min(exact_total_log, float(optimized.fun))
    return {
        "log_count_upper": selected_log,
        "log2_count_upper": selected_log / LOG2,
        "support_tilt": math.exp(-math.exp(float(optimized.x))),
        "optimizer_success": bool(optimized.success),
        "optimizer_message": str(optimized.message),
        "used_exact_total": exact_total_log <= float(optimized.fun),
    }


def load_inner_caps(path: Path) -> dict[float, list[dict[str, object]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    grouped: dict[float, list[dict[str, object]]] = {}
    for row in payload["rows"]:
        grouped.setdefault(float(row["relative_binary_weight"]), []).append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: int(row["support_left"]))
        caps = [float(row["outward_interval_cap_log2"]) * LOG2 for row in rows]
        if any(caps[index] < caps[index + 1] for index in range(len(caps) - 1)):
            raise AssertionError("inner interval caps are not decreasing")
    return grouped


def logdiffexp(left: float, right: float) -> float:
    if right == -math.inf:
        return left
    if not left > right:
        raise ValueError("logdiffexp requires left > right")
    return left + math.log1p(-math.exp(right - left))


def shell_first_moment(
    model: LocalMomentModel,
    occupation: int,
    inner_rows: list[dict[str, object]],
    method: str,
) -> dict[str, object]:
    maximum_support = min(
        32 * (occupation + 2), int(inner_rows[-1]["support_right"])
    )
    used_rows = [
        row for row in inner_rows if int(row["support_left"]) <= maximum_support
    ]
    caps = [float(row["outward_interval_cap_log2"]) * LOG2 for row in used_rows]
    cumulative: list[dict[str, object]] = []
    exact_total_log = occupation * math.log(FIELD_SIZE - 1)
    for index, row in enumerate(used_rows):
        cutoff = min(maximum_support, int(row["support_right"]))
        if index == len(used_rows) - 1:
            result = {
                "log_count_upper": exact_total_log,
                "log2_count_upper": exact_total_log / LOG2,
                "support_tilt": 1.0,
                "optimizer_success": True,
                "optimizer_message": "exact total at maximum support",
                "used_exact_total": True,
            }
        else:
            result = optimize_cumulative_count(
                model, occupation, cutoff, method
            )
        result["cutoff"] = cutoff
        cumulative.append(result)

    terms = []
    for index in range(len(used_rows) - 1):
        coefficient = logdiffexp(caps[index], caps[index + 1])
        terms.append(coefficient + float(cumulative[index]["log_count_upper"]))
    terms.append(caps[-1] + exact_total_log)
    per_position_set_log = float(logsumexp(terms))
    position_log = (
        math.lgamma(DATA_POSITIONS + 1)
        - math.lgamma(occupation + 1)
        - math.lgamma(DATA_POSITIONS - occupation + 1)
    )
    return {
        "method": method,
        "occupation": occupation,
        "maximum_packet_support": maximum_support,
        "intervals_used": len(used_rows),
        "position_choices_log2": position_log / LOG2,
        "per_position_set_log2_first_moment": per_position_set_log / LOG2,
        "shell_log2_first_moment": (position_log + per_position_set_log) / LOG2,
        "cumulative_bounds": cumulative,
    }


def pointwise_inner_logs(
    inner_rows: list[dict[str, object]], maximum_support: int
) -> tuple[list[int], list[float]]:
    supports = list(range(18, maximum_support + 1))
    logs = []
    row_index = 0
    for support in supports:
        while support > int(inner_rows[row_index]["support_right"]):
            row_index += 1
        row = inner_rows[row_index]
        z = float(row["left_outward"]["verified_weight_tilt_decimal"])
        radius = float(row["left_outward"]["verified_coefficient_radius_decimal"])
        value = log_common_bound(
            TARGET_PACKETS,
            TARGET_G,
            TARGET_SIGMA,
            support,
            int(row["distance"]),
            z,
            radius,
        )
        logs.append(min(0.0, value))
    # A decreasing majorant permits summation by parts even if independent
    # interval tilts create a small upward jump at a boundary.
    for index in range(len(logs) - 2, -1, -1):
        logs[index] = max(logs[index], logs[index + 1])
    return supports, logs


def pointwise_shell_first_moment(
    model: LocalMomentModel,
    occupation: int,
    inner_rows: list[dict[str, object]],
    method: str,
    cumulative_cache: dict[tuple[int, str, int], dict[str, object]],
) -> dict[str, object]:
    maximum_support = min(
        32 * (occupation + 2), int(inner_rows[-1]["support_right"])
    )
    supports, inner_logs = pointwise_inner_logs(inner_rows, maximum_support)
    exact_total_log = occupation * math.log(FIELD_SIZE - 1)
    cumulative = []
    for support in supports:
        key = (occupation, method, support)
        if support == maximum_support:
            result = {
                "log_count_upper": exact_total_log,
                "log2_count_upper": exact_total_log / LOG2,
                "support_tilt": 1.0,
                "optimizer_success": True,
                "optimizer_message": "exact total at maximum support",
                "used_exact_total": True,
            }
        elif key in cumulative_cache:
            result = cumulative_cache[key]
        else:
            result = optimize_cumulative_count(
                model, occupation, support, method
            )
            cumulative_cache[key] = result
        cumulative.append(result)

    terms: list[tuple[int, float]] = []
    for index in range(len(supports) - 1):
        if inner_logs[index] <= inner_logs[index + 1] + 1e-14:
            continue
        coefficient = logdiffexp(inner_logs[index], inner_logs[index + 1])
        terms.append(
            (
                supports[index],
                coefficient + float(cumulative[index]["log_count_upper"]),
            )
        )
    terms.append((supports[-1], inner_logs[-1] + exact_total_log))
    position_log = (
        math.lgamma(DATA_POSITIONS + 1)
        - math.lgamma(occupation + 1)
        - math.lgamma(DATA_POSITIONS - occupation + 1)
    )
    dominant_support, dominant_log = max(terms, key=lambda item: item[1])
    per_position_set_log = float(logsumexp([value for _support, value in terms]))
    return {
        "method": method,
        "occupation": occupation,
        "support_count": len(supports),
        "minimum_packet_support": supports[0],
        "maximum_packet_support": supports[-1],
        "position_choices_log2": position_log / LOG2,
        "per_position_set_log2_first_moment": per_position_set_log / LOG2,
        "shell_log2_first_moment": (position_log + per_position_set_log) / LOG2,
        "dominant_cumulative_cutoff": dominant_support,
        "dominant_term_log2_with_positions": (dominant_log + position_log) / LOG2,
        "largest_terms": [
            {
                "cumulative_cutoff": support,
                "log2_with_positions": (value + position_log) / LOG2,
            }
            for support, value in sorted(terms, key=lambda item: item[1], reverse=True)[:8]
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--inner-receipt", type=Path, default=DEFAULT_INNER_RECEIPT)
    parser.add_argument(
        "--occupations", nargs="+", type=int, default=list(DEFAULT_OCCUPATIONS)
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spectrum = load_spectrum(args.spectrum)
    model = LocalMomentModel(spectrum)
    inner = load_inner_caps(args.inner_receipt)
    rows = []
    cumulative_cache: dict[tuple[int, str, int], dict[str, object]] = {}
    for delta, inner_rows in sorted(inner.items()):
        for occupation in args.occupations:
            if not 1 <= occupation <= DATA_POSITIONS:
                raise ValueError("occupation must lie in [1,16384]")
            methods = [
                shell_first_moment(model, occupation, inner_rows, method)
                for method in ("cauchy", "parity_worst")
            ]
            for method_row in methods:
                method_row["pointwise_refinement"] = pointwise_shell_first_moment(
                    model,
                    occupation,
                    inner_rows,
                    str(method_row["method"]),
                    cumulative_cache,
                )
            best = min(
                methods,
                key=lambda row: float(
                    row["pointwise_refinement"]["shell_log2_first_moment"]
                ),
            )
            row = {
                "relative_binary_weight": delta,
                "occupation": occupation,
                "best_method": best["method"],
                "best_shell_log2_first_moment": best["pointwise_refinement"][
                    "shell_log2_first_moment"
                ],
                "methods": methods,
            }
            rows.append(row)
            print(
                f"delta,{delta:.6f},r,{occupation},best,{row['best_method']},"
                f"shell_log2,{row['best_shell_log2_first_moment']:.6f},"
                f"cauchy,{methods[0]['pointwise_refinement']['shell_log2_first_moment']:.6f},"
                f"parity_worst,{methods[1]['pointwise_refinement']['shell_log2_first_moment']:.6f},"
                f"dominant_h,{best['pointwise_refinement']['dominant_cumulative_cutoff']}",
                flush=True,
            )

    payload = {
        "schema": "riffle-randomstepconv-goal03-v1",
        "construction": "Riffle RandomStepConv g=4 sigma=20",
        "evidence": {
            "holder_and_cauchy_inequalities": "EXACT",
            "bch_spectrum": "EXACT_VALIDATED",
            "packet_support_laws": "EXACT_INTEGER_VALIDATED",
            "inner_caps": "OUTWARD_ROUNDED_FROM_GOAL02",
            "pointwise_inner_evaluation_at_valid_tilts": "FLOATING_DIAGNOSTIC",
            "moment_optimization": "FLOATING_DIAGNOSTIC",
        },
        "parameters": {
            "field_bits": FIELD_BITS,
            "data_positions": DATA_POSITIONS,
            "bch_block_bits": BLOCK_BITS,
            "packets_per_bch_block": PACKETS_PER_BLOCK,
            "occupations": list(args.occupations),
        },
        "spectrum_validation": {
            "total": sum(spectrum.values()),
            "minimum_nonzero_weight": min(weight for weight in spectrum if weight),
            "weight_classes": len(spectrum),
            "complement_symmetric": True,
        },
        "coefficient_matrix_gate": {
            "determinant_for_positions_i_j": "alpha_j-alpha_i",
            "nonzero_for_distinct_positions": True,
            "status": "EXACT_FIELD_ARGUMENT",
        },
        "rows": rows,
        "scope": (
            "This is a floating-point low-shell diagnostic built from exact "
            "outer inequalities and the certified Goal 02 inner caps. It is not "
            "an outward-rounded end-to-end distance certificate."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
