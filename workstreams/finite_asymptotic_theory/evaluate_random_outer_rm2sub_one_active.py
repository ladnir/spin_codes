#!/usr/bin/env python3
"""Binary64 one-active diagnostic for random outer blocks and RM2Sub-S19.

The script imports the audited constituent data and the existing two-state
entrywise transfer implementation.  It changes only the local outer spectrum:
each block is a uniform random linear injection.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (  # noqa: E402
    LOG2,
    load_nonzero_spectrum,
    load_outer_spectrum_logs,
    log_choose,
    logsumexp,
    splitstate_impulse_matrices,
)
from analyze_riffle_splitstate_preaddmul_one_active import (  # noqa: E402
    marked_region,
    weight_conditioned_log_moments,
)


DEFAULT_ACTIVATION = REPO_ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_b_kernel_spectrum.json"
)
DEFAULT_LIVE_SPECTRUM = REPO_ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_a_spectrum_audit.json"
)
DEFAULT_COMPARISON_SPECTRUM = REPO_ROOT / (
    "constructions/riffle_parityshear12_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/receipts/"
    "parityfanout31x33_outer256_d38_expected_spectrum.json"
)


def log_two_power_minus_one(bits: int) -> float:
    return bits * LOG2 + math.log1p(-math.ldexp(1.0, -bits))


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(resolved)


def random_injection_spectrum_logs(length: int, dimension: int) -> np.ndarray:
    if not 0 < dimension <= length:
        raise ValueError("outer dimension must lie in [1, outer length]")
    result = np.full(length + 1, -math.inf, dtype=np.float64)
    result[0] = 0.0
    density = log_two_power_minus_one(dimension) - log_two_power_minus_one(length)
    for weight in range(1, length + 1):
        result[weight] = density + log_choose(length, weight)
    return result


def load_nonactivation(path: Path, step_bits: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = np.full(step_bits + 1, np.nan, dtype=np.float64)
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > step_bits:
            continue
        key = "uniform_support_average_distinct_upper_bound"
        if key not in row:
            key = "uniform_256_support_average_distinct_upper_bound"
        result[weight] = float(row[key])
    missing = np.flatnonzero(np.isnan(result))
    if missing.size:
        raise ValueError(f"activation table is missing weight {int(missing[0])}")
    return result


def tilt_grid(minimum: float, maximum: float, spacing: float) -> np.ndarray:
    if not minimum < maximum:
        raise ValueError("tilt minimum must be less than tilt maximum")
    if spacing <= 0.0:
        raise ValueError("tilt spacing must be positive")
    count = math.floor((maximum - minimum) / spacing + 1e-12)
    return minimum + spacing * np.arange(count + 1, dtype=np.float64)


def log_hamming_ball(length: int, radius: int) -> float:
    if not 0 <= radius <= length:
        raise ValueError("Hamming-ball radius must lie in [0, length]")
    if radius > length // 2:
        raise ValueError("diagnostic recurrence expects radius at most length/2")
    term = 1.0
    normalized_sum = 1.0
    for weight in range(radius, 0, -1):
        term *= weight / (length - weight + 1)
        normalized_sum += term
        if term == 0.0:
            break
    return log_choose(length, radius) + math.log(normalized_sum)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.outer_blocks <= 0:
        raise ValueError("outer-block count must be positive")
    if args.outer_blocks % args.step_bits:
        raise ValueError("step size must divide the number of outer blocks")

    outer_bits = args.outer_bits
    outer_dimension = args.outer_dimension
    output_bits = args.outer_blocks * outer_bits
    target_distance = math.floor(args.relative_distance * output_bits)
    epochs_per_region = args.outer_blocks // args.step_bits

    spectrum = random_injection_spectrum_logs(outer_bits, outer_dimension)
    nonactivation = load_nonactivation(args.activation, args.step_bits)
    live_spectrum = load_nonzero_spectrum(
        args.live_spectrum, args.step_bits, args.state_bits
    )

    best_by_weight = np.full(outer_bits + 1, math.inf, dtype=np.float64)
    best_tilt_by_weight = np.full(outer_bits + 1, math.nan, dtype=np.float64)
    best_common = math.inf
    best_common_tilt = math.nan
    density_log = (
        log_two_power_minus_one(outer_dimension)
        - log_two_power_minus_one(outer_bits)
    )

    for log_surprisal in tilt_grid(
        args.tilt_minimum, args.tilt_maximum, args.tilt_spacing
    ):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        impulse = splitstate_impulse_matrices(
            z=z,
            step_bits=args.step_bits,
            state_bits=args.state_bits,
            constituent_distance=args.constituent_distance,
            live_moment_order=3,
            live_model="support-averaged-preaddmul",
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        inactive_region, active_region = marked_region(
            impulse[0], impulse[1], epochs_per_region
        )
        moments = weight_conditioned_log_moments(
            outer_bits=outer_bits,
            zero_row=inactive_region,
            active_row=active_region,
        )

        tilted = moments + target_distance * surprisal
        improved = tilted < best_by_weight
        best_by_weight[improved] = tilted[improved]
        best_tilt_by_weight[improved] = log_surprisal

        nonzero_terms = np.asarray(
            [
                log_choose(outer_bits, weight) + float(moments[weight])
                for weight in range(1, outer_bits + 1)
            ],
            dtype=np.float64,
        )
        common = (
            math.log(args.outer_blocks)
            + density_log
            + float(logsumexp(nonzero_terms))
            + target_distance * surprisal
        )
        if common < best_common:
            best_common = common
            best_common_tilt = float(log_surprisal)

    rows: list[dict[str, float | int]] = []
    aggregate_terms = []
    for weight in range(1, outer_bits + 1):
        inner_log = min(0.0, float(best_by_weight[weight]))
        pointwise = (
            math.log(args.outer_blocks)
            + float(spectrum[weight])
            + inner_log
        )
        aggregate_terms.append(pointwise)
        rows.append(
            {
                "outer_weight": weight,
                "best_log_surprisal": float(best_tilt_by_weight[weight]),
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": pointwise / LOG2,
                "pointwise_margin_bits": -pointwise / LOG2,
            }
        )

    aggregate = float(logsumexp(np.asarray(aggregate_terms)))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    log_beta = (
        log_two_power_minus_one(outer_dimension)
        - math.log1p(-math.ldexp(1.0, -outer_bits))
    )
    all_active_log_upper = (
        args.outer_blocks * log_beta
        - output_bits * LOG2
        + log_hamming_ball(output_bits, target_distance)
    )
    payload: dict[str, object] = {
        "schema": "random-outer-rm2sub-one-active-diagnostic-v1",
        "construction": "random outer x factored region shuffle x RM2Sub-S19",
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": outer_dimension,
            "outer_blocks": args.outer_blocks,
            "output_bits": output_bits,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "epochs_per_region": epochs_per_region,
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
            "activation": display_path(args.activation),
            "live_spectrum": display_path(args.live_spectrum),
            "tilt_minimum": args.tilt_minimum,
            "tilt_maximum": args.tilt_maximum,
            "tilt_spacing": args.tilt_spacing,
        },
        "per_weight_aggregate_log2_upper": aggregate / LOG2,
        "per_weight_aggregate_margin_bits": -aggregate / LOG2,
        "common_tilt_log2_upper": best_common / LOG2,
        "common_tilt_margin_bits": -best_common / LOG2,
        "best_common_log_surprisal": best_common_tilt,
        "dominant_row": dominant,
        "all_active": {
            "log2_upper": all_active_log_upper / LOG2,
            "margin_bits": -all_active_log_upper / LOG2,
            "method": "nonzero-row relaxation plus inner bijectivity",
        },
        "weight_rows": rows,
        "scope": (
            "Binary64 diagnostic. The random-injection outer spectrum is exact "
            "in form. The RM2Sub constituent receipts are exact inputs, but the "
            "two-state transfer remains an entrywise envelope and the numerical "
            "optimization is not outward rounded."
        ),
    }

    if args.comparison_spectrum is not None:
        comparison_spectrum = load_outer_spectrum_logs(
            args.comparison_spectrum, outer_bits
        )
        comparison_rows = []
        comparison_terms = []
        for weight in range(1, outer_bits + 1):
            if not math.isfinite(float(comparison_spectrum[weight])):
                continue
            inner_log = min(0.0, float(best_by_weight[weight]))
            pointwise = (
                math.log(args.outer_blocks)
                + float(comparison_spectrum[weight])
                + inner_log
            )
            comparison_terms.append(pointwise)
            comparison_rows.append(
                {
                    "outer_weight": weight,
                    "pointwise_log2_upper": pointwise / LOG2,
                    "pointwise_margin_bits": -pointwise / LOG2,
                }
            )
        comparison_aggregate = float(
            logsumexp(np.asarray(comparison_terms, dtype=np.float64))
        )
        comparison_dominant = max(
            comparison_rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
        )
        payload["comparison"] = {
            "spectrum": display_path(args.comparison_spectrum),
            "aggregate_log2_upper": comparison_aggregate / LOG2,
            "aggregate_margin_bits": -comparison_aggregate / LOG2,
            "margin_difference_vs_random_bits": (
                (-comparison_aggregate + aggregate) / LOG2
            ),
            "dominant_row": comparison_dominant,
        }

    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--outer-dimension", type=int, default=128)
    parser.add_argument("--outer-blocks", type=int, default=8192)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=19)
    parser.add_argument("--constituent-distance", type=int, default=48)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--live-spectrum", type=Path, default=DEFAULT_LIVE_SPECTRUM)
    parser.add_argument(
        "--comparison-spectrum",
        type=Path,
        default=DEFAULT_COMPARISON_SPECTRUM,
        help="optional structured spectrum evaluated with the same inner tilts",
    )
    parser.add_argument("--tilt-minimum", type=float, default=-10.0)
    parser.add_argument("--tilt-maximum", type=float, default=-6.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = evaluate(args)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote,{args.output}")
    print(
        "per_weight_margin_bits,"
        f"{payload['per_weight_aggregate_margin_bits']:.9f}"
    )
    print(f"common_tilt_margin_bits,{payload['common_tilt_margin_bits']:.9f}")
    print(f"all_active_margin_bits,{payload['all_active']['margin_bits']:.9f}")
    print(
        "dominant_outer_weight,"
        f"{payload['dominant_row']['outer_weight']},"
        "pointwise_margin_bits,"
        f"{payload['dominant_row']['pointwise_margin_bits']:.9f}"
    )
    if "comparison" in payload:
        comparison = payload["comparison"]
        print(
            "comparison_margin_bits,"
            f"{comparison['aggregate_margin_bits']:.9f},"
            "comparison_minus_random_bits,"
            f"{comparison['margin_difference_vs_random_bits']:.9f}"
        )


if __name__ == "__main__":
    main()
