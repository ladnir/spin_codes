#!/usr/bin/env python3
"""Replace the Riffle bulk worst-density surcharge by full-spectrum moments.

The regular-bulk evaluator bounds the bad-event probability under independent
uniform ambient words for each active outer block.  If Q is the distribution
induced by the supplied outer spectrum after its output-coordinate permutation,
write L=dQ/dU.  For r active blocks and p>1, Holder gives

    E_U[1_bad product_i L_i]
      <= E_U[L^p]^(r/p) P_U[bad]^((p-1)/p).

Optimizing p uses the complete spectrum.  The old maximum density ratio is the
p=infinity endpoint.  This remains a floating-point diagnostic because the
input bulk probability bound and spectrum receipt are floating point.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    load_outer_spectrum_logs,
    log_choose,
    logsumexp,
    log_two_power_minus_one,
)


def holder_grid() -> np.ndarray:
    # Dense in log2(p-1), including the neighborhood of p=1 and a practical
    # approximation to the infinity endpoint.
    return 1.0 + np.exp2(np.linspace(-16.0, 16.0, 513))


def evaluate(
    bulk_path: Path,
    spectrum_path: Path,
    one_active_path: Path | None,
) -> dict[str, object]:
    bulk = json.loads(bulk_path.read_text(encoding="utf-8"))
    parameters = bulk["parameters"]
    outer_bits = int(parameters["outer_bits"])
    dimension = outer_bits // 2
    outer_blocks = int(parameters["outer_blocks"])
    spectrum = load_outer_spectrum_logs(spectrum_path, outer_bits)

    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    weights: list[int] = []
    log_reference_probabilities: list[float] = []
    log_likelihoods: list[float] = []
    for weight in range(1, outer_bits + 1):
        log_multiplicity = float(spectrum[weight])
        if not math.isfinite(log_multiplicity):
            continue
        log_shell = log_choose(outer_bits, weight)
        weights.append(weight)
        log_reference_probabilities.append(log_shell - outer_bits * LOG2)
        log_likelihoods.append(
            log_multiplicity - log_shell - reference_log_density
        )

    log_reference = np.asarray(log_reference_probabilities, dtype=np.float64)
    log_likelihood = np.asarray(log_likelihoods, dtype=np.float64)
    normalization = float(logsumexp(log_reference + log_likelihood))
    if abs(normalization) > 1e-8:
        raise AssertionError(
            f"outer likelihood does not normalize: log expectation={normalization}"
        )

    p_values = holder_grid()
    log_moments = np.asarray(
        [
            logsumexp(log_reference + float(p) * log_likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )
    max_log_likelihood = float(np.max(log_likelihood))
    max_weight = weights[int(np.argmax(log_likelihood))]

    rows = []
    contribution_logs = []
    for source in bulk["occupation_rows"]:
        occupation = int(source["active_regular_outer_blocks"])
        inner_log = float(source["inner_log2_upper"]) * LOG2
        location_log = log_choose(outer_blocks, occupation)
        base_log = location_log + occupation * log_mass

        corrections = (
            occupation * log_moments / p_values
            + (1.0 - 1.0 / p_values) * inner_log
        )
        best_index = int(np.argmin(corrections))
        best_correction = float(corrections[best_index])
        best_p: float | str = float(p_values[best_index])

        infinity_correction = occupation * max_log_likelihood + inner_log
        if infinity_correction < best_correction:
            best_correction = infinity_correction
            best_p = "infinity"

        contribution = base_log + best_correction
        contribution_logs.append(contribution)
        rows.append(
            {
                "active_outer_blocks": occupation,
                "reference_inner_log2_upper": inner_log / LOG2,
                "best_holder_p": best_p,
                "holder_correction_log2": best_correction / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "pointwise_margin_bits": -contribution / LOG2,
            }
        )

    multi_log = float(logsumexp(np.asarray(contribution_logs)))
    combined_logs = [multi_log]
    one_active_log2 = None
    if one_active_path is not None:
        one_active = json.loads(one_active_path.read_text(encoding="utf-8"))
        one_active_log2 = float(one_active["aggregate_log2_upper"])
        combined_logs.append(one_active_log2 * LOG2)
    combined_log = float(logsumexp(np.asarray(combined_logs)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:20]

    return {
        "schema": "riffle-bulk-full-spectrum-holder-v1",
        "inputs": {
            "bulk_uniform_reference": str(bulk_path),
            "outer_spectrum": str(spectrum_path),
            "one_active": str(one_active_path) if one_active_path else None,
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "active_block_range": [
                int(rows[0]["active_outer_blocks"]),
                int(rows[-1]["active_outer_blocks"]),
            ],
            "holder_grid_size": int(len(p_values)),
            "holder_p_min": float(p_values[0]),
            "holder_p_max": float(p_values[-1]),
        },
        "spectrum_change_of_measure": {
            "log2_max_likelihood_ratio": max_log_likelihood / LOG2,
            "maximizing_weight": max_weight,
            "normalization_log2_error": normalization / LOG2,
        },
        "multi_active_log2_upper": multi_log / LOG2,
        "multi_active_margin_bits": -multi_log / LOG2,
        "one_active_log2_upper": one_active_log2,
        "combined_log2_upper": combined_log / LOG2,
        "combined_margin_bits": -combined_log / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Full-spectrum Holder change-of-measure diagnostic over the "
            "uniform-reference inner bounds. It avoids the single worst-density "
            "surcharge but inherits all assumptions and nearest-binary64 "
            "arithmetic of the supplied receipts; it is not an outward-rounded "
            "certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--one-active", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args.bulk, args.outer_spectrum, args.one_active)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"multi_active_margin,{result['multi_active_margin_bits']:.9f}")
    print(f"combined_margin,{result['combined_margin_bits']:.9f}")
    for row in result["dominant_rows"][:5]:
        print(
            "dominant,"
            f"{row['active_outer_blocks']},"
            f"p,{row['best_holder_p']},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
