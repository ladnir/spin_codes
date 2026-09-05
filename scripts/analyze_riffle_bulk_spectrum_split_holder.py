#!/usr/bin/env python3
"""Two-piece full-spectrum Holder bound for the Riffle bulk calculation.

Split the outer likelihood ratio into disjoint tail and body supports. For each
active-block occupation r, expand the product by the exact number j of tail
words. Apply Holder separately to every j term, optimize its exponent, and sum
the resulting bounds. This prevents a rare low-weight shell from setting the
density surcharge for every active outer block.
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


def holder_grid(points: int) -> np.ndarray:
    if points < 3:
        raise ValueError("Holder grid needs at least three points")
    return 1.0 + np.exp2(np.linspace(-16.0, 16.0, points))


def class_log_moments(
    p_values: np.ndarray,
    log_reference: np.ndarray,
    log_likelihood: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    if not np.any(selected):
        return np.full_like(p_values, -math.inf)
    reference = log_reference[selected]
    likelihood = log_likelihood[selected]
    return np.asarray(
        [
            logsumexp(reference + float(p) * likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )


def evaluate(
    bulk_path: Path,
    spectrum_path: Path,
    one_active_path: Path | None,
    tail_endpoint_width: int,
    grid_points: int,
) -> dict[str, object]:
    bulk = json.loads(bulk_path.read_text(encoding="utf-8"))
    parameters = bulk["parameters"]
    outer_bits = int(parameters["outer_bits"])
    dimension = outer_bits // 2
    outer_blocks = int(parameters["outer_blocks"])
    if not 1 <= tail_endpoint_width < outer_bits // 2:
        raise ValueError("tail endpoint width must be below half the outer block")
    spectrum = load_outer_spectrum_logs(spectrum_path, outer_bits)

    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    weights = []
    log_reference_probabilities = []
    log_likelihoods = []
    log_multiplicities = []
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
        log_multiplicities.append(log_multiplicity)

    weights_array = np.asarray(weights, dtype=np.int64)
    log_reference = np.asarray(log_reference_probabilities, dtype=np.float64)
    log_likelihood = np.asarray(log_likelihoods, dtype=np.float64)
    log_multiplicity = np.asarray(log_multiplicities, dtype=np.float64)
    normalization = float(logsumexp(log_reference + log_likelihood))
    if abs(normalization) > 1e-8:
        raise AssertionError(
            f"outer likelihood does not normalize: log expectation={normalization}"
        )

    tail = np.minimum(weights_array, outer_bits - weights_array) <= tail_endpoint_width
    body = ~tail
    if not np.any(tail) or not np.any(body):
        raise ValueError("tail cutoff must leave both classes nonempty")
    tail_log_count = float(logsumexp(log_multiplicity[tail]))
    body_log_count = float(logsumexp(log_multiplicity[body]))

    p_values = holder_grid(grid_points)
    inverse_p = 1.0 / p_values
    event_exponents = 1.0 - inverse_p
    tail_moments = class_log_moments(
        p_values, log_reference, log_likelihood, tail
    )
    body_moments = class_log_moments(
        p_values, log_reference, log_likelihood, body
    )

    rows = []
    occupation_logs = []
    for source in bulk["occupation_rows"]:
        occupation = int(source["active_regular_outer_blocks"])
        inner_log = float(source["inner_log2_upper"]) * LOG2
        tail_counts = np.arange(occupation + 1, dtype=np.float64)
        body_counts = occupation - tail_counts
        log_class_choices = np.asarray(
            [log_choose(occupation, int(j)) for j in tail_counts],
            dtype=np.float64,
        )

        # Each column is one Holder exponent; each row is an exact number of
        # tail outer words. The minimum is taken independently for every row.
        corrections = (
            tail_counts[:, None] * tail_moments[None, :] * inverse_p[None, :]
            + body_counts[:, None] * body_moments[None, :] * inverse_p[None, :]
            + inner_log * event_exponents[None, :]
        )
        best_indices = np.argmin(corrections, axis=1)
        best_corrections = corrections[
            np.arange(occupation + 1), best_indices
        ]
        class_logs = log_class_choices + best_corrections
        class_sum = float(logsumexp(class_logs))
        location_log = log_choose(outer_blocks, occupation)
        contribution = location_log + occupation * log_mass + class_sum
        occupation_logs.append(contribution)

        dominant_tail_count = int(np.argmax(class_logs))
        rows.append(
            {
                "active_outer_blocks": occupation,
                "reference_inner_log2_upper": inner_log / LOG2,
                "dominant_tail_blocks": dominant_tail_count,
                "dominant_holder_p": float(
                    p_values[best_indices[dominant_tail_count]]
                ),
                "class_sum_log2": class_sum / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
                "pointwise_margin_bits": -contribution / LOG2,
            }
        )

    multi_log = float(logsumexp(np.asarray(occupation_logs)))
    combined_logs = [multi_log]
    one_active_log2 = None
    if one_active_path is not None:
        one_active = json.loads(one_active_path.read_text(encoding="utf-8"))
        one_active_log2 = float(one_active["aggregate_log2_upper"])
        combined_logs.append(one_active_log2 * LOG2)
    combined_log = float(logsumexp(np.asarray(combined_logs)))
    dominant_rows = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:20]

    tail_probability_log2 = (tail_log_count - log_mass) / LOG2
    body_max_index = int(np.argmax(np.where(body, log_likelihood, -math.inf)))
    tail_max_index = int(np.argmax(np.where(tail, log_likelihood, -math.inf)))
    return {
        "schema": "riffle-bulk-spectrum-split-holder-v1",
        "inputs": {
            "bulk_uniform_reference": str(bulk_path),
            "outer_spectrum": str(spectrum_path),
            "one_active": str(one_active_path) if one_active_path else None,
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "tail_endpoint_width": tail_endpoint_width,
            "active_block_range": [
                int(rows[0]["active_outer_blocks"]),
                int(rows[-1]["active_outer_blocks"]),
            ],
            "holder_grid_size": int(len(p_values)),
            "holder_p_min": float(p_values[0]),
            "holder_p_max": float(p_values[-1]),
        },
        "spectrum_split": {
            "tail_log2_expected_count": tail_log_count / LOG2,
            "tail_log2_probability_given_nonzero_message": tail_probability_log2,
            "body_log2_expected_count": body_log_count / LOG2,
            "tail_max_log2_likelihood_ratio": float(
                log_likelihood[tail_max_index] / LOG2
            ),
            "tail_maximizing_weight": int(weights_array[tail_max_index]),
            "body_max_log2_likelihood_ratio": float(
                log_likelihood[body_max_index] / LOG2
            ),
            "body_maximizing_weight": int(weights_array[body_max_index]),
            "normalization_log2_error": normalization / LOG2,
        },
        "multi_active_log2_upper": multi_log / LOG2,
        "multi_active_margin_bits": -multi_log / LOG2,
        "one_active_log2_upper": one_active_log2,
        "combined_log2_upper": combined_log / LOG2,
        "combined_margin_bits": -combined_log / LOG2,
        "dominant_rows": dominant_rows,
        "occupation_rows": rows,
        "scope": (
            "Two-piece full-spectrum Holder diagnostic. Every exact tail-count "
            "class is bounded separately before summation. It inherits the "
            "uniform-reference inner assumptions and nearest-binary64 arithmetic "
            "of the source receipts; it is not an outward-rounded certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--one-active", type=Path)
    parser.add_argument("--tail-endpoint-width", type=int, required=True)
    parser.add_argument("--holder-grid-points", type=int, default=129)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        args.bulk,
        args.outer_spectrum,
        args.one_active,
        args.tail_endpoint_width,
        args.holder_grid_points,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    split = result["spectrum_split"]
    print(
        "split,"
        f"tail_count_log2,{split['tail_log2_expected_count']:.9f},"
        f"tail_probability_log2,{split['tail_log2_probability_given_nonzero_message']:.9f},"
        f"body_eta_log2,{split['body_max_log2_likelihood_ratio']:.9f}"
    )
    print(f"multi_active_margin,{result['multi_active_margin_bits']:.9f}")
    print(f"combined_margin,{result['combined_margin_bits']:.9f}")
    for row in result["dominant_rows"][:5]:
        print(
            "dominant,"
            f"{row['active_outer_blocks']},"
            f"tail_blocks,{row['dominant_tail_blocks']},"
            f"p,{row['dominant_holder_p']:.9f},"
            f"margin,{row['pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
