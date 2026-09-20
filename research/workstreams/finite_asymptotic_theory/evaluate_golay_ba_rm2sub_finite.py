#!/usr/bin/env python3
"""Finite binary64 diagnostic for Golay--BA-3/RM2Sub-S19.

The BA spectrum is the exact ensemble expectation in combinatorial form,
evaluated in binary64 log space.  The one-active inner calculation uses the
finite RM2Sub entrywise transfer.  The result is diagnostic, not an outward
certificate for a reused sampled BA realization.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from analyze_golay_ba_rm2sub_joint import expected_ba_log_spectrum  # noqa: E402
from evaluate_random_outer_rm2sub_one_active import (  # noqa: E402
    DEFAULT_ACTIVATION,
    DEFAULT_LIVE_SPECTRUM,
    load_nonactivation,
    log_hamming_ball,
    tilt_grid,
)
from analyze_riffle_bitshuffle_splitstate_regular_bulk import (  # noqa: E402
    LOG2,
    load_nonzero_spectrum,
    log_choose,
    splitstate_impulse_matrices,
)
from analyze_riffle_splitstate_preaddmul_one_active import (  # noqa: E402
    marked_region,
    weight_conditioned_log_moments,
)


def logadd(values: list[float]) -> float:
    if not values:
        return -math.inf
    return float(logsumexp(np.asarray(values, dtype=np.float64)))


def one_active(
    *,
    spectrum: np.ndarray,
    outer_bits: int,
    outer_blocks: int,
    distance: int,
    step_bits: int,
    state_bits: int,
    constituent_distance: int,
    activation: Path,
    live_spectrum_path: Path,
    tilt_minimum: float,
    tilt_maximum: float,
    tilt_spacing: float,
) -> dict[str, object]:
    epochs_per_region = outer_blocks // step_bits
    nonactivation = load_nonactivation(activation, step_bits)
    live_spectrum = load_nonzero_spectrum(
        live_spectrum_path, step_bits, state_bits
    )
    best_by_weight = np.full(outer_bits + 1, math.inf, dtype=np.float64)
    best_tilt = np.full(outer_bits + 1, math.nan, dtype=np.float64)

    for log_surprisal in tilt_grid(
        tilt_minimum, tilt_maximum, tilt_spacing
    ):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        impulse = splitstate_impulse_matrices(
            z=z,
            step_bits=step_bits,
            state_bits=state_bits,
            constituent_distance=constituent_distance,
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
        candidate = moments + distance * surprisal
        improved = candidate < best_by_weight
        best_by_weight[improved] = candidate[improved]
        best_tilt[improved] = log_surprisal

    rows = []
    terms = []
    for weight in range(1, outer_bits + 1):
        if not math.isfinite(float(spectrum[weight])):
            continue
        inner = min(0.0, float(best_by_weight[weight]))
        term = math.log(outer_blocks) + float(spectrum[weight]) + inner
        terms.append(term)
        rows.append(
            {
                "outer_weight": weight,
                "log2_expected_ba_multiplicity": (
                    float(spectrum[weight]) / LOG2
                ),
                "best_log_surprisal": float(best_tilt[weight]),
                "inner_log2_upper": inner / LOG2,
                "pointwise_log2_upper": term / LOG2,
                "pointwise_margin_bits": -term / LOG2,
            }
        )
    aggregate = logadd(terms)
    dominant = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "aggregate_log2_upper": aggregate / LOG2,
        "aggregate_margin_bits": -aggregate / LOG2,
        "dominant_rows": dominant,
        "weight_rows": rows,
        "method": (
            "Expected BA spectrum over its two interleavers, followed by the "
            "finite two-state RM2Sub entrywise transfer and per-weight tilt."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=216)
    parser.add_argument("--outer-blocks", type=int, default=9728)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=19)
    parser.add_argument("--constituent-distance", type=int, default=48)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument(
        "--live-spectrum", type=Path, default=DEFAULT_LIVE_SPECTRUM
    )
    parser.add_argument("--tilt-minimum", type=float, default=-10.0)
    parser.add_argument("--tilt-maximum", type=float, default=-6.0)
    parser.add_argument("--tilt-spacing", type=float, default=0.01)
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKSTREAM / "golay_ba3_rm2sub_finite_k20_d11.json",
    )
    parser.add_argument(
        "--conditioned-spectrum-output",
        type=Path,
        default=WORKSTREAM / "golay_ba3_B216_conditioned_spectrum_upper.json",
    )
    args = parser.parse_args()

    if args.outer_bits % 24:
        raise ValueError("outer bits must be divisible by 24")
    if args.outer_blocks % args.step_bits:
        raise ValueError("step bits must divide the outer-block count")
    parent_output_bits = args.outer_bits * args.outer_blocks
    parent_dimension = parent_output_bits // 2
    if args.message_bits > parent_dimension:
        raise ValueError("requested message dimension exceeds the parent")
    distance = math.floor(args.relative_distance * parent_output_bits)

    spectrum = expected_ba_log_spectrum(args.outer_bits)
    lower = math.ceil(0.104 * args.outer_bits)
    upper = math.floor(0.896 * args.outer_bits)
    low_tail = logadd(
        [float(value) for value in spectrum[1:lower] if math.isfinite(value)]
    )
    high_tail = logadd(
        [
            float(value)
            for value in spectrum[upper + 1 :]
            if math.isfinite(value)
        ]
    )
    tail = float(np.logaddexp(low_tail, high_tail))
    good_probability_lower = 1.0 - min(1.0, math.exp(tail))
    if good_probability_lower <= 0.0:
        raise ValueError("the Markov lower bound for the good event vanished")
    conditioning_cost = -math.log(good_probability_lower)

    active = one_active(
        spectrum=spectrum,
        outer_bits=args.outer_bits,
        outer_blocks=args.outer_blocks,
        distance=distance,
        step_bits=args.step_bits,
        state_bits=args.state_bits,
        constituent_distance=args.constituent_distance,
        activation=args.activation,
        live_spectrum_path=args.live_spectrum,
        tilt_minimum=args.tilt_minimum,
        tilt_maximum=args.tilt_maximum,
        tilt_spacing=args.tilt_spacing,
    )
    conditional_one_active_terms = [
        float(row["pointwise_log2_upper"]) * LOG2
        for row in active["weight_rows"]
        if lower <= int(row["outer_weight"]) <= upper
    ]
    conditional_one_active_log = (
        logadd(conditional_one_active_terms) - math.log(good_probability_lower)
    )

    random_margin = (
        parent_output_bits
        - args.message_bits
        - log_hamming_ball(parent_output_bits, distance) / LOG2
    )
    result = {
        "schema": "golay-ba3-rm2sub-finite-k20-v1",
        "status": "BINARY64_FINITE_DIAGNOSTIC",
        "parameters": {
            "requested_message_bits": args.message_bits,
            "parent_output_bits": parent_output_bits,
            "parent_dimension": parent_dimension,
            "shortened_input_coordinates": (
                parent_dimension - args.message_bits
            ),
            "effective_rate": args.message_bits / parent_output_bits,
            "outer_bits": args.outer_bits,
            "outer_dimension": args.outer_bits // 2,
            "outer_blocks": args.outer_blocks,
            "golay_blocks_per_outer": args.outer_bits // 24,
            "step_bits": args.step_bits,
            "inner_epochs": parent_output_bits // args.step_bits,
            "relative_distance": args.relative_distance,
            "distance": distance,
        },
        "outer_selection": {
            "asymptotic_weight_window": [lower, upper],
            "expected_low_tail_words_log2": low_tail / LOG2,
            "expected_high_tail_words_log2": high_tail / LOG2,
            "expected_total_tail_words_log2": tail / LOG2,
            "markov_no_tail_failure_probability_upper": min(1.0, math.exp(tail)),
            "good_event_probability_lower": good_probability_lower,
            "expected_rejection_sampling_trials_upper": 1.0 / good_probability_lower,
            "conditional_expected_spectrum_cost_bits": conditioning_cost / LOG2,
        },
        "one_active": active,
        "one_active_conditioned_tail_free": {
            "aggregate_log2_upper": conditional_one_active_log / LOG2,
            "aggregate_margin_bits": -conditional_one_active_log / LOG2,
            "probability_space": (
                "One BA draw conditioned on the recorded tail-free event, "
                "one independent coordinate permutation, the region "
                "permutations, and the RM2Sub multipliers."
            ),
            "derivation": (
                "Delete the forbidden outer weights and divide each remaining "
                "expected multiplicity by the good-event probability lower bound."
            ),
        },
        "ideal_random_rate_comparator": {
            "finite_hamming_ball_margin_bits": random_margin,
            "interpretation": (
                "First-moment margin for a uniform random binary linear code "
                "with the shortened dimension and parent output length."
            ),
        },
        "spectrum": [
            {
                "weight": weight,
                "log2_expected_multiplicity": (
                    None
                    if not math.isfinite(float(spectrum[weight]))
                    else float(spectrum[weight]) / LOG2
                ),
            }
            for weight in range(args.outer_bits + 1)
        ],
        "limitations": [
            "Nearest-binary64 arithmetic is not outward rounded.",
            "The one-active result averages over the sampled BA outer.",
            "A reused BA outer requires higher spectrum moments for Q>=2.",
            "No finite union over all occupations is included.",
            "The shortened subspace can only decrease the bad-codeword count; its exact block profile is not modeled.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    conditioned_spectrum = {
        "schema": "golay-ba3-conditioned-expected-spectrum-upper-v1",
        "status": "BINARY64_EXPECTATION_UPPER_FROM_MARKOV_CONDITIONING",
        "outer_bits": args.outer_bits,
        "outer_dimension": args.outer_bits // 2,
        "conditioning_event": (
            f"no nonzero BA word outside weights {lower}..{upper}"
        ),
        "good_event_probability_lower": good_probability_lower,
        "scope": (
            "For one BA draw conditioned on the stated event, each central "
            "expected multiplicity is at most its unconditional expectation "
            "divided by the good-event probability lower bound. Products of "
            "these expectations apply only when BA draws are independent "
            "across outer rows."
        ),
        "spectrum": [
            {
                "weight": weight,
                "log2_expected_multiplicity": (
                    0.0
                    if weight == 0
                    else (
                        float(spectrum[weight] + conditioning_cost) / LOG2
                        if lower <= weight <= upper
                        and math.isfinite(float(spectrum[weight]))
                        else None
                    )
                ),
            }
            for weight in range(args.outer_bits + 1)
        ],
    }
    args.conditioned_spectrum_output.write_text(
        json.dumps(conditioned_spectrum, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "conditioned_spectrum_output": str(
                    args.conditioned_spectrum_output
                ),
                "one_active_margin_bits": active["aggregate_margin_bits"],
                "conditional_one_active_margin_bits": (
                    -conditional_one_active_log / LOG2
                ),
                "tail_failure_probability_upper": result[
                    "outer_selection"
                ]["markov_no_tail_failure_probability_upper"],
                "ideal_random_margin_bits": random_margin,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
