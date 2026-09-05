#!/usr/bin/env python3
"""Dense-occupation diagnostic from selected BA distance and RandomStepConv.

Fix one repeated [240,120] constituent of minimum distance at least d. For
each active row, delete all but d ones after its independent coordinate
permutation. Adding input ones can only decrease the RandomStepConv output
moment. The retained d-subset is uniform. A Bernoulli envelope removes its
fixed-weight conditioning and makes the 240 regions independent.

Nearest binary64 log arithmetic is diagnostic, not an outward certificate.
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
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from analyze_riffle_striped_random_outer import (  # noqa: E402
    log_choose,
    log_matmul_batch,
    log_matrix_entries,
    log_matrix_power_moment,
    log_two_power_minus_one,
)
from evaluate_selected_ba_randomstepconv_g1_q2_64 import (  # noqa: E402
    B,
    D,
    K,
    L,
    LOG2,
    N,
    envelope,
    selected_spectrum,
    step_matrices,
)


DEFAULT_OUTPUT = (
    WORKSTREAM / "selected_ba240_randomstepconv_g1_s30_dense_q65_8832_d11.json"
)
TARGET_FULL_ROWS = (1 << 20) // K
TARGET_PARTIAL_BITS = (1 << 20) % K


def log_uniform_candidate_regions(
    *, zero: np.ndarray, candidate: np.ndarray, maximum_q: int
) -> np.ndarray:
    """Return log entrywise averages for every exact candidate count."""
    log_zero = log_matrix_entries(tuple(zero))
    log_candidate = log_matrix_entries(tuple(candidate))
    current = np.full((maximum_q + 1, 2, 2), -math.inf)
    current[0, 0, 0] = 0.0
    current[0, 1, 1] = 0.0

    for completed in range(L):
        next_maximum = min(completed + 1, maximum_q)
        updated = np.full_like(current, -math.inf)

        zero_products = log_matmul_batch(
            current[: next_maximum + 1], log_zero
        )
        degrees = np.arange(next_maximum + 1, dtype=np.float64)
        zero_weights = ((completed + 1.0) - degrees) / (completed + 1.0)
        positive = zero_weights > 0.0
        updated[: next_maximum + 1][positive] = (
            zero_products[positive] + np.log(zero_weights[positive])[:, None, None]
        )

        if next_maximum:
            candidate_products = log_matmul_batch(
                current[:next_maximum], log_candidate
            )
            selected_degrees = np.arange(1, next_maximum + 1, dtype=np.float64)
            candidate_terms = candidate_products + np.log(
                selected_degrees / (completed + 1.0)
            )[:, None, None]
            updated[1 : next_maximum + 1] = np.logaddexp(
                updated[1 : next_maximum + 1], candidate_terms
            )
        current = updated
        if completed and completed % 2048 == 0:
            print(f"region_dp,{completed},{L}", flush=True)
    return current


def logadd(*values: float) -> float:
    maximum = max(values)
    if maximum == -math.inf:
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def logdiffexp(high: float, low: float) -> float:
    """Return log(exp(high)-exp(low)) for high strictly above low."""
    if not high > low:
        raise ArithmeticError("positive global-activation mass was lost")
    return high + math.log(-math.expm1(low - high))


def nonzero_region_power_log_moment(
    log_region: np.ndarray, candidate_probability: float, candidates: int
) -> float:
    """Power one region while excluding the globally all-zero input."""
    log_all_zero_region = candidates * math.log1p(-candidate_probability)
    log_activate_zero = logdiffexp(
        float(log_region[0, 0]), log_all_zero_region
    )
    log_activate_live = float(log_region[0, 1])

    inactive = 0.0
    active_zero = -math.inf
    active_live = -math.inf
    for _ in range(B):
        next_zero = logadd(
            inactive + log_activate_zero,
            active_zero + float(log_region[0, 0]),
            active_live + float(log_region[1, 0]),
        )
        next_live = logadd(
            inactive + log_activate_live,
            active_zero + float(log_region[0, 1]),
            active_live + float(log_region[1, 1]),
        )
        inactive += log_all_zero_region
        active_zero, active_live = next_zero, next_live
    return logadd(active_zero, active_live)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    occupations = list(range(args.q_min, args.q_max + 1))
    best = np.full(args.q_max + 1, math.inf)
    witnesses: list[dict[str, float] | None] = [None] * (args.q_max + 1)
    message_log = log_two_power_minus_one(K)
    log_spectrum = None
    selection_metadata = None
    if args.outer_envelope == "selected-spectrum":
        log_spectrum, tail_log, good_lower, selection_factor = selected_spectrum(
            args.lower_weight, args.upper_weight
        )
        selection_metadata = {
            "allowed_outer_weights": [args.lower_weight, args.upper_weight],
            "tail_word_expectation_log2": tail_log / LOG2,
            "tail_free_probability_lower": good_lower,
            "common_spectrum_factor": selection_factor,
            "common_spectrum_factor_log2": math.log2(selection_factor),
        }
    elif args.outer_envelope == "information-set":
        occupations = list(
            range(args.q_min, min(args.q_max, TARGET_FULL_ROWS) + 1)
        )
        log_spectrum = np.full(B + 1, -math.inf)
        for weight in range(1, K + 1):
            log_spectrum[weight] = log_choose(K, weight)
        selection_metadata = {
            "information_set_coordinates": K,
            "effective_spectrum": "A_w=C(K,w) for 1<=w<=K",
            "justification": (
                "Restriction to any K-coordinate information set is a "
                "bijection. After the independent row-coordinate permutation, "
                "each retained weight-w support is uniform among B coordinates."
            ),
        }
    elif args.outer_envelope == "random-spectrum":
        occupations = list(
            range(args.q_min, min(args.q_max, TARGET_FULL_ROWS) + 1)
        )
        log_spectrum = np.full(B + 1, -math.inf)
        message_log = log_two_power_minus_one(K)
        for weight in range(1, B + 1):
            log_spectrum[weight] = (
                message_log + log_choose(B, weight) - B * LOG2
            )
        selection_metadata = {
            "effective_spectrum": (
                "A_w=(2^K-1) C(B,w)/2^B for 1<=w<=B"
            ),
            "status": "reference spectrum; fixed-code realization is separate",
        }

    p_values = np.arange(
        args.p_min, args.p_max + 0.5 * args.p_step, args.p_step
    )
    u_values = np.arange(
        args.grid_min,
        args.grid_max + 0.5 * args.grid_step,
        args.grid_step,
    )
    pair_count = len(p_values) * len(u_values)
    pair_index = 0
    for p_raw in p_values:
        p = float(p_raw)
        if log_spectrum is None:
            exact_weight_log_probability = (
                log_choose(B, args.minimum_weight)
                + args.minimum_weight * math.log(p)
                + (B - args.minimum_weight) * math.log1p(-p)
            )
            per_row_log = message_log - exact_weight_log_probability
        else:
            per_row_log, _ = envelope(log_spectrum, p)
        for u_raw in u_values:
            pair_index += 1
            u = float(u_raw)
            surprisal = math.exp(u)
            z = math.exp(-surprisal)
            zero, active = step_matrices(z, args.memory_bits)
            candidate = (1.0 - p) * zero + p * active
            log_regions = log_uniform_candidate_regions(
                zero=zero,
                candidate=candidate,
                maximum_q=args.q_max,
            )
            correction = D * surprisal
            for q in occupations:
                entries = log_regions[q].reshape(4)
                if args.exclude_global_zero:
                    inner_moment = nonzero_region_power_log_moment(
                        log_regions[q], p, q
                    )
                else:
                    scale = float(np.max(entries))
                    normalized = tuple(
                        math.exp(float(value - scale)) for value in entries
                    )
                    inner_moment = (
                        log_matrix_power_moment(normalized, B) + B * scale
                    )
                inner = inner_moment + correction
                inner = min(0.0, inner)
                row_universe = (
                    TARGET_FULL_ROWS
                    if args.outer_envelope in ("information-set", "random-spectrum")
                    else L
                )
                if q > row_universe:
                    continue
                partial_overcount = (
                    TARGET_PARTIAL_BITS * LOG2
                    if args.outer_envelope in ("information-set", "random-spectrum")
                    else 0.0
                )
                value = (
                    log_choose(row_universe, q)
                    + q * per_row_log
                    + partial_overcount
                    + inner
                )
                if value < best[q]:
                    best[q] = value
                    witnesses[q] = {
                        "candidate_probability": p,
                        "log_surprisal": u,
                        "z": z,
                        "per_row_envelope_log2": per_row_log / LOG2,
                        "inner_log2_upper": inner / LOG2,
                    }
            print(f"pair,{pair_index},{pair_count},p,{p:.6g},u,{u:.6g}", flush=True)

    rows = [
        {
            "active_outer_rows": q,
            "pointwise_log2_upper": float(best[q]) / LOG2,
            **(witnesses[q] or {}),
        }
        for q in occupations
    ]
    total = float(logsumexp(np.asarray([best[q] for q in occupations])))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "selected-ba240-randomstepconv-g1-dense-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim_scope": {
            "outer": (
                "one fixed [240,120] linear constituent of minimum distance "
                f"at least {args.minimum_weight}, repeated in every row"
            ),
            "occupations": [occupations[0], occupations[-1]],
            "log2_expected_bad_upper": total / LOG2,
            "margin_bits": -total / LOG2,
            "dominant": dominant,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "minimum_constituent_distance": args.minimum_weight,
            "outer_envelope": args.outer_envelope,
            "globally_zero_reference_input_excluded": args.exclude_global_zero,
            "target_full_rows": TARGET_FULL_ROWS,
            "target_partial_row_message_bits": TARGET_PARTIAL_BITS,
        },
        "constituent_selection": selection_metadata,
        "bound": {
            "monotonicity": (
                "Deleting input ones can only increase the tilted output moment."
            ),
            "deletion": (
                "After the uniform row-coordinate permutation, retain a uniform "
                "minimum_weight-subset of each nonzero constituent support."
            ),
            "bernoulli_envelope": (
                "A uniform minimum_weight-subset equals iid Bernoulli-p region "
                "incidences conditioned on exact weight minimum_weight. Dropping "
                "that conditioning costs the reciprocal binomial mass per row."
            ),
            "region_transfer": (
                "The coefficient recurrence exactly averages the Q candidate "
                "positions over the uniform region permutation."
            ),
            "global_zero_exclusion": (
                "A three-state region recurrence omits the globally all-zero "
                "reference input without alternating subtraction. Reference "
                "tuples with some, but not all, zero rows remain included."
            ),
        },
        "grid": {
            "candidate_probability": [args.p_min, args.p_max, args.p_step],
            "log_surprisal": [args.grid_min, args.grid_max, args.grid_step],
            "evaluated_pairs": pair_count,
        },
        "occupation_rows": rows,
        "limitations": [
            "The calculation uses nearest binary64 log arithmetic.",
            "The witness grid is diagnostic and may be sharpened.",
            "The constituent-selection or construction proof is separate.",
            "Occupations below the displayed range are not covered.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument(
        "--outer-envelope",
        choices=(
            "random-spectrum",
            "information-set",
            "selected-spectrum",
            "minimum-distance",
        ),
        default="selected-spectrum",
    )
    parser.add_argument("--minimum-weight", type=int, default=25)
    parser.add_argument("--lower-weight", type=int, default=25)
    parser.add_argument("--upper-weight", type=int, default=215)
    parser.add_argument("--q-min", type=int, default=65)
    parser.add_argument("--q-max", type=int, default=L)
    parser.add_argument("--p-min", type=float, default=25 / 240)
    parser.add_argument("--p-max", type=float, default=25 / 240)
    parser.add_argument("--p-step", type=float, default=0.01)
    parser.add_argument("--grid-min", type=float, default=-4.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--include-global-zero",
        dest="exclude_global_zero",
        action="store_false",
        help="retain the old loose all-zero reference contribution",
    )
    parser.set_defaults(exclude_global_zero=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim_scope"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
