#!/usr/bin/env python3
"""Fast Holder bound for an EBCH body and its unique all-one tail word.

The per-block permutation makes a word of weight w uniform on its Hamming
shell.  Every nonzero outer word except the all-one word is placed in the
body and enters through its complete weight spectrum.  The unique all-one
word is the tail.  At a fixed Holder exponent, the sum over every possible
number of body and tail blocks is a binomial middle sum, so evaluating all
outer occupations costs O(M P), rather than O(M^2 P).

The supplied bulk receipt must contain uniform-reference inner bounds indexed
by the total number of active outer blocks.  Occupation one is taken from the
separate exact one-active receipt.
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


def log_expm1(value: float) -> float:
    """Return log(exp(value)-1) accurately for positive value."""
    if value <= 0.0:
        return -math.inf
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def log_subtract(left: float, right: float) -> float:
    """Return log(exp(left)-exp(right)) for left >= right."""
    if right == -math.inf:
        return left
    if right >= left:
        if right - left < 1e-11:
            return -math.inf
        raise ArithmeticError("negative logarithmic subtraction")
    return left + math.log1p(-math.exp(right - left))


def log_binomial_middle(active: int, body_log: float, tail_log: float) -> float:
    """Log sum C(active,k) body^(active-k) tail^k, 1 <= k < active."""
    if active < 2:
        return -math.inf
    high = max(body_log, tail_log)
    low = min(body_log, tail_log)
    log_ratio = low - high
    # Factor out exp(active*high).  The remaining middle sum is
    # (1+z)^active - 1 - z^active, where z=exp(log_ratio).
    z = math.exp(log_ratio)
    all_non_high = active * math.log1p(z)
    log_all_non_high = log_expm1(all_non_high)
    low_endpoint = active * log_ratio
    middle_factor = log_subtract(log_all_non_high, low_endpoint)
    return active * high + middle_factor


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    bulks = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in args.bulk
    ]
    if not bulks:
        raise ValueError("at least one bulk receipt is required")
    bulk = bulks[0]
    one_active = json.loads(args.one_active.read_text(encoding="utf-8"))
    parameters = bulk["parameters"]
    common_parameters = dict(parameters)
    for grid_parameter in ("active_blocks", "log_surprisals"):
        common_parameters.pop(grid_parameter, None)
    for candidate in bulks[1:]:
        candidate_parameters = dict(candidate["parameters"])
        for grid_parameter in ("active_blocks", "log_surprisals"):
            candidate_parameters.pop(grid_parameter, None)
        if candidate_parameters != common_parameters:
            raise ValueError("bulk receipts do not have identical parameters")
    outer_bits = int(parameters["outer_bits"])
    outer_blocks = int(parameters["outer_blocks"])
    dimension = outer_bits // 2
    spectrum = load_outer_spectrum_logs(args.outer_spectrum, outer_bits)

    if args.no_tail:
        body_weight_range = range(1, outer_bits + 1)
    else:
        all_one_log_count = float(spectrum[outer_bits])
        if not math.isfinite(all_one_log_count) or abs(all_one_log_count) > 1e-12:
            raise ValueError("the all-one tail must have multiplicity exactly one")
        body_weight_range = range(1, outer_bits)

    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    body_reference: list[float] = []
    body_likelihood: list[float] = []
    body_weights: list[int] = []
    for weight in body_weight_range:
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(outer_bits, weight)
        body_reference.append(shell - outer_bits * LOG2)
        body_likelihood.append(
            multiplicity - shell - reference_log_density
        )
        body_weights.append(weight)
    if not body_weights:
        raise ValueError("outer spectrum has no non-tail body")

    envelope_by_probability = {
        float(row["probability"]): float(row["log2_mass"]) * LOG2
        for row in bulk["envelope"]["bernoulli_envelopes"]
    }
    inner_by_occupation: dict[int, float] = {}
    for candidate in bulks:
        for row in candidate["occupation_rows"]:
            occupation = int(row["active_regular_outer_blocks"])
            probability = float(row["best_candidate_probability"])
            try:
                envelope_log = envelope_by_probability[probability]
            except KeyError as error:
                raise ValueError(
                    f"missing Bernoulli envelope for probability {probability}"
                ) from error
            # The bulk receipt stores occupation*envelope plus the moment
            # under the selected Bernoulli reference.  Holder reweights from
            # that reference using the exact body spectrum, so remove the
            # envelope before applying the likelihood-ratio moment.
            value = (
                float(row["inner_log2_upper"]) * LOG2
                - occupation * envelope_log
            )
            previous = inner_by_occupation.get(occupation)
            if previous is not None and abs(previous - value) > 1e-10:
                raise ValueError(
                    f"bulk receipts disagree at occupation {occupation}"
                )
            inner_by_occupation[occupation] = value
    if any(a not in inner_by_occupation for a in range(2, outer_blocks + 1)):
        raise ValueError("bulk receipt must contain occupations 2..outer_blocks")

    p_values = 1.0 + np.exp2(
        np.linspace(args.holder_log2_min, args.holder_log2_max,
                    args.holder_grid_points)
    )
    inverse_p = 1.0 / p_values
    q_values = 1.0 - inverse_p
    reference = np.asarray(body_reference, dtype=np.float64)
    likelihood = np.asarray(body_likelihood, dtype=np.float64)
    body_log_moments = np.asarray(
        [
            logsumexp(reference + float(p) * likelihood)
            for p in p_values
        ],
        dtype=np.float64,
    )
    # One body block contributes the active-message multiplicity log_mass and
    # its Holder spectrum moment.  One all-one tail block contributes only its
    # exact shell-conditioning penalty; its multiplicity is one.
    body_terms = log_mass + body_log_moments * inverse_p
    tail_terms = q_values * outer_bits * LOG2

    body_rows: list[dict[str, object]] = []
    mixed_rows: list[dict[str, object]] = []
    body_logs: list[float] = []
    mixed_logs: list[float] = []
    for active in range(2, outer_blocks + 1):
        inner_log = inner_by_occupation[active]
        location = log_choose(outer_blocks, active)

        body_candidates = (
            location + active * body_terms + q_values * inner_log
        )
        body_index = int(np.argmin(body_candidates))
        body_value = float(body_candidates[body_index])
        body_logs.append(body_value)
        body_rows.append(
            {
                "active_blocks": active,
                "best_holder_p": float(p_values[body_index]),
                "reference_inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": body_value / LOG2,
                "pointwise_margin_bits": -body_value / LOG2,
            }
        )

        if not args.no_tail:
            mixed_candidates = np.empty_like(p_values)
            for index, (body_term, tail_term, q) in enumerate(
                zip(body_terms, tail_terms, q_values)
            ):
                mixed_candidates[index] = (
                    location
                    + q * inner_log
                    + log_binomial_middle(
                        active, float(body_term), float(tail_term)
                    )
                )
            mixed_index = int(np.argmin(mixed_candidates))
            mixed_value = float(mixed_candidates[mixed_index])
            mixed_logs.append(mixed_value)
            mixed_rows.append(
                {
                    "active_blocks": active,
                    "best_holder_p": float(p_values[mixed_index]),
                    "reference_inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": mixed_value / LOG2,
                    "pointwise_margin_bits": -mixed_value / LOG2,
                }
            )

    one_body_logs: list[float] = []
    one_tail_logs: list[float] = []
    for row in one_active["weight_rows"]:
        value = float(row["pointwise_log2_upper"]) * LOG2
        if not args.no_tail and int(row["outer_weight"]) == outer_bits:
            one_tail_logs.append(value)
        else:
            one_body_logs.append(value)
    one_body_log = float(logsumexp(np.asarray(one_body_logs)))
    one_tail_log = (
        float(logsumexp(np.asarray(one_tail_logs)))
        if one_tail_logs else -math.inf
    )
    body_aggregate = float(
        logsumexp(np.asarray([one_body_log, *body_logs]))
    )
    mixed_aggregate = (
        float(logsumexp(np.asarray(mixed_logs)))
        if mixed_logs else -math.inf
    )

    dominant_body = sorted(
        body_rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    dominant_mixed = sorted(
        mixed_rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    return {
        "schema": "riffle-allone-tail-body-holder-v1",
        "inputs": {
            "bulk_uniform_reference": [str(path) for path in args.bulk],
            "outer_spectrum": str(args.outer_spectrum),
            "one_active": str(args.one_active),
        },
        "parameters": {
            "outer_bits": outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "holder_grid_points": args.holder_grid_points,
            "holder_log2_min": args.holder_log2_min,
            "holder_log2_max": args.holder_log2_max,
        },
        "split": {
            "body_minimum_weight": min(body_weights),
            "body_maximum_weight": max(body_weights),
            "body_weight_count": len(body_weights),
            "tail_weight": None if args.no_tail else outer_bits,
            "tail_multiplicity": 0 if args.no_tail else 1,
        },
        "one_active": {
            "body_log2_upper": one_body_log / LOG2,
            "body_margin_bits": -one_body_log / LOG2,
            "tail_log2_upper": one_tail_log / LOG2,
            "tail_margin_bits": -one_tail_log / LOG2,
        },
        "body_only_log2_upper": body_aggregate / LOG2,
        "body_only_margin_bits": -body_aggregate / LOG2,
        "mixed_log2_upper": mixed_aggregate / LOG2,
        "mixed_margin_bits": -mixed_aggregate / LOG2,
        "dominant_body_rows": dominant_body,
        "dominant_mixed_rows": dominant_mixed,
        "body_occupation_rows": body_rows,
        "mixed_occupation_rows": mixed_rows,
        "scope": (
            "Rigorous Holder bounds for body-only occupations and every mixed "
            "occupation containing at least one non-all-one body word and at "
            "least one all-one tail word. The complete supplied outer weight "
            "spectrum is used for the body. The all-one shell and multiplicity "
            "are exact. A single Holder exponent is optimized after summing "
            "all body/tail counts at each total occupation, which is valid but "
            "can be weaker than optimizing each count separately. Tail-only "
            "occupations of two or more blocks are outside this receipt. "
            "Inherits the supplied inner-model assumptions and uses nearest-"
            "binary64 arithmetic without outward rounding."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, nargs="+", required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--one-active", type=Path, required=True)
    parser.add_argument("--holder-grid-points", type=int, default=513)
    parser.add_argument("--holder-log2-min", type=float, default=-16.0)
    parser.add_argument("--holder-log2-max", type=float, default=16.0)
    parser.add_argument(
        "--no-tail",
        action="store_true",
        help="place every nonzero spectrum shell in the Holder body",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"body_margin_bits,{result['body_only_margin_bits']:.9f}")
    print(f"mixed_margin_bits,{result['mixed_margin_bits']:.9f}")
    for kind in ("body", "mixed"):
        rows = result[f"dominant_{kind}_rows"]
        for row in rows[:3]:
            print(
                f"{kind},active,{row['active_blocks']},"
                f"p,{row['best_holder_p']:.9f},"
                f"margin,{row['pointwise_margin_bits']:.9f}"
            )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
