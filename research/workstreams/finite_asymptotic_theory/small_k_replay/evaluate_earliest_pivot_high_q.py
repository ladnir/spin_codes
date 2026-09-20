#!/usr/bin/env python3
"""Higher-occupation bound from the earliest guaranteed outer pivot.

For each spectrum band, delete a word to the band's lower weight.  After the
uniform row permutation, the retained support is a uniform fixed-weight
subset of the B coordinate regions.  Across all active rows, keep only the
earliest retained input bit.  RandomStepConv monotonicity makes both deletion
steps conservative.

The exact distribution of the earliest occupied coordinate region follows
from hypergeometric survival probabilities.  The proof then moves the kept
bit to the final position of that region.  The one-pulse tilted moment is
monotone under this delay.  The script checks that numerical monotonicity for
every witness but does not outward-certify it.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as inner  # noqa: E402
from small_k_replay.evaluate_banded_spectrum_high_q import (  # noqa: E402
    compositions,
    log_multinomial,
    parse_bands,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "earliest_pivot_high_q_diagnostic.json"


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_difference(first: float, second: float) -> float:
    """Return log(exp(first)-exp(second)) for first>second."""
    if second == -math.inf:
        return first
    if not second < first:
        raise ArithmeticError("nonpositive survival difference")
    return first + math.log1p(-math.exp(second - first))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constituent", choices=sorted(CONSTITUENTS), default="rm49")
    parser.add_argument("--message-exponent", type=int, default=13)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--band-starts", type=int, nargs="+", default=[32, 96, 192, 256])
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 0.7, 0.8, 1.0],
    )
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    constituent = CONSTITUENTS[args.constituent]
    spectrum = load_spectrum(constituent)
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("constituent dimension must divide total message bits")
    outer_rows = message_bits // constituent.dimension
    if not 1 <= args.minimum_occupation <= outer_rows:
        parser.error("minimum occupation must lie in [1,outer_rows]")
    output_bits = 2 * message_bits
    distance_cutoff = (
        args.distance_numerator * output_bits + args.distance_denominator - 1
    ) // args.distance_denominator
    block_bits = constituent.block_bits
    bands = parse_bands(args.band_starts, block_bits)

    active_bands = []
    band_rows = []
    for lower, upper in bands:
        mass = sum(
            count for weight, count in spectrum.items()
            if weight and lower <= weight <= upper
        )
        if not mass:
            continue
        active_bands.append((lower, upper, mass))
        band_rows.append({
            "lower_weight": lower,
            "upper_weight": upper,
            "exact_nonzero_mass": str(mass),
            "exact_nonzero_mass_log2": math.log2(mass),
        })
    if sum(item[2] for item in active_bands) != (1 << constituent.dimension) - 1:
        raise AssertionError("bands do not cover every nonzero codeword")

    all_compositions = [
        composition
        for occupation in range(args.minimum_occupation, outer_rows + 1)
        for composition in compositions(occupation, len(active_bands))
    ]
    # survival[type, r] is the probability that one row's retained support
    # avoids the first r coordinate regions.
    survival = np.full((len(active_bands), block_bits + 1), -math.inf)
    for index, (lower, _upper, _mass) in enumerate(active_bands):
        denominator = log_choose(block_bits, lower)
        for prefix in range(block_bits + 1):
            survival[index, prefix] = log_choose(
                block_bits - prefix, lower
            ) - denominator

    best = {composition: math.inf for composition in all_compositions}
    witnesses = {composition: math.nan for composition in all_compositions}
    monotonicity_checks = []
    moment_rows = []
    surprisals = []
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        surprisals.append(surprisal)
        z = math.exp(-surprisal)
        zero, active = inner.step_matrices(z, args.memory_bits)
        log_zero = inner.log_entries(zero)
        log_active = inner.log_entries(active)
        log_moments = np.empty(block_bits, dtype=np.float64)
        for region in range(block_bits):
            pulse_time = (region + 1) * outer_rows - 1
            suffix = output_bits - pulse_time - 1
            complete = inner.log_matmul(log_active, inner.log_power(log_zero, suffix))
            log_moments[region] = np.logaddexp(complete[0, 0], complete[0, 1])
        monotone = bool(np.all(log_moments[1:] >= log_moments[:-1] - 1e-12))
        if not monotone:
            raise AssertionError("one-pulse moment is not monotone under delay")
        monotonicity_checks.append({
            "log_surprisal": log_surprisal,
            "nondecreasing_under_delay": True,
        })
        moment_rows.append(log_moments)

    moment_matrix = np.asarray(moment_rows)
    surprisal_array = np.asarray(surprisals)
    for composition in all_compositions:
        occupation = sum(composition)
        log_survivals = np.zeros(block_bits + 1, dtype=np.float64)
        for index, count in enumerate(composition):
            if count:
                log_survivals += count * survival[index]
        first_survival = log_survivals[:-1]
        second_survival = log_survivals[1:]
        log_first = np.full(block_bits, -math.inf)
        terminal = np.isfinite(first_survival) & ~np.isfinite(second_survival)
        log_first[terminal] = first_survival[terminal]
        interior = np.isfinite(first_survival) & np.isfinite(second_survival)
        log_first[interior] = first_survival[interior] + np.log1p(
            -np.exp(second_survival[interior] - first_survival[interior])
        )
        log_average_moments = logsumexp(
            moment_matrix + log_first[None, :], axis=1
        )
        counts = (outer_rows - occupation, *composition)
        outer = log_multinomial(counts) + sum(
            count * math.log(active_bands[index][2])
            for index, count in enumerate(composition)
        )
        values = outer + np.minimum(
            0.0, log_average_moments + distance_cutoff * surprisal_array
        )
        winner = int(np.argmin(values))
        best[composition] = float(values[winner])
        witnesses[composition] = args.log_surprisals[winner]

    occupation_rows = []
    occupation_values = []
    for occupation in range(args.minimum_occupation, outer_rows + 1):
        same_q = [
            composition for composition in all_compositions
            if sum(composition) == occupation
        ]
        values = [best[composition] for composition in same_q]
        aggregate = float(logsumexp(values))
        occupation_values.append(aggregate)
        worst = max(same_q, key=lambda composition: best[composition])
        occupation_rows.append({
            "occupation": occupation,
            "composition_count": len(same_q),
            "log2_upper_diagnostic": aggregate / LOG2,
            "margin_bits_diagnostic": -aggregate / LOG2,
            "largest_pointwise_composition": list(worst),
            "largest_pointwise_log2": best[worst] / LOG2,
            "largest_pointwise_log_surprisal": witnesses[worst],
        })
    total = float(logsumexp(occupation_values))
    payload = {
        "schema": "fixed-outer-earliest-pivot-randomstepconv-high-q-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "constituent": constituent.name,
            "message_exponent": args.message_exponent,
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "memory_bits": args.memory_bits,
            "distance_cutoff": distance_cutoff,
            "minimum_occupation": args.minimum_occupation,
        },
        "bands": band_rows,
        "claim": {
            "covered_occupations": [args.minimum_occupation, outer_rows],
            "aggregate_log2_upper_diagnostic": total / LOG2,
            "aggregate_margin_bits_diagnostic": -total / LOG2,
            "closes_40_bits_diagnostic": total < -40 * LOG2,
            "occupation_rows": occupation_rows,
        },
        "checked_transfer_property": monotonicity_checks,
        "proof_reduction": [
            "Delete each word to its band's lower weight.",
            "After row permutation, each retained support is a uniform fixed-weight subset of coordinate regions.",
            "Keep only the earliest retained input across all active rows.",
            "Move that input to the final position of its region; the one-pulse tilted moment increases under this delay.",
            "Average the resulting one-pulse moment with the exact hypergeometric law of the earliest occupied region.",
        ],
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "The finite tilt list is not asserted optimal.",
            "A written outward proof must establish one-pulse delay monotonicity exactly.",
            "The inner is RandomStepConv, not RM2Sub.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
