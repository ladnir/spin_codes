#!/usr/bin/env python3
"""Bound higher occupations with a positive banded spectrum mixture.

Partition the exact constituent spectrum into bands.  For each word in a
band, delete active coordinates until the band lower endpoint remains.
RandomStepConv monotonicity makes this deletion conservative.  An independent
uniform row permutation maps the retained support to a uniform subset.

Represent each uniform endpoint subset as a Bernoulli row conditioned on its
weight.  The resulting positive mixture dominates the complete constituent
counting measure.  A multiset dynamic program averages heterogeneous row
types exactly inside each transposed region.

Nearest binary64 log arithmetic makes the result diagnostic.
"""

from __future__ import annotations

import argparse
import functools
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
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "banded_spectrum_high_q_diagnostic.json"


def log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first, *rest)


def log_multinomial(parts: tuple[int, ...]) -> float:
    total = sum(parts)
    return math.lgamma(total + 1) - sum(math.lgamma(part + 1) for part in parts)


def multiset_product_averager(matrices: tuple[np.ndarray, ...]):
    """Return a cached evaluator for average products of matrix multisets."""
    @functools.lru_cache(maxsize=None)
    def recurse(state: tuple[int, ...]) -> np.ndarray:
        total = sum(state)
        if total == 0:
            return inner.log_identity()
        result = np.full((2, 2), -math.inf)
        for index, count in enumerate(state):
            if count == 0:
                continue
            following = list(state)
            following[index] -= 1
            term = inner.log_matmul(matrices[index], recurse(tuple(following)))
            term += math.log(count / total)
            result = np.logaddexp(result, term)
        return result

    return recurse


def parse_bands(values: list[int], block_bits: int) -> list[tuple[int, int]]:
    starts = sorted(set(values))
    if not starts or starts[0] <= 0 or starts[-1] > block_bits:
        raise ValueError("band starts must lie in [1,block_bits]")
    return [
        (start, starts[index + 1] - 1 if index + 1 < len(starts) else block_bits)
        for index, start in enumerate(starts)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constituent", choices=sorted(CONSTITUENTS), default="rm49")
    parser.add_argument("--message-exponent", type=int, default=11)
    parser.add_argument("--memory-bits", type=int, default=18)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--band-starts", type=int, nargs="+", default=[32, 48, 64, 96, 128, 192, 256])
    parser.add_argument(
        "--reference-probability",
        type=float,
        help="use one Bernoulli probability for every band instead of lower_weight/block_bits",
    )
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 0.7, 0.8, 1.0],
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
    if bands[0][0] > constituent.minimum_distance:
        parser.error("bands do not cover the constituent minimum distance")

    band_rows = []
    active_bands = []
    for lower, upper in bands:
        mass = sum(
            count for weight, count in spectrum.items()
            if weight and lower <= weight <= upper
        )
        if mass == 0:
            continue
        probability = (
            lower / block_bits
            if args.reference_probability is None
            else args.reference_probability
        )
        if not 0 < probability < 1:
            parser.error("reference-probability must lie in (0,1)")
        log_point_probability = (
            log_choose(block_bits, lower)
            + lower * math.log(probability)
            + (block_bits - lower) * math.log1p(-probability)
        )
        log_majorant = math.log(mass) - log_point_probability
        active_bands.append((lower, upper, mass, probability, log_majorant))
        band_rows.append({
            "lower_weight": lower,
            "upper_weight": upper,
            "exact_nonzero_mass": str(mass),
            "exact_nonzero_mass_log2": math.log2(mass),
            "bernoulli_probability": probability,
            "conditioning_probability_log2": log_point_probability / LOG2,
            "counting_majorant_log2": log_majorant / LOG2,
        })
    covered_mass = sum(item[2] for item in active_bands)
    if covered_mass != (1 << constituent.dimension) - 1:
        raise AssertionError("active bands do not cover every nonzero codeword")

    all_compositions = [
        composition
        for occupation in range(args.minimum_occupation, outer_rows + 1)
        for composition in compositions(occupation, len(active_bands))
    ]
    best = {composition: math.inf for composition in all_compositions}
    witnesses = {composition: math.nan for composition in all_compositions}

    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = inner.step_matrices(z, args.memory_bits)
        log_matrices = [inner.log_entries(zero)]
        for _lower, _upper, _mass, probability, _majorant in active_bands:
            candidate = (1.0 - probability) * zero + probability * active
            log_matrices.append(inner.log_entries(candidate))
        matrices = tuple(log_matrices)
        average_product = multiset_product_averager(matrices)
        for composition in all_compositions:
            occupation = sum(composition)
            counts = (outer_rows - occupation, *composition)
            region = average_product(counts)
            complete = inner.log_power(region, block_bits)
            moment = float(np.logaddexp(complete[0, 0], complete[0, 1]))
            log_no_candidate_input = block_bits * sum(
                count * math.log1p(-active_bands[index][3])
                for index, count in enumerate(composition)
            )
            if log_no_candidate_input >= moment:
                if log_no_candidate_input - moment > 1e-10:
                    raise ArithmeticError("zero-input atom exceeds total moment")
                log_live_moment = -math.inf
            else:
                log_live_moment = moment + math.log1p(
                    -math.exp(log_no_candidate_input - moment)
                )
            outer = log_multinomial(counts) + sum(
                count * active_bands[index][4]
                for index, count in enumerate(composition)
            )
            conditional = float(
                np.logaddexp(
                    log_no_candidate_input,
                    log_live_moment + distance_cutoff * surprisal,
                )
            )
            value = outer + min(0.0, conditional)
            if value < best[composition]:
                best[composition] = value
                witnesses[composition] = log_surprisal

    occupation_rows = []
    occupation_values = []
    for occupation in range(args.minimum_occupation, outer_rows + 1):
        selected = [
            value for composition, value in best.items()
            if sum(composition) == occupation
        ]
        aggregate = float(logsumexp(selected))
        occupation_values.append(aggregate)
        worst = max(
            (composition for composition in best if sum(composition) == occupation),
            key=lambda composition: best[composition],
        )
        occupation_rows.append({
            "occupation": occupation,
            "composition_count": len(selected),
            "log2_upper_diagnostic": aggregate / LOG2,
            "margin_bits_diagnostic": -aggregate / LOG2,
            "largest_pointwise_composition": list(worst),
            "largest_pointwise_log2": best[worst] / LOG2,
            "largest_pointwise_log_surprisal": witnesses[worst],
        })
    total = float(logsumexp(occupation_values))
    payload = {
        "schema": "fixed-outer-banded-spectrum-randomstepconv-high-q-v1",
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
            "band_count": len(active_bands),
            "common_reference_probability": args.reference_probability,
        },
        "bands": band_rows,
        "claim": {
            "covered_occupations": [args.minimum_occupation, outer_rows],
            "aggregate_log2_upper_diagnostic": total / LOG2,
            "aggregate_margin_bits_diagnostic": -total / LOG2,
            "closes_40_bits_diagnostic": total < -40 * LOG2,
            "occupation_rows": occupation_rows,
        },
        "proof_reduction": [
            "Partition the exact nonzero constituent counting measure by the displayed weight bands.",
            "Within each band, delete active coordinates to the lower endpoint; RandomStepConv monotonicity makes this conservative.",
            "After the uniform row permutation, the retained support is a uniform endpoint-weight subset.",
            "Dominate that support by the displayed Bernoulli law times its exact conditioning reciprocal.",
            "Expand the positive mixture across active rows and average each heterogeneous row composition exactly inside every region.",
            "Split the all-zero candidate-input atom exactly and apply the Chernoff factor only to the positive-input remainder.",
        ],
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "The finite tilt list is not asserted optimal.",
            "The inner is RandomStepConv, not RM2Sub.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
