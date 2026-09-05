#!/usr/bin/env python3
"""Build deterministic shell caps from primal tails and dual distance.

The setup samples one BA constituent and repeats it in every outer row.  A
first-moment event excludes low and high primal weights and low dual weights.
On that event, the dual distance makes the weight of a uniform codeword match
the Binomial(B,1/2) moments through degree d_perp-1.  The Christoffel bound
then caps each remaining shell using the optimal square of a polynomial of
degree floor((d_perp-1)/2).
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from evaluate_single_random_constituent_highprob_bands import band_majorant


WORKSTREAM = Path(__file__).resolve().parent
SPECTRA = WORKSTREAM / "ebch128x4_ba0_16_B512_expected_spectra.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba12_christoffel_caps_probe.json"
LENGTH = 512
DIMENSION = 256
LOG2 = math.log(2.0)


def spectrum_logs(stage: dict[str, object]) -> list[float]:
    result = [-math.inf] * (LENGTH + 1)
    for row in stage["spectrum"]:  # type: ignore[index]
        value = row["log2_expected_multiplicity"]
        if value is not None:
            result[int(row["weight"])] = float(value)
    return result


def logadd2(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(2.0 ** (left - maximum) + 2.0 ** (right - maximum))


def excluded_low(log_spectrum: list[float], failure_bits: float) -> tuple[int, float]:
    cumulative = -math.inf
    accepted = -math.inf
    first_allowed = 1
    for weight in range(1, LENGTH + 1):
        trial = logadd2(cumulative, log_spectrum[weight])
        if trial > -failure_bits:
            break
        cumulative = trial
        accepted = cumulative
        first_allowed = weight + 1
    return first_allowed, accepted


def excluded_high(log_spectrum: list[float], failure_bits: float) -> tuple[int, float]:
    cumulative = -math.inf
    accepted = -math.inf
    last_allowed = LENGTH
    for weight in range(LENGTH, 0, -1):
        trial = logadd2(cumulative, log_spectrum[weight])
        if trial > -failure_bits:
            break
        cumulative = trial
        accepted = cumulative
        last_allowed = weight - 1
    return last_allowed, accepted


def krawtchouk_prefix(length: int, weight: int, maximum_degree: int) -> list[int]:
    values = [1]
    if maximum_degree == 0:
        return values
    values.append(length - 2 * weight)
    for degree in range(1, maximum_degree):
        numerator = (
            (length - 2 * weight) * values[degree]
            - (length - degree + 1) * values[degree - 1]
        )
        following, remainder = divmod(numerator, degree + 1)
        if remainder:
            raise AssertionError("Krawtchouk recurrence lost integrality")
        values.append(following)
    return values


def christoffel_cap(weight: int, degree: int) -> tuple[int, Fraction]:
    values = krawtchouk_prefix(LENGTH, weight, degree)
    kernel = sum(
        (Fraction(value * value, math.comb(LENGTH, index)) for index, value in enumerate(values)),
        Fraction(0),
    )
    upper = Fraction(1 << DIMENSION, 1) / kernel
    return min(1 << DIMENSION, upper.numerator // upper.denominator), kernel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=int, default=12)
    parser.add_argument("--per-event-failure-bits", type=float, default=44.0)
    parser.add_argument("--low-upper", type=int, default=79)
    parser.add_argument("--central-upper", type=int, default=432)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    payload = json.loads(SPECTRA.read_text(encoding="utf-8"))
    primal_stage = payload["stages"][args.stages]
    dual_stage = payload["dual_stages"][args.stages]
    if int(primal_stage["accumulators"]) != args.stages:
        raise AssertionError("primal stage index mismatch")
    if int(dual_stage["accumulators"]) != args.stages:
        raise AssertionError("dual stage index mismatch")
    primal = spectrum_logs(primal_stage)
    dual = spectrum_logs(dual_stage)
    first_allowed, low_failure = excluded_low(primal, args.per_event_failure_bits)
    last_allowed, high_failure = excluded_high(primal, args.per_event_failure_bits)
    dual_distance, dual_failure = excluded_low(dual, args.per_event_failure_bits)
    degree = (dual_distance - 1) // 2

    caps = np.zeros(LENGTH + 1, dtype=object)
    kernels = []
    for weight in range(first_allowed, last_allowed + 1):
        cap, kernel = christoffel_cap(weight, degree)
        caps[weight] = cap
        kernels.append(
            {
                "weight": weight,
                "cap": str(cap),
                "cap_log2": math.log2(cap) if cap else -math.inf,
                "christoffel_kernel_numerator": str(kernel.numerator),
                "christoffel_kernel_denominator": str(kernel.denominator),
            }
        )

    support = [weight for weight in range(1, LENGTH + 1) if int(caps[weight])]
    low = band_majorant(caps, min(support), args.low_upper, block_bits=LENGTH)
    central = band_majorant(
        caps, args.low_upper + 1, args.central_upper, block_bits=LENGTH
    )
    high = band_majorant(caps, args.central_upper + 1, max(support), block_bits=LENGTH)
    event_logs = np.asarray([low_failure, high_failure, dual_failure]) * LOG2
    event_union = float(logsumexp(event_logs) / LOG2)
    result = {
        "schema": "ebch128x4-ba-christoffel-caps-probe-v1",
        "status": "EXACT_CAPS_FROM_BINARY64_EXPECTED_TAILS",
        "construction": {
            "base": "direct sum of four fixed extended BCH [128,64,22] codes",
            "accumulator_stages": args.stages,
            "sampler": "one independent uniform 512-coordinate permutation before each accumulator",
            "reuse": "the sampled constituent is repeated in every SPIN outer row",
        },
        "event": {
            "per_tail_target_bits": args.per_event_failure_bits,
            "first_allowed_primal_weight": first_allowed,
            "last_allowed_primal_weight": last_allowed,
            "dual_distance": dual_distance,
            "primal_low_failure_log2_upper": low_failure,
            "primal_high_failure_log2_upper": high_failure,
            "dual_low_failure_log2_upper": dual_failure,
            "union_failure_log2_upper": event_union,
            "union_margin_bits": -event_union,
        },
        "christoffel": {
            "degree": degree,
            "claim": "A_w/2^K <= 1 / sum_{j=0}^m K_j(w)^2/binom(B,j)",
            "reason": "dual distance makes all polynomials through degree 2m have the Binomial(B,1/2) expectation; optimize the evaluation inequality over degree-m polynomials",
        },
        "bands": {"low": low, "central": central, "high": high},
        "caps": kernels,
        "limitations": [
            "The Christoffel caps are exact rational consequences of the stated primal and dual distances.",
            "The expected primal and dual tail sums are inherited from a nearest-binary64 spectrum recurrence.",
            "The event union and Bernoulli band optimization use nearest binary64 arithmetic.",
            "No structured routing or RandomStepConv transfer is included.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event": result["event"], "bands": result["bands"]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
