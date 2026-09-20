#!/usr/bin/env python3
"""Outward verifier for the conditional BA-128 random-cap event.

The only hypothesis is the stated singular-value bound for one permuted
accumulator.  Given that hypothesis, this verifier reconstructs the certified
random-[512,256] cap vector and proves with Arb intervals that a fixed full-rank
[512,256] base code followed by 128 independent permute-accumulate stages
obeys those caps except with probability below 2^-40.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx

from certify_single_random_constituent_dense_outward import exact_caps


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ba128_pair_singular_random_cap_event_outward_B512.json"

BLOCK_BITS = 512
DIMENSION = 256
STAGES = 128
SINGULAR_NUMERATOR = 63
SINGULAR_DENOMINATOR = 1000
# Central shells subtract two roughly 2^510 terms to expose a variance-scale
# remainder.  Use enough precision that this cancellation is certified rather
# than hidden by interval radius.
PRECISION_BITS = 1024


def point(numerator: int, denominator: int = 1) -> arb:
    return arb(numerator) / denominator


def upper_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.upper()), math.inf)


def lower_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.lower()), -math.inf)


def main() -> None:
    ctx.prec = PRECISION_BITS
    caps, _ = exact_caps()
    cap_list = [int(cap) for cap in caps]

    ambient_nonzero_int = (1 << BLOCK_BITS) - 1
    base_nonzero_int = (1 << DIMENSION) - 1
    ambient_pairs_int = ambient_nonzero_int * (ambient_nonzero_int - 1)
    base_pairs_int = base_nonzero_int * (base_nonzero_int - 1)

    ambient_nonzero = arb(ambient_nonzero_int)
    base_nonzero = arb(base_nonzero_int)
    ambient_pairs = arb(ambient_pairs_int)
    base_pairs = arb(base_pairs_int)
    chi_one = ambient_nonzero / base_nonzero - 1
    chi_pair = ambient_pairs / base_pairs - 1
    contraction = point(SINGULAR_NUMERATOR, SINGULAR_DENOMINATOR) ** STAGES

    failure = arb(0)
    shell_rows = []
    worst_tail_log2_upper = None
    worst_weight = None
    log_two = arb(2).log()
    for weight in range(1, BLOCK_BITS + 1):
        shell_size_int = math.comb(BLOCK_BITS, weight)
        shell_size = arb(shell_size_int)
        beta = shell_size / ambient_nonzero
        pi_pair = shell_size * (shell_size - 1) / ambient_pairs
        epsilon_one = contraction * (chi_one * beta * (1 - beta)).sqrt()
        epsilon_pair = contraction * (chi_pair * pi_pair * (1 - pi_pair)).sqrt()

        marginal_lower = beta - epsilon_one
        marginal_upper = beta + epsilon_one
        pair_upper = pi_pair + epsilon_pair
        if not bool(marginal_lower.lower() > 0):
            raise ArithmeticError(f"unresolved positive marginal lower bound at weight {weight}")
        if not bool(marginal_upper.upper() < 1):
            raise ArithmeticError(f"unresolved marginal upper bound below one at weight {weight}")
        if not bool(pair_upper.upper() < 1):
            raise ArithmeticError(f"unresolved pair upper bound below one at weight {weight}")

        mean_lower = base_nonzero * marginal_lower
        mean_upper = base_nonzero * marginal_upper
        variance_upper = (
            base_pairs * pair_upper + mean_upper - mean_lower * mean_lower
        )
        if not bool(variance_upper.lower() > 0):
            raise ArithmeticError(f"unresolved positive variance bound at weight {weight}")

        cap = cap_list[weight]
        if cap == 0:
            tail = mean_upper
            method = "Markov-at-one"
        else:
            deviation = arb(cap + 1) - mean_upper
            if not bool(deviation.lower() > 0):
                raise ArithmeticError(f"cap does not exceed mean at weight {weight}")
            tail = variance_upper / (variance_upper + deviation * deviation)
            method = "Cantelli"
        if not bool(tail.lower() > 0) or not bool(tail.upper() < 1):
            raise ArithmeticError(f"unresolved tail probability at weight {weight}")
        failure += tail
        tail_log2 = tail.log() / log_two
        tail_log2_upper = tail_log2.upper()
        if worst_tail_log2_upper is None or bool(
            tail_log2_upper > worst_tail_log2_upper
        ):
            worst_tail_log2_upper = tail_log2_upper
            worst_weight = weight
        shell_rows.append(
            {
                "weight": weight,
                "cap": cap,
                "method": method,
                "mean_lower": str(mean_lower.lower()),
                "mean_upper": str(mean_upper.upper()),
                "variance_upper": str(variance_upper.upper()),
                "tail_log2_upper": str(tail_log2_upper),
            }
        )

    failure_log2 = failure.log() / log_two
    margin = -failure_log2
    if not bool(margin.lower() > 40):
        raise AssertionError("BA-128 cap-event margin does not exceed 40 bits")

    cap_bytes = ",".join(str(cap) for cap in cap_list).encode("ascii")
    payload = {
        "schema": "ba128-pair-singular-random-cap-event-outward-v1",
        "status": "CONDITIONAL_OUTWARD_CERTIFICATE",
        "claim": {
            "cap_event_failure_upper": upper_float(failure),
            "cap_event_failure_log2_upper": str(failure_log2.upper()),
            "cap_event_failure_log2_upper_float": upper_float(failure_log2),
            "cap_event_margin_bits_lower": str(margin.lower()),
            "cap_event_margin_bits_lower_float": lower_float(margin),
            "comparison_to_2^-40": True,
            "worst_single_shell_weight": worst_weight,
            "worst_single_shell_log2_upper": str(worst_tail_log2_upper),
        },
        "parameters": {
            "block_bits": BLOCK_BITS,
            "dimension": DIMENSION,
            "accumulator_stages": STAGES,
            "singular_value_upper_hypothesis": "63/1000",
            "arithmetic_precision_bits": PRECISION_BITS,
            "cap_source": "certify_single_random_constituent_dense_outward.exact_caps",
            "cap_vector_sha256": hashlib.sha256(cap_bytes).hexdigest(),
        },
        "probability_space": (
            "a fixed full-rank binary [512,256] base code; 128 independent "
            "uniform coordinate permutations, each followed by prefix accumulation"
        ),
        "hypothesis": (
            "in stationary L2, the second singular values of the one-word and "
            "ordered-distinct-pair type operators for one permuted accumulator "
            "are each at most 63/1000"
        ),
        "proof_accounting": {
            "marginal": "stationary-L2 contraction from the uniform distribution on the nonzero base codewords",
            "second_factorial_moment": "stationary-L2 contraction from the uniform distribution on ordered distinct base-codeword pairs",
            "tail": "Markov for zero caps and Cantelli for positive caps",
            "union": "sum over all 512 nonzero Hamming-weight shells",
        },
        "shells": shell_rows,
        "limitations": [
            "The singular-value hypothesis is not proved by this verifier.",
            "This proves the outer cap event only; the existing cap-conditional distance receipts are reused separately.",
            "The 128-stage count is a universal fallback and is not an optimized implementation proposal.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"claim": payload["claim"], "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()
