#!/usr/bin/env python3
"""Test random-constituent spectrum caps for a conditionally mixed BA-t code.

The input cap vector is the exact vector used by the certified random-[512,256]
outer transfer.  This program assumes an upper bound on the second singular
value of both the one-word and ordered-distinct-pair accumulator channels.
It then applies stationary-L2 contraction and Cantelli's inequality to bound
the probability that BA-t violates any entry of that fixed cap vector.

The singular-value premise is an open proof obligation.  All arithmetic after
that premise is evaluated with high-precision mpmath arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp

from certify_single_random_constituent_dense_outward import exact_caps


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ba_pair_singular_random_cap_event_B512.json"


def log2(value: mp.mpf) -> mp.mpf:
    return mp.log(value, 2) if value > 0 else mp.ninf


def shell_tail(
    *,
    block_bits: int,
    dimension: int,
    stages: int,
    singular_upper: mp.mpf,
    weight: int,
    cap: int,
) -> dict[str, object]:
    """Return a conditional upper bound on Pr[A_weight > cap]."""

    ambient_nonzero = mp.mpf(2) ** block_bits - 1
    base_nonzero = mp.mpf(2) ** dimension - 1
    ambient_pairs = ambient_nonzero * (ambient_nonzero - 1)
    base_pairs = base_nonzero * (base_nonzero - 1)
    shell_size = mp.mpf(math.comb(block_bits, weight))

    beta = shell_size / ambient_nonzero
    pi_pair = shell_size * (shell_size - 1) / ambient_pairs
    chi_one = ambient_nonzero / base_nonzero - 1
    chi_pair = ambient_pairs / base_pairs - 1
    contraction = singular_upper**stages
    epsilon_one = contraction * mp.sqrt(chi_one * beta * (1 - beta))
    epsilon_pair = contraction * mp.sqrt(chi_pair * pi_pair * (1 - pi_pair))

    marginal_lower = max(mp.mpf(0), beta - epsilon_one)
    marginal_upper = min(mp.mpf(1), beta + epsilon_one)
    pair_upper = min(mp.mpf(1), pi_pair + epsilon_pair)
    mean_lower = base_nonzero * marginal_lower
    mean_upper = base_nonzero * marginal_upper

    # A(A-1) counts ordered pairs of distinct constituent codewords.  Hence
    # E[A(A-1)] <= base_pairs * pair_upper.  For every admissible mean,
    # Var(A) = E[A(A-1)] + E[A] - E[A]^2 is at most the expression below.
    variance_upper = max(
        mp.mpf(0), base_pairs * pair_upper + mean_upper - mean_lower**2
    )

    threshold = mp.mpf(cap + 1)
    if cap == 0:
        tail_upper = min(mp.mpf(1), mean_upper)
        method = "Markov-at-one"
    elif threshold <= mean_upper:
        tail_upper = mp.mpf(1)
        method = "trivial"
    else:
        deviation = threshold - mean_upper
        tail_upper = min(
            mp.mpf(1), variance_upper / (variance_upper + deviation**2)
        )
        method = "Cantelli"

    return {
        "weight": weight,
        "cap": cap,
        "method": method,
        "log2_mean_lower": float(log2(mean_lower)),
        "log2_mean_upper": float(log2(mean_upper)),
        "log2_variance_upper": float(log2(variance_upper)),
        "log2_tail_upper": float(log2(tail_upper)),
        "log2_one_word_error": float(log2(epsilon_one)),
        "log2_pair_error": float(log2(epsilon_pair)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--minimum-stages", type=int, default=120)
    parser.add_argument("--maximum-stages", type=int, default=160)
    parser.add_argument("--singular-upper", default="0.063")
    parser.add_argument("--target-failure-bits", type=float, default=40.0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if (args.block_bits, args.dimension) != (512, 256):
        parser.error("the imported certified cap vector is for [512,256]")
    if not 0 <= args.minimum_stages <= args.maximum_stages:
        parser.error("invalid stage interval")

    caps, random_event_failure = exact_caps()
    with mp.workdps(220):
        singular_upper = mp.mpf(args.singular_upper)
        rows = []
        first_passing = None
        for stages in range(args.minimum_stages, args.maximum_stages + 1):
            shells = [
                shell_tail(
                    block_bits=args.block_bits,
                    dimension=args.dimension,
                    stages=stages,
                    singular_upper=singular_upper,
                    weight=weight,
                    cap=caps[weight],
                )
                for weight in range(1, args.block_bits + 1)
            ]
            failure = mp.fsum(
                mp.power(2, mp.mpf(shell["log2_tail_upper"]))
                for shell in shells
                if shell["log2_tail_upper"] != float("-inf")
            )
            worst = max(shells, key=lambda shell: shell["log2_tail_upper"])
            row = {
                "accumulator_stages": stages,
                "event_failure_log2_upper": float(log2(failure)),
                "worst_shell": worst,
            }
            rows.append(row)
            if (
                first_passing is None
                and row["event_failure_log2_upper"]
                <= -args.target_failure_bits
            ):
                first_passing = row
            print(
                f"stages,{stages},event_failure_bits,"
                f"{-row['event_failure_log2_upper']:.12f},"
                f"worst_weight,{worst['weight']}",
                flush=True,
            )

        payload = {
            "schema": "ba-pair-singular-random-cap-event-v1",
            "status": "CONDITIONAL_HIGH_PRECISION_BOUND",
            "parameters": {
                "block_bits": args.block_bits,
                "dimension": args.dimension,
                "stage_interval": [args.minimum_stages, args.maximum_stages],
                "one_and_pair_second_singular_value_upper_hypothesis": args.singular_upper,
                "target_event_failure_bits": args.target_failure_bits,
                "cap_source": "certify_single_random_constituent_dense_outward.exact_caps",
            },
            "probability_space": (
                "a fixed full-rank binary [512,256] base code; t independent "
                "uniform coordinate permutations, each followed by prefix accumulation"
            ),
            "random_code_cap_event_failure_log2_upper": float(
                log2(mp.mpf(random_event_failure.numerator) / random_event_failure.denominator)
            ),
            "first_passing": first_passing,
            "rows": rows,
            "theorem_interface": {
                "hypothesis": (
                    "the one-word and ordered-distinct-pair type operators of one "
                    "permuted accumulator have second singular value at most lambda "
                    "in their stationary L2 spaces"
                ),
                "conclusion": (
                    "the BA-t constituent obeys the exact certified random-[512,256] "
                    "spectrum cap vector except with the stated union-bound probability"
                ),
            },
            "limitations": [
                "The singular-value hypothesis is an open proof obligation.",
                "The receipt uses high-precision floating-point arithmetic, not outward intervals.",
                "Passing this event gate reuses the conditional distance transfer for the same cap vector; it does not itself reverify that transfer.",
                "The stage count is a universal fallback and does not exploit the base-code spectrum.",
            ],
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"first_passing": first_passing}, indent=2))
        print(f"output={args.output}")


if __name__ == "__main__":
    main()
