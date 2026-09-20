#!/usr/bin/env python3
"""Conditional BA-t variance sweep from a pair maximal-correlation bound.

Assume the rank-two pair-type channel of one permuted accumulator has second
singular value at most lambda.  The same bound then applies to its one-word
factor.  Chi-square contraction from an arbitrary binary [B,K] base code
gives a genus-two-free shell-variance bound after t independent stages.

The singular-value premise is not proved here.  All subsequent inequalities
are evaluated with high-precision arithmetic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ba_pair_singular_mixing_sweep_B512.json"


def log2(value: mp.mpf) -> mp.mpf:
    return mp.log(value, 2) if value > 0 else mp.ninf


def shell_bound(
    *,
    block_bits: int,
    dimension: int,
    stages: int,
    singular_upper: mp.mpf,
    weight: int,
) -> dict[str, object]:
    ambient_nonzero = mp.mpf(2) ** block_bits - 1
    codewords = mp.mpf(2) ** dimension
    base_nonzero = codewords - 1
    ambient_pairs = ambient_nonzero * (ambient_nonzero - 1)
    base_pairs = base_nonzero * (base_nonzero - 1)
    shell_size_int = mp.binomial(block_bits, weight)
    shell_size = mp.mpf(shell_size_int)

    beta = shell_size / ambient_nonzero
    pi_pair = shell_size * (shell_size - 1) / ambient_pairs
    chi_one = ambient_nonzero / base_nonzero - 1
    chi_pair = ambient_pairs / base_pairs - 1
    contraction = singular_upper**stages
    epsilon_one = contraction * mp.sqrt(chi_one * beta * (1 - beta))
    epsilon_pair = contraction * mp.sqrt(chi_pair * pi_pair * (1 - pi_pair))

    marginal_lower = max(mp.mpf(0), beta - epsilon_one)
    mean_lower = base_nonzero * marginal_lower
    stationary_dependence = abs(pi_pair - beta * beta)
    marginal_shift = epsilon_one * (2 * beta + epsilon_one)
    correlation_error = epsilon_pair + stationary_dependence + marginal_shift
    if mean_lower > 0:
        mixing_factor = 1 + base_pairs * correlation_error / mean_lower
    else:
        mixing_factor = mp.inf
    support_factor = shell_size
    factor = min(mixing_factor, support_factor)
    return {
        "weight": weight,
        "log2_ambient_shell_size": float(log2(shell_size)),
        "log2_stationary_marginal": float(log2(beta)),
        "log2_contraction": float(log2(contraction)),
        "log2_one_word_error": float(log2(epsilon_one)),
        "log2_pair_error": float(log2(epsilon_pair)),
        "mixing_variance_factor_log2_upper": (
            None if not mp.isfinite(mixing_factor) else float(log2(mixing_factor))
        ),
        "support_variance_factor_log2_upper": float(log2(support_factor)),
        "combined_variance_factor_log2_upper": float(log2(factor)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--minimum-stages", type=int, default=100)
    parser.add_argument("--maximum-stages", type=int, default=150)
    parser.add_argument("--singular-upper", default="0.063")
    parser.add_argument("--target-factor-bits", type=float, default=9.0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if not 0 < args.dimension < args.block_bits:
        parser.error("dimension must lie strictly between zero and block length")
    if not 0 <= args.minimum_stages <= args.maximum_stages:
        parser.error("invalid stage interval")

    with mp.workdps(220):
        singular_upper = mp.mpf(args.singular_upper)
        rows = []
        first_passing = None
        for stages in range(args.minimum_stages, args.maximum_stages + 1):
            shells = [
                shell_bound(
                    block_bits=args.block_bits,
                    dimension=args.dimension,
                    stages=stages,
                    singular_upper=singular_upper,
                    weight=weight,
                )
                for weight in range(1, args.block_bits + 1)
            ]
            worst = max(
                shells,
                key=lambda row: float(row["combined_variance_factor_log2_upper"]),
            )
            row = {
                "accumulator_stages": stages,
                "worst_weight": worst["weight"],
                "worst_variance_factor_log2_upper": worst[
                    "combined_variance_factor_log2_upper"
                ],
                "worst_shell": worst,
            }
            rows.append(row)
            if (
                first_passing is None
                and float(row["worst_variance_factor_log2_upper"])
                <= args.target_factor_bits
            ):
                first_passing = row
            print(
                f"stages,{stages},worst_weight,{row['worst_weight']},"
                f"variance_factor_bits,{row['worst_variance_factor_log2_upper']:.12f}",
                flush=True,
            )

        payload = {
            "schema": "ba-pair-singular-mixing-sweep-v1",
            "status": "CONDITIONAL_HIGH_PRECISION_BOUND",
            "parameters": {
                "block_bits": args.block_bits,
                "dimension": args.dimension,
                "stage_interval": [args.minimum_stages, args.maximum_stages],
                "pair_second_singular_value_upper_hypothesis": args.singular_upper,
                "target_variance_factor_bits": args.target_factor_bits,
            },
            "theorem_interface": {
                "hypothesis": "the rank-two pair-type operator of one permuted accumulator has second singular value at most lambda in stationary L2",
                "conclusion": "for every fixed binary [B,K] base code and every shell, Var(A_w after t stages) <= F_t E[A_w after t stages]",
                "probability_space": "t independent uniform coordinate permutations, each followed by prefix accumulation",
            },
            "first_passing": first_passing,
            "rows": rows,
            "limitations": [
                "The pair singular-value hypothesis is an open proof obligation.",
                "The calculation does not claim that the binary64 one-word singular value proves the pair bound.",
                "The shell bound is universal and conservative; it does not use the EBCH spectrum or dual distance.",
                "The result is a variance lemma only; the existing cap and SPIN transfer verifiers must still be rerun for the selected stage count.",
            ],
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"first_passing": first_passing}, indent=2))
        print(f"output={args.output}")


if __name__ == "__main__":
    main()
