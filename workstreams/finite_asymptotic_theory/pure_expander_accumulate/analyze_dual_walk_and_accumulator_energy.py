#!/usr/bin/env python3
"""Analyze the dual random walk and accumulator Fourier energy.

The calculation is exact in integers except for the reported logarithms.  It
supports the proposed dual-kernel variance route; it is not itself a variance
or spectrum certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "dual_walk_accumulator_energy_K256_B512_r33.json"


def log2_integer(value: int) -> float:
    if value <= 0:
        return -math.inf
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def log2_sum(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return -math.inf
    maximum = max(finite)
    return maximum + math.log2(sum(2.0 ** (value - maximum) for value in finite))


def krawtchouk_row(length: int, shell_weight: int) -> list[int]:
    values = [math.comb(length, shell_weight)]
    values.append((length - 2 * shell_weight) * values[0] // length)
    for input_weight in range(1, length):
        numerator = (
            (length - 2 * shell_weight) * values[-1]
            - input_weight * values[-2]
        )
        denominator = length - input_weight
        if numerator % denominator:
            raise AssertionError("Krawtchouk recurrence lost integrality")
        values.append(numerator // denominator)
    return values


def accumulator_joint_counts(length: int) -> list[dict[int, int]]:
    """Return counts by wt(t) and wt(A^{-T}t)."""
    result: list[dict[int, int]] = []
    for weight in range(length + 1):
        counts: dict[int, int] = {}
        if weight == 0:
            counts[0] = 1
        else:
            zeros = length - weight
            maximum_runs = min(weight, zeros + 1)
            for runs in range(1, maximum_runs + 1):
                one_compositions = math.comb(weight - 1, runs - 1)
                # The string starts in one.  Internal zero runs are positive;
                # the final zero run may be empty.
                if runs - 1 <= zeros:
                    derivative_weight = 2 * runs - 1
                    counts[derivative_weight] = counts.get(derivative_weight, 0) + (
                        one_compositions * math.comb(zeros, runs - 1)
                    )
                # The string starts in zero.  The initial and internal zero
                # runs are positive; the final zero run may be empty.
                if runs <= zeros:
                    derivative_weight = 2 * runs
                    counts[derivative_weight] = counts.get(derivative_weight, 0) + (
                        one_compositions * math.comb(zeros, runs)
                    )
        if sum(counts.values()) != math.comb(length, weight):
            raise AssertionError("accumulator joint count has the wrong marginal")
        result.append(counts)
    return result


def walk_mixing(message_bits: int, right_degree: int) -> list[dict[str, object]]:
    denominator = math.comb(message_bits, right_degree)
    log_terms = []
    for weight in range(1, message_bits):
        numerator = abs(krawtchouk(message_bits, right_degree, weight))
        log_bias = (
            -math.inf
            if numerator == 0
            else log2_integer(numerator) - log2_integer(denominator)
        )
        log_terms.append((log2_integer(math.comb(message_bits, weight)), log_bias))
    rows = []
    for steps in range(1, 2 * message_bits + 1):
        spectral_mass = log2_sum(
            [multiplicity + steps * log_bias for multiplicity, log_bias in log_terms]
        )
        rows.append(
            {
                "steps": steps,
                "nontrivial_spectral_mass_log2": spectral_mass,
            }
        )
    return rows


def energy_rows(
    output_bits: int,
    shell_weights: list[int],
    joint_counts: list[dict[int, int]],
    cutoffs: list[int],
) -> list[dict[str, object]]:
    rows = []
    for shell_weight in shell_weights:
        kraw = krawtchouk_row(output_bits, shell_weight)
        energy_by_subset_weight = []
        for subset_weight, counts in enumerate(joint_counts):
            energy = sum(
                count * kraw[derivative_weight] ** 2
                for derivative_weight, count in counts.items()
            )
            energy_by_subset_weight.append(energy)
        total = sum(energy_by_subset_weight)
        expected_total = (1 << output_bits) * math.comb(output_bits, shell_weight)
        if total != expected_total:
            raise AssertionError("accumulator Fourier energy violates Parseval")
        cumulative = []
        for cutoff in cutoffs:
            low = sum(energy_by_subset_weight[1 : cutoff + 1])
            cumulative.append(
                {
                    "subset_weight_upper": cutoff,
                    "subset_weight_lower": 1,
                    "energy_fraction_log2": log2_integer(low) - log2_integer(total),
                }
            )
        rows.append(
            {
                "shell_weight": shell_weight,
                "total_energy_log2": log2_integer(total),
                "cumulative_low_subset_weight_energy": cumulative,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=512)
    parser.add_argument("--right-degree", type=int, default=33)
    parser.add_argument(
        "--shell-weights",
        type=int,
        nargs="+",
        default=[42, 80, 128, 192, 256, 320, 384, 432, 470],
    )
    parser.add_argument(
        "--cutoffs", type=int, nargs="+", default=[2, 4, 8, 16, 24, 32, 48, 64]
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits != 2 * args.message_bits:
        parser.error("the present report expects output_bits = 2 * message_bits")
    if not 1 <= args.right_degree <= args.message_bits:
        parser.error("invalid right degree")
    if any(not 0 <= weight <= args.output_bits for weight in args.shell_weights):
        parser.error("invalid shell weight")
    if any(not 0 <= cutoff <= args.output_bits for cutoff in args.cutoffs):
        parser.error("invalid cutoff")

    joint_counts = accumulator_joint_counts(args.output_bits)
    mixing = walk_mixing(args.message_bits, args.right_degree)
    energies = energy_rows(
        args.output_bits, args.shell_weights, joint_counts, args.cutoffs
    )
    payload = {
        "schema": "pure-ea-dual-walk-accumulator-energy-v1",
        "status": "EXACT_INTEGER_IDENTITIES_WITH_BINARY64_LOGS",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "shell_weights": args.shell_weights,
            "cutoffs": args.cutoffs,
        },
        "walk_mixing": mixing,
        "accumulator_fourier_energy": energies,
        "scope": [
            "The Krawtchouk values, joint counts, and Parseval checks are exact integers.",
            "Reported base-two logarithms are nearest binary64 diagnostics.",
            "These data do not bound the full covariance quadratic form.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for threshold in (0, -10, -20, -40, -80):
        first = next(
            (
                row["steps"]
                for row in mixing
                if row["nontrivial_spectral_mass_log2"] <= threshold
            ),
            None,
        )
        print(f"mixing_threshold_log2={threshold},first_steps={first}")
    for row in energies:
        low = row["cumulative_low_subset_weight_energy"][-1]
        print(
            f"shell={row['shell_weight']},cutoff={low['subset_weight_upper']},"
            f"energy_log2_fraction={low['energy_fraction_log2']:.6f}"
        )


if __name__ == "__main__":
    main()
