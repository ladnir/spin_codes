#!/usr/bin/env python3
"""Hypothetical BA shell caps under an explicit variance-inflation lemma.

For each shell, assume

    Var(A_w) <= F E[A_w].

For every shell, the program takes the smaller valid cap supplied by Markov's
inequality or Cantelli's inequality.  The resulting caps and two-band
Bernoulli majorants quantify the exact concentration theorem that would let
the one-shot random-constituent proof transfer to EBCH128x4--BA.

The variance inequality is an input hypothesis, not a proved fact.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from evaluate_single_random_constituent_highprob_bands import band_majorant


WORKSTREAM = Path(__file__).resolve().parent
SPECTRA = WORKSTREAM / "ebch128x4_ba0_16_B512_expected_spectra.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
B = 512
K = 256


def load_means(stage: int) -> list[mp.mpf]:
    payload = json.loads(SPECTRA.read_text(encoding="utf-8"))
    source = next(row for row in payload["stages"] if int(row["accumulators"]) == stage)
    result = [mp.mpf(0)] * (B + 1)
    for row in source["spectrum"]:
        value = row["log2_expected_multiplicity"]
        if value is not None:
            result[int(row["weight"])] = mp.power(2, mp.mpf(value))
    result[0] = mp.mpf(0)
    return result


def caps_from_variance(
    means: list[mp.mpf],
    delta: mp.mpf,
    inflation: mp.mpf,
    forced_zero_through: int,
) -> tuple[np.ndarray, mp.mpf, dict[str, int]]:
    caps = np.zeros(B + 1, dtype=object)
    failure = mp.mpf(0)
    methods = {"forced_zero": 0, "markov": 0, "cantelli": 0, "support": 0}
    for weight in range(1, B + 1):
        mean = means[weight]
        if weight <= forced_zero_through:
            caps[weight] = 0
            failure += mean
            methods["forced_zero"] += 1
            continue
        maximum = min((1 << K) - 1, math.comb(B, weight))

        markov_cap = min(maximum, max(0, int(mp.ceil(mean / delta)) - 1))
        if markov_cap == maximum:
            markov_failure = mp.mpf(0)
        else:
            markov_failure = mean / (markov_cap + 1)

        if mean == 0:
            cantelli_cap = 0
            cantelli_failure = mp.mpf(0)
        else:
            variance = inflation * mean
            deviation = mp.sqrt(variance * (1 - delta) / delta)
            cantelli_cap = min(maximum, max(0, int(mp.ceil(mean + deviation)) - 1))
            if cantelli_cap == maximum:
                cantelli_failure = mp.mpf(0)
            else:
                actual = mp.mpf(cantelli_cap + 1) - mean
                if actual <= 0:
                    raise ArithmeticError("Cantelli cap does not exceed the mean")
                cantelli_failure = variance / (variance + actual * actual)

        choices = (
            (markov_cap, markov_failure, "markov"),
            (cantelli_cap, cantelli_failure, "cantelli"),
        )
        cap, shell_failure, method = min(
            choices, key=lambda row: (row[0], row[1])
        )
        if cap == maximum and shell_failure == 0:
            method = "support"
        caps[weight] = cap
        failure += shell_failure
        methods[method] += 1
    return caps, failure, methods


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument(
        "--variance-inflation-bits",
        type=float,
        nargs="+",
        default=[0.0, 2.0, 4.0, 6.0, 9.0, 12.0],
    )
    parser.add_argument("--per-shell-failure-bits", type=float, default=51.0)
    parser.add_argument("--low-upper", type=int, default=79)
    parser.add_argument("--central-upper", type=int, default=432)
    parser.add_argument("--forced-zero-through", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    rows = []
    with mp.workdps(200):
        delta = mp.power(2, -args.per_shell_failure_bits)
        for stage in args.stages:
            means = load_means(stage)
            for inflation_bits in args.variance_inflation_bits:
                inflation = mp.power(2, inflation_bits)
                caps, failure, cap_methods = caps_from_variance(
                    means, delta, inflation, args.forced_zero_through
                )
                support = [weight for weight in range(1, B + 1) if int(caps[weight])]
                low = band_majorant(caps, min(support), args.low_upper, block_bits=B)
                central = band_majorant(
                    caps, args.low_upper + 1, args.central_upper, block_bits=B
                )
                high = band_majorant(
                    caps, args.central_upper + 1, max(support), block_bits=B
                )
                defect_probability = min(
                    float(low["value_probability"]),
                    1.0 - float(high["value_probability"]),
                )
                defect_log2 = 1.0 + max(
                    float(low["log_majorant"]), float(high["log_majorant"])
                ) / math.log(2.0)
                row = {
                    "accumulator_stages": stage,
                    "variance_inflation_bits": inflation_bits,
                    "variance_hypothesis": f"Var(A_w) <= 2^{inflation_bits} E[A_w] for every nonzero shell",
                    "event_failure_log2_upper": float(mp.log(failure, 2)),
                    "cap_method_shell_counts": cap_methods,
                    "cap_support": [min(support), max(support)],
                    "cap_mass_log2": math.log2(sum(int(value) for value in caps)),
                    "caps": [int(value) for value in caps],
                    "low": low,
                    "central": central,
                    "high": high,
                    "merged_defect": {
                        "value_probability": defect_probability,
                        "log2_majorant": defect_log2,
                    },
                }
                rows.append(row)
                print(
                    f"stage,{stage},variance_bits,{inflation_bits:.3f},"
                    f"failure,{row['event_failure_log2_upper']:.6f},"
                    f"support,{row['cap_support']},"
                    f"defect_majorant,{defect_log2:.6f},"
                    f"central_majorant,{float(central['log_majorant']) / math.log(2.0):.6f}",
                    flush=True,
                )

    result = {
        "schema": "ebch128x4-ba-variance-cap-requirements-v1",
        "status": "CONDITIONAL_BINARY64_CAP_DIAGNOSTIC",
        "parameters": {
            "outer": "four EBCH [128,64,22] blocks followed by the stated number of accumulator stages",
            "length": B,
            "dimension": K,
            "per_shell_failure_bits": args.per_shell_failure_bits,
            "bands": [[None, args.low_upper], [args.low_upper + 1, args.central_upper], [args.central_upper + 1, None]],
            "forced_zero_through": args.forced_zero_through,
        },
        "rows": rows,
        "limitations": [
            "No variance-inflation hypothesis in this receipt is proved.",
            "Expected BA spectra and band optimization use nearest binary64 arithmetic.",
            "The caps have not yet been passed through the complete all-occupation transfer verifier.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
