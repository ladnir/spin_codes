#!/usr/bin/env python3
"""Expected spectrum of a regional expander followed by accumulators.

The sparse map sends 256 inputs to 512 intermediate coordinates using one
independent edge per input and per region.  Uniform interleavers precede the
accumulator stages.  All calculations are binary64 diagnostics.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp


LN2 = math.log(2.0)
DEFAULT_EXPANDER_ROOT = Path(
    r"C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes"
)


def log2_comb(n: int, k: np.ndarray | int) -> np.ndarray | float:
    return (gammaln(n + 1) - gammaln(np.asarray(k) + 1)
            - gammaln(n - np.asarray(k) + 1)) / LN2


def apply_accumulator(log_spectrum: np.ndarray) -> np.ndarray:
    block_size = len(log_spectrum) - 1
    result = np.full(block_size + 1, -math.inf)
    result[0] = log_spectrum[0]
    for h in range(1, block_size + 1):
        if not math.isfinite(float(log_spectrum[h])):
            continue
        down = h // 2
        up = (h + 1) // 2
        lo = up
        hi = block_size - down
        w = np.arange(lo, hi + 1)
        log_probability = (
            log2_comb(block_size - w, down)
            + log2_comb(w - 1, up - 1)
            - log2_comb(block_size, h)
        )
        result[lo : hi + 1] = np.logaddexp2(
            result[lo : hi + 1], log_spectrum[h] + log_probability
        )
    return result


def weighted_log2(log_spectrum: np.ndarray, marker: float) -> float:
    weights = np.arange(1, len(log_spectrum))
    terms = log_spectrum[1:] + weights * math.log2(marker)
    return float(logsumexp(terms * LN2) / LN2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expander-root", type=Path, default=DEFAULT_EXPANDER_ROOT)
    parser.add_argument("--k", type=int, default=256)
    parser.add_argument("--accumulators", type=int, default=4)
    parser.add_argument(
        "--region-lengths", type=int, nargs="+", default=(37,) * 8 + (36,) * 6
    )
    parser.add_argument("--markers", type=float, nargs="+", default=(0.05, 0.1, 0.2, 0.5))
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(args.expander_root / "scripts"))
    from regular_ec_diagnostic import regional_shell_distribution  # pylint: disable=import-error,import-outside-toplevel

    n = sum(args.region_lengths)
    spectrum = np.full(n + 1, -math.inf)
    kernel_logs: list[float] = []
    for support in range(args.k + 1):
        distribution = np.array([1.0])
        for length in args.region_lengths:
            shell = regional_shell_distribution(length, support)
            local = np.zeros(length + 1)
            for weight, probability in shell.items():
                local[weight] = probability
            distribution = np.convolve(distribution, local)
        nz = np.flatnonzero(distribution > 0.0)
        values = float(log2_comb(args.k, support)) + np.log2(distribution[nz])
        spectrum[nz] = np.logaddexp2(spectrum[nz], values)
        if support and distribution[0] > 0.0:
            kernel_logs.append(
                float(log2_comb(args.k, support)) + math.log2(distribution[0])
            )

    kernel_log2 = float(logsumexp(np.asarray(kernel_logs) * LN2) / LN2)

    stages: list[dict[str, object]] = []
    for stage in range(args.accumulators + 1):
        total_log2 = float(logsumexp(spectrum * LN2) / LN2)
        dense_random = {
            str(marker): args.k + n * (math.log2(1.0 + marker) - 1.0)
            for marker in args.markers
        }
        stages.append(
            {
                "accumulator_stages": stage,
                "expected_nonzero_kernel_words_log2": kernel_log2,
                "expected_spectrum_log2_by_weight": [
                    None if not math.isfinite(float(value)) else float(value)
                    for value in spectrum
                ],
                "weighted_positive_spectrum_log2": {
                    str(marker): weighted_log2(spectrum, marker)
                    for marker in args.markers
                },
                "random_linear_weighted_spectrum_log2": dense_random,
                "first_weight_expected_at_least_one": next(
                    (w for w in range(1, n + 1) if spectrum[w] >= 0.0), None
                ),
                "total_expected_codewords_log2": total_log2,
            }
        )
        if stage < args.accumulators:
            spectrum = apply_accumulator(spectrum)

    payload = {
        "schema": "block-expand-accumulate-expected-spectrum-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.k,
            "output_bits": n,
            "left_degree": len(args.region_lengths),
            "region_lengths_in_order": args.region_lengths,
        },
        "construction": "regional left-regular sparse map, then independent uniform interleaver and prefix accumulator at every stage",
        "stages": stages,
        "limitations": [
            "The calculation is binary64 and is not an outward certificate.",
            "It gives the ensemble expected spectrum, not concentration for one reused constituent.",
        ],
    }
    encoded = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
        print(
            json.dumps(
                {
                    "status": payload["status"],
                    "output": str(args.output),
                    "stages": len(stages),
                },
                indent=2,
            )
        )
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
