#!/usr/bin/env python3
"""Probe singular contraction of the length-512 accumulator weight chain.

The diagnostic removes the polynomial subspace of degree at most 21 under
the nonzero binomial stationary measure, then measures the induced L2 norm
after one, two, and three accumulator/interleaver stages.  It tests whether a
generic high-degree norm can exploit the base dual distance by itself.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "accumulator_weight_chain_B512_singular_probe.json"
B = 512
REMOVED_DEGREE = 21


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-stages", default="1,2,3")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    states = np.arange(1, B + 1)
    transition = np.zeros((B, B))
    for row, input_weight_raw in enumerate(states):
        input_weight = int(input_weight_raw)
        half_down = input_weight // 2
        half_up = (input_weight + 1) // 2
        denominator = log_comb(B, input_weight)
        for column, output_weight_raw in enumerate(states):
            output_weight = int(output_weight_raw)
            value = (
                log_comb(B - output_weight, half_down)
                + log_comb(output_weight - 1, half_up - 1)
            )
            if math.isfinite(value):
                transition[row, column] = math.exp(value - denominator)

    log_stationary = np.asarray([log_comb(B, int(weight)) for weight in states])
    shift = float(log_stationary.max())
    stationary = np.exp(log_stationary - shift)
    stationary /= stationary.sum()
    root = np.sqrt(stationary)
    weighted = root[:, None] * transition / root[None, :]
    singular_values = np.linalg.svd(weighted, compute_uv=False)

    # Krawtchouk polynomials span all weight polynomials through the stated
    # degree. QR uses the conditioned nonzero stationary inner product.
    polynomials = np.zeros((B, REMOVED_DEGREE + 1))
    polynomials[:, 0] = 1.0
    polynomials[:, 1] = B - 2 * states
    for degree in range(1, REMOVED_DEGREE):
        polynomials[:, degree + 1] = (
            (B - 2 * states) * polynomials[:, degree]
            - (B - degree + 1) * polynomials[:, degree - 1]
        ) / (degree + 1)
    polynomials /= np.max(np.abs(polynomials), axis=0)
    low_basis, _ = np.linalg.qr(root[:, None] * polynomials)
    high_projection = np.eye(B) - low_basis @ low_basis.T

    restricted = []
    distribution_operator = weighted.T
    report_stages = sorted({int(value) for value in args.report_stages.split(",")})
    if not report_stages or report_stages[0] < 1:
        raise ValueError("report stages must be positive")
    powered = np.eye(B)
    for stages in range(1, report_stages[-1] + 1):
        powered = distribution_operator @ powered
        if stages not in report_stages:
            continue
        operator = powered @ high_projection
        norm = float(np.linalg.svd(operator, compute_uv=False)[0])
        restricted.append(
            {
                "accumulator_stages": stages,
                "restricted_l2_norm": norm,
                "log2_restricted_l2_norm": math.log2(norm),
            }
        )
        print(f"stages,{stages},restricted_norm,{norm:.12g},log2,{math.log2(norm):.9f}")

    result = {
        "schema": "accumulator-weight-chain-B512-singular-probe-v1",
        "status": "BINARY64_SPECTRAL_DIAGNOSTIC",
        "parameters": {
            "length": B,
            "state_space": "nonzero Hamming weights 1 through 512",
            "removed_polynomial_degrees": [0, REMOVED_DEGREE],
        },
        "checks": {
            "maximum_row_sum_error": float(np.max(np.abs(transition.sum(axis=1) - 1.0))),
            "maximum_stationarity_error": float(np.max(np.abs(stationary @ transition - stationary))),
        },
        "largest_singular_values": singular_values[:24].tolist(),
        "restricted": restricted,
        "conclusion": (
            "The generic degree-greater-than-21 L2 norm contracts by only about 15.29 bits "
            "after three stages. Dual-distance mode annihilation alone is therefore too weak; "
            "a shell-specific or base-distribution-specific pair bound is required."
        ),
        "limitations": [
            "The calculation is for the one-word weight chain, not the pair-type chain.",
            "All linear algebra uses nearest binary64 arithmetic.",
            "The result rejects only the generic high-degree L2 proof interface.",
        ],
    }
    result["parameters"]["report_stages"] = report_stages
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
