#!/usr/bin/env python3
"""Test an exact column-Hölder reduction for Q=L at distance 11%."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from certify_golay_ba_rm2sub_finite_q2_64 import (  # noqa: E402
    B,
    DISTANCE,
    L,
    impulse_matrices_upper,
    load_activation_and_live,
)


OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_qL_column_holder_d11.json"
STEP_BITS = 128
EPOCHS_PER_REGION = L // STEP_BITS
OUTER_NONZERO_WORDS = (1 << (B // 2)) - 1
BIT_PROBABILITY = (1 << (B // 2 - 1)) / OUTER_NONZERO_WORDS


def log_choose(total: int, selected: np.ndarray) -> np.ndarray:
    return (
        gammaln(total + 1.0)
        - gammaln(selected + 1.0)
        - gammaln(total - selected + 1.0)
    )


def region_shell_matrices(impulses: list[tuple[float, float, float, float]]) -> np.ndarray:
    """Return the averaged 2x2 region matrix for every input weight."""
    impulse_array = np.asarray(impulses, dtype=np.float64).reshape(STEP_BITS + 1, 2, 2)
    current = impulse_array.copy()
    completed_bits = STEP_BITS
    for epoch in range(1, EPOCHS_PER_REGION):
        next_bits = completed_bits + STEP_BITS
        updated = np.zeros((next_bits + 1, 2, 2), dtype=np.float64)
        source = np.arange(completed_bits + 1, dtype=np.float64)
        log_source_shell = log_choose(completed_bits, source)
        for added in range(STEP_BITS + 1):
            destinations = source.astype(np.int64) + added
            log_weight = (
                math.lgamma(STEP_BITS + 1.0)
                - math.lgamma(added + 1.0)
                - math.lgamma(STEP_BITS - added + 1.0)
                + log_source_shell
                - log_choose(next_bits, destinations.astype(np.float64))
            )
            weights = np.exp(log_weight)
            products = current @ impulse_array[added]
            updated[destinations] += products * weights[:, None, None]
        current = updated
        completed_bits = next_bits
    if current.shape[0] != L + 1:
        raise ArithmeticError("region recurrence length mismatch")
    return current


def evaluate(z: float) -> dict[str, object]:
    activation, live = load_activation_and_live()
    regions = region_shell_matrices(impulse_matrices_upper(z, activation, live))
    weights = np.arange(L + 1, dtype=np.float64)
    log_binomial = (
        log_choose(L, weights)
        + weights * math.log(BIT_PROBABILITY)
        + (L - weights) * math.log1p(-BIT_PROBABILITY)
    )

    def objective(log_ratio: float) -> float:
        ratio = math.exp(log_ratio)
        row_zero = regions[:, 0, 0] + regions[:, 0, 1] * ratio
        row_live = regions[:, 1, 0] / ratio + regions[:, 1, 1]
        norms = np.maximum(row_zero, row_live)
        # A zero here can be binary64 underflow.  Replacing it by the least
        # positive subnormal is conservative for this diagnostic upper path.
        norms = np.maximum(norms, np.nextafter(0.0, math.inf))
        log_expectation = float(logsumexp(log_binomial + B * np.log(norms)))
        terminal_factor = max(1.0, 1.0 / ratio)
        return (
            L * math.log(OUTER_NONZERO_WORDS)
            + math.log(terminal_factor)
            + log_expectation
            - DISTANCE * math.log(z)
        )

    result = minimize_scalar(
        objective,
        bounds=(-20.0, 20.0),
        method="bounded",
        options={"xatol": 1e-11},
    )
    value = float(result.fun)
    ratio = math.exp(float(result.x))
    return {
        "z_binary64_hex_exact": z.hex(),
        "positive_norm_ratio_v1_over_v0": ratio,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(result.success),
        "region_matrix_entries_nonnegative": bool(np.min(regions) >= 0.0),
        "region_matrix_maximum": float(np.max(regions)),
        "region_matrix_minimum_positive": float(np.min(regions[regions > 0.0])),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=float, nargs="+", default=(0.123595,))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows = []
    for z in args.z:
        if not 0.0 < z < 1.0:
            raise ValueError("z must lie in (0,1)")
        row = evaluate(z)
        rows.append(row)
        print(json.dumps(row, indent=2), flush=True)
    best = min(rows, key=lambda row: float(row["log2_upper"]))
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-qL-column-holder-d11-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "occupation": L,
            "distance": DISTANCE,
            "outer_nonzero_words_per_row": OUTER_NONZERO_WORDS,
            "one_coordinate_probability_for_uniform_nonzero_outer_word": BIT_PROBABILITY,
        },
        "method": {
            "region_transfer": "exact finite uniform-shell 2x2 envelope",
            "region_product": "common positive weighted norm",
            "column_dependence": "Holder with exponent B across B regions",
            "shell_mixture": "summed implicitly from exact one-coordinate balance",
        },
        "best": best,
        "rows": rows,
        "limitations": [
            "The calculation uses nearest binary64 and is not an outward certificate.",
            "The column-Hölder reduction requires a written proof audit.",
            "This receipt covers Q=L only.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
