#!/usr/bin/env python3
"""Binary64 column-Hölder diagnostic for repeated EBCH128 at Q=L.

The fixed [128,64] constituent is linear and has full coordinate support.
Consequently, a uniformly selected nonzero constituent word has a one in
each fixed coordinate with probability 2^63/(2^64-1).  Independent
coordinate permutations preserve that marginal.  The 128 transposed-region
moments may be dependent; Hölder's inequality bounds their product using
only the common one-region marginal.  This avoids a pointwise spectrum
envelope at the all-active endpoint.

This is a nearest-binary64 witness search.  It is not an outward certificate.
"""

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
    impulse_matrices_upper,
    load_activation_and_live,
)


OUTPUT = WORKSTREAM / "ebch128_rm2sub_s19_k20_d11_qL_column_holder_diagnostic.json"
OUTER_BITS = 128
OUTER_DIMENSION = 64
MESSAGE_BITS = 1 << 20
OUTER_BLOCKS = MESSAGE_BITS // OUTER_DIMENSION
OUTPUT_BITS = OUTER_BITS * OUTER_BLOCKS
DISTANCE = math.floor(0.11 * OUTPUT_BITS)
STEP_BITS = 128
EPOCHS_PER_REGION = OUTER_BLOCKS // STEP_BITS
OUTER_NONZERO_WORDS = (1 << OUTER_DIMENSION) - 1
BIT_PROBABILITY = (1 << (OUTER_DIMENSION - 1)) / OUTER_NONZERO_WORDS


def log_choose(total: int, selected: np.ndarray) -> np.ndarray:
    return (
        gammaln(total + 1.0)
        - gammaln(selected + 1.0)
        - gammaln(total - selected + 1.0)
    )


def region_shell_matrices(
    impulses: list[tuple[float, float, float, float]],
) -> np.ndarray:
    """Return the averaged two-state region matrix for every input weight."""
    impulse_array = np.asarray(impulses, dtype=np.float64).reshape(
        STEP_BITS + 1, 2, 2
    )
    current = impulse_array.copy()
    completed_bits = STEP_BITS
    for _ in range(1, EPOCHS_PER_REGION):
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
            products = current @ impulse_array[added]
            updated[destinations] += (
                products * np.exp(log_weight)[:, None, None]
            )
        current = updated
        completed_bits = next_bits
    return current


def evaluate(z: float) -> dict[str, object]:
    activation, live = load_activation_and_live()
    regions = region_shell_matrices(
        impulse_matrices_upper(z, activation, live)
    )
    weights = np.arange(OUTER_BLOCKS + 1, dtype=np.float64)
    log_binomial = (
        log_choose(OUTER_BLOCKS, weights)
        + weights * math.log(BIT_PROBABILITY)
        + (OUTER_BLOCKS - weights) * math.log1p(-BIT_PROBABILITY)
    )

    def objective(log_ratio: float) -> float:
        ratio = math.exp(log_ratio)
        row_zero = regions[:, 0, 0] + regions[:, 0, 1] * ratio
        row_live = regions[:, 1, 0] / ratio + regions[:, 1, 1]
        norms = np.maximum(row_zero, row_live)
        norms = np.maximum(norms, np.nextafter(0.0, math.inf))
        holder_log = float(
            logsumexp(log_binomial + OUTER_BITS * np.log(norms))
        )
        terminal_factor = max(1.0, 1.0 / ratio)
        return (
            OUTER_BLOCKS * math.log(OUTER_NONZERO_WORDS)
            + math.log(terminal_factor)
            + holder_log
            - DISTANCE * math.log(z)
        )

    result = minimize_scalar(
        objective,
        bounds=(-20.0, 20.0),
        method="bounded",
        options={"xatol": 1e-11},
    )
    value = float(result.fun)
    return {
        "z_binary64_hex_exact": z.hex(),
        "positive_norm_ratio_v1_over_v0": math.exp(float(result.x)),
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(result.success),
        "region_matrix_entries_nonnegative": bool(np.min(regions) >= 0.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--z",
        type=float,
        nargs="+",
        default=(0.10, 0.11, 0.12, 0.13, 0.14, 0.15),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows = []
    for z in args.z:
        row = evaluate(z)
        rows.append(row)
        print(json.dumps(row), flush=True)
    best = min(rows, key=lambda row: float(row["log2_upper"]))
    payload = {
        "schema": "ebch128-rm2sub-s19-k20-qL-column-holder-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": MESSAGE_BITS,
            "output_bits": OUTPUT_BITS,
            "distance": DISTANCE,
            "outer_bits": OUTER_BITS,
            "outer_dimension": OUTER_DIMENSION,
            "outer_blocks": OUTER_BLOCKS,
            "occupation": OUTER_BLOCKS,
            "outer_nonzero_words_per_active_block": OUTER_NONZERO_WORDS,
            "one_coordinate_probability": BIT_PROBABILITY,
        },
        "method": {
            "region_transfer": "exact finite uniform-shell two-state envelope",
            "region_product": "common positive weighted norm",
            "column_dependence": "Holder with exponent 128",
            "outer_property": "full-support linear [128,64] code",
            "arithmetic": "nearest binary64",
        },
        "best": best,
        "rows": rows,
        "limitations": [
            "The calculation is not outward rounded.",
            "The receipt covers only occupation Q=L.",
            "The two-state epoch transfer is an entrywise upper envelope.",
        ],
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
