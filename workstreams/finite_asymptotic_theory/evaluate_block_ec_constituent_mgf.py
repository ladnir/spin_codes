#!/usr/bin/env python3
"""Compare a small two-sided EC constituent with a random linear code.

This is a binary64 diagnostic.  It evaluates the exact regional first-moment
law of the rate-half two-sided regular EC ensemble at fixed output markers.
The default parent has dimensions [510,255], degrees 10/5, and memory 15.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


DEFAULT_EXPANDER_ROOT = Path(
    r"C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expander-root", type=Path, default=DEFAULT_EXPANDER_ROOT)
    parser.add_argument("--k", type=int, default=255)
    parser.add_argument("--left-degree", type=int, default=10)
    parser.add_argument("--right-degree", type=int, default=5)
    parser.add_argument("--memory", type=int, default=15)
    parser.add_argument(
        "--markers",
        type=float,
        nargs="+",
        default=(0.03, 0.05, 0.1, 0.2, 0.5),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scripts = args.expander_root / "scripts"
    sys.path.insert(0, str(scripts))
    from binary_biregular_diagnostic import (  # pylint: disable=import-error,import-outside-toplevel
        biregular_parity_shell_counts,
    )
    from expander_bounds import log_weight_mgf  # pylint: disable=import-error,import-outside-toplevel
    from regular_ec_diagnostic import (  # pylint: disable=import-error,import-outside-toplevel
        uniform_slice_transfer_matrices,
    )

    k = args.k
    d_left = args.left_degree
    d_right = args.right_degree
    if k % d_right:
        raise ValueError("right degree must divide k")
    region_length = k // d_right
    n = d_left * region_length
    rows: list[dict[str, float]] = []

    for marker in args.markers:
        if not 0.0 < marker < 1.0:
            raise ValueError("markers must lie in (0,1)")
        slices = uniform_slice_transfer_matrices(
            length=region_length,
            max_weight=region_length,
            output_marker=marker,
            memory=args.memory,
        )
        support_logs: list[float] = []
        for support in range(1, k + 1):
            denominator, counts = biregular_parity_shell_counts(
                left_vertices=k,
                right_degree=d_right,
                support_size=support,
            )
            region = np.zeros_like(slices[0])
            for weight, count in counts.items():
                region += (count / denominator) * slices[weight]
            log_moment = log_weight_mgf(
                region,
                d_left,
                args.memory,
            )
            support_logs.append(math.lgamma(k + 1) - math.lgamma(support + 1)
                                - math.lgamma(k - support + 1) + log_moment)

        ec_log2 = float(logsumexp(support_logs) / math.log(2.0))
        random_log2 = (
            math.log2(2.0**k - 1.0)
            + n * (math.log2(1.0 + marker) - 1.0)
        )
        worst_support = int(np.argmax(np.asarray(support_logs))) + 1
        rows.append(
            {
                "marker": marker,
                "ec_nonzero_mgf_log2": ec_log2,
                "random_linear_nonzero_mgf_log2": random_log2,
                "ec_minus_random_bits": ec_log2 - random_log2,
                "largest_ec_support": worst_support,
                "largest_ec_support_log2": support_logs[worst_support - 1]
                / math.log(2.0),
            }
        )

    payload = {
        "schema": "block-ec-constituent-mgf-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": k,
            "output_bits": n,
            "left_degree": d_left,
            "right_degree": d_right,
            "regions": d_left,
            "region_length": region_length,
            "memory": args.memory,
        },
        "quantity": "E[sum over nonzero messages of marker^output_weight]",
        "rows": rows,
        "limitations": [
            "The calculation is binary64 and is not an outward certificate.",
            "It is a first moment and does not prove concentration of one reused constituent.",
            "A non-power-of-two parameter set needs a wrapper or a separate power-of-two regional variant.",
        ],
    }
    encoded = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
