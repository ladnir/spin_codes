#!/usr/bin/env python3
"""Compare a modeled/expected local spectrum with a uniform random injection.

This is a binary64 diagnostic.  It does not authenticate the input spectrum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_RECEIPT = Path(
    "constructions/"
    "riffle_parityshear12_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s19/receipts/"
    "parityfanout31x33_outer256_d38_expected_spectrum.json"
)


def log2_binomial(n: int, k: int) -> float:
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def log2_nonzero_count(bits: int) -> float:
    return bits + math.log2(1.0 - math.ldexp(1.0, -bits))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", nargs="?", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument(
        "--minimum-weight",
        type=int,
        default=None,
        help="also report the maximum over weights at least this value",
    )
    parser.add_argument(
        "--maximum-weight",
        type=int,
        default=None,
        help="limit the reported restricted maximum to weights at most this value",
    )
    args = parser.parse_args()

    data = json.loads(args.receipt.read_text(encoding="utf-8"))
    parameters = data["parameters"]
    block_length = int(parameters["outer_bits"])
    dimension = int(parameters["outer_dimension"])

    log_numerator = log2_nonzero_count(dimension)
    log_denominator = log2_nonzero_count(block_length)
    comparisons: list[dict[str, float | int]] = []

    for cell in data["spectrum"]:
        weight = int(cell["weight"])
        if weight == 0:
            continue
        structured = float(cell["log2_expected_multiplicity"])
        random_reference = (
            log_numerator
            + log2_binomial(block_length, weight)
            - log_denominator
        )
        comparisons.append(
            {
                "weight": weight,
                "structured_log2": structured,
                "random_log2": random_reference,
                "excess_bits": structured - random_reference,
            }
        )

    maximum = max(comparisons, key=lambda cell: float(cell["excess_bits"]))
    output: dict[str, object] = {
        "status": "binary64 diagnostic; input spectrum not authenticated",
        "receipt": str(args.receipt),
        "block_length": block_length,
        "dimension": dimension,
        "maximum_excess_bits": maximum["excess_bits"],
        "maximum_excess_weight": maximum["weight"],
        "maximum_excess_per_output_bit": (
            float(maximum["excess_bits"]) / block_length
        ),
    }

    if args.minimum_weight is not None or args.maximum_weight is not None:
        restricted = [
            cell
            for cell in comparisons
            if (
                args.minimum_weight is None
                or int(cell["weight"]) >= args.minimum_weight
            )
            and (
                args.maximum_weight is None
                or int(cell["weight"]) <= args.maximum_weight
            )
        ]
        if not restricted:
            raise ValueError("no listed spectrum weight satisfies the requested range")
        restricted_maximum = max(
            restricted, key=lambda cell: float(cell["excess_bits"])
        )
        output["restricted_minimum_weight"] = args.minimum_weight
        output["restricted_maximum_weight"] = args.maximum_weight
        output["restricted_maximum_excess_bits"] = restricted_maximum[
            "excess_bits"
        ]
        output["restricted_maximum_excess_weight"] = restricted_maximum["weight"]

    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
