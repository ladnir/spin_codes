#!/usr/bin/env python3
"""Arb verifier for the split-low repair of the rejected dense type."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import flint
from flint import arb, ctx

from certify_ebch32_parityfanout_dense_cells_outward import (
    B,
    DISTANCE,
    L,
    DenseOutwardVerifier,
    definitely_negative,
    definitely_positive,
    exact_float,
)
from diagnose_ebch32_parityfanout_split_low_interval_cover import (
    BANDS,
    CELL_CONTRIBUTION_MARGIN,
    EXTERNAL_COUNTS,
)


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_split_low_rejected_point_cover_d11.json"
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_split_low_rejected_point_cover_outward.json"
DIAGNOSTIC_GENERATOR = WORKSTREAM / "diagnose_ebch32_parityfanout_split_low_interval_cover.py"
POINT_OPTIMIZER = WORKSTREAM / "diagnose_ebch32_parityfanout_split_low_point.py"
SETUP_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_ba_setup.py"
SETUP_RECEIPT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
SPECTRUM_SOURCE = WORKSTREAM / "ebch32_16_delta8_spectrum.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_interval(verifier: DenseOutwardVerifier, row):
    probabilities = [
        exact_float(float(value))
        for value in row["reference_type_probabilities"]
    ]
    values = [float(value) for value in row["reference_value_probabilities"]]
    if len(probabilities) != 5 or len(values) != 4:
        raise ArithmeticError("split-low witness dimension mismatch")
    for probability in probabilities:
        if not definitely_positive(probability):
            raise ArithmeticError("reference type probability is not positive")
    probability_sum = sum(probabilities, arb(0))
    probabilities = [probability / probability_sum for probability in probabilities]
    order_float = float(row["holder_order"])
    order = exact_float(order_float)
    if not definitely_positive(order - 1):
        raise ArithmeticError("Renyi order is not above one")
    surprisal_float = float(row["surprisal"])
    surprisal = exact_float(surprisal_float)
    band_moments = [
        verifier.band_log_moment(band, value, order_float)
        for band, value in zip(BANDS, values)
    ]
    bit_probability = sum(
        (
            probabilities[index + 1] * exact_float(values[index])
            for index in range(4)
        ),
        arb(0),
    )
    inner = verifier.inner_log_moment(bit_probability, surprisal_float)
    conjugate = (order - 1) / order
    convex_coefficient = arb(1) - B * conjugate
    if not definitely_negative(convex_coefficient):
        raise ArithmeticError("convexity coefficient sign is unresolved")

    endpoint_rows = []
    vertex_values = []
    for first_low_count in sorted(
        {
            int(row["lower_first_low_count"]),
            int(row["upper_first_low_count"]),
        }
    ):
        active_counts = [
            arb(first_low_count),
            arb(EXTERNAL_COUNTS[0] - first_low_count),
            arb(EXTERNAL_COUNTS[1]),
            arb(EXTERNAL_COUNTS[2]),
        ]
        occupation = sum(active_counts, arb(0))
        counts = [arb(L) - occupation, *active_counts]
        log_type_count = verifier.log_multinomial(counts)
        log_conditioning = log_type_count + sum(
            (
                count * probability.log()
                for count, probability in zip(counts, probabilities)
            ),
            arb(0),
        )
        raw = inner + DISTANCE * surprisal - B * log_conditioning
        if not definitely_negative(raw):
            raise ArithmeticError("raw reference exponent is not proved negative")
        value = (
            log_type_count
            + sum(
                (
                    active_counts[index] * band_moments[index]
                    for index in range(4)
                ),
                arb(0),
            )
            / order
            + conjugate * raw
        )
        vertex_values.append(value)
        endpoint_rows.append(
            {
                "first_low_count": first_low_count,
                "raw_upper": str(raw.upper()),
                "exponent_upper": str(value.upper()),
            }
        )
    maximum = vertex_values[0]
    for value in vertex_values[1:]:
        if bool(value.upper() > maximum.upper()):
            maximum = value
    integer_count = int(row["integer_count"])
    log_contribution = arb(integer_count).log() + maximum
    threshold = log_contribution + CELL_CONTRIBUTION_MARGIN * verifier.log2
    if not definitely_negative(threshold):
        raise ArithmeticError("interval contribution does not pass the configured threshold")
    return {
        "lower_first_low_count": int(row["lower_first_low_count"]),
        "upper_first_low_count": int(row["upper_first_low_count"]),
        "integer_count": integer_count,
        "log_contribution_upper": str(log_contribution.upper()),
        "log2_contribution_upper_display": float(log_contribution.upper()) / math.log(2.0),
        "cell_threshold_check": True,
        "endpoint_rows": endpoint_rows,
    }, log_contribution


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    ctx.prec = args.precision_bits
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source.get("status") != "COMPLETE_BINARY64_DIAGNOSTIC":
        raise ValueError("split-low diagnostic cover is incomplete")
    intervals = sorted(
        source["intervals"], key=lambda row: int(row["lower_first_low_count"])
    )
    expected = 0
    for row in intervals:
        if int(row["lower_first_low_count"]) != expected:
            raise ArithmeticError("split-low intervals have a gap or overlap")
        expected = int(row["upper_first_low_count"]) + 1
    if expected != EXTERNAL_COUNTS[0] + 1:
        raise ArithmeticError("split-low intervals do not cover the full range")

    verifier = DenseOutwardVerifier(args.precision_bits)
    rows = []
    aggregate = arb(0)
    for index, row in enumerate(intervals, 1):
        result, log_contribution = verify_interval(verifier, row)
        aggregate += log_contribution.exp()
        rows.append(result)
        print(
            f"interval={index}/{len(intervals)},range="
            f"[{result['lower_first_low_count']},{result['upper_first_low_count']}],"
            f"outward_log2={result['log2_contribution_upper_display']:.9f}",
            flush=True,
        )
    aggregate_log2 = aggregate.log() / verifier.log2
    claim = {
        "full_internal_split_range_verified": True,
        "verified_intervals": len(rows),
        "covered_first_low_counts": [0, EXTERNAL_COUNTS[0]],
        "all_intervals_pass_configured_threshold": all(
            bool(row["cell_threshold_check"]) for row in rows
        ),
        "aggregate_log2_sum_upper": str(aggregate_log2.upper()),
        "aggregate_margin_bits_lower_display": -float(aggregate_log2.upper()),
    }
    payload = {
        "schema": "ebch32-parityfanout31x33-b256-split-low-rejected-point-arb-v1",
        "status": "COMPLETE_OUTWARD_SPLIT_LOW_POINT_CERTIFICATE",
        "claim": claim,
        "arithmetic": {
            "library": "python-flint Arb",
            "library_version": flint.__version__,
            "precision_bits": args.precision_bits,
            "witness_interpretation": "JSON binary64 values are converted to exact integer ratios",
        },
        "external_three_band_counts": list(EXTERNAL_COUNTS),
        "split_bands": [list(band) for band in BANDS],
        "rows": rows,
        "dependencies": [
            {"path": INPUT.name, "sha256": sha256(INPUT)},
            {"path": Path(__file__).name, "sha256": sha256(Path(__file__).resolve())},
            {"path": DIAGNOSTIC_GENERATOR.name, "sha256": sha256(DIAGNOSTIC_GENERATOR)},
            {"path": POINT_OPTIMIZER.name, "sha256": sha256(POINT_OPTIMIZER)},
            {"path": SETUP_VERIFIER.name, "sha256": sha256(SETUP_VERIFIER)},
            {"path": SETUP_RECEIPT.name, "sha256": sha256(SETUP_RECEIPT)},
            {"path": SPECTRUM_SOURCE.name, "sha256": sha256(SPECTRUM_SOURCE)},
        ],
        "scope": (
            "This receipt certifies only the external three-band type "
            "(3344,1,21), refined by the split of weights 24..100."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(claim, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
