#!/usr/bin/env python3
"""Arb verifier for the 10% one-band occupation-interval cover."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

from flint import arb

for key, value in (
    ("SPIN_EBCH_DISTANCE", "209715"),
    ("SPIN_EBCH_LOWER_WEIGHT", "24"),
    ("SPIN_EBCH_PARITY_FANOUT", "1"),
):
    if key in os.environ and os.environ[key] != value:
        raise RuntimeError(f"{key} must equal {value}")
    os.environ[key] = value

from certify_ebch32_parityfanout_dense_cells_outward import (
    DenseOutwardVerifier,
    definitely_negative,
    definitely_positive,
    exact_float,
)


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10.json"
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10_outward.json"
B = 256
L = 8192
N = B * L
DISTANCE = 209_715
MIN_OCCUPATION = 65
BAND = (24, 232)
CELL_MARGIN_BITS = 44


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    parameters = source["parameters"]
    if int(parameters["distance"]) != DISTANCE:
        raise ArithmeticError("diagnostic distance mismatch")
    if list(parameters["occupation_interval"]) != [MIN_OCCUPATION, L]:
        raise ArithmeticError("diagnostic occupation range mismatch")
    if list(parameters["outer_weight_band"]) != list(BAND):
        raise ArithmeticError("diagnostic band mismatch")

    verifier = DenseOutwardVerifier(256)
    if verifier.envelope is None:
        raise ArithmeticError("RM2Sub envelope did not load")
    log2 = verifier.log2
    rows = []
    aggregate = arb(0)
    cursor = MIN_OCCUPATION

    for index, row in enumerate(source["intervals"], 1):
        lower, upper = map(int, row["occupation_interval"])
        if lower != cursor or upper < lower or upper > L:
            raise ArithmeticError("interval cover has a gap, overlap, or invalid endpoint")
        cursor = upper + 1

        probabilities = [
            exact_float(float(value))
            for value in row["reference_type_probabilities"]
        ]
        values = [float(value) for value in row["reference_value_probabilities"]]
        if len(probabilities) != 2 or len(values) != 1:
            raise ArithmeticError("one-band witness dimension mismatch")
        probability_sum = sum(probabilities, arb(0))
        if not definitely_positive(probability_sum):
            raise ArithmeticError("reference probability sum is not positive")
        probabilities = [probability / probability_sum for probability in probabilities]
        for probability in probabilities:
            if not definitely_positive(probability):
                raise ArithmeticError("reference type probability is not positive")

        order_float = float(row["holder_order"])
        order = exact_float(order_float)
        if not definitely_positive(order - 1):
            raise ArithmeticError("Renyi order is not above one")
        surprisal_float = float(row["surprisal"])
        surprisal = exact_float(surprisal_float)
        if not definitely_positive(surprisal):
            raise ArithmeticError("surprisal is not positive")
        value_probability = exact_float(values[0])
        band_moment = verifier.band_log_moment(BAND, values[0], order_float)
        bit_probability = probabilities[1] * value_probability
        inner = verifier.inner_log_moment(bit_probability, surprisal_float)
        conjugate = (order - 1) / order
        convex_coefficient = arb(1) - B * conjugate
        if not definitely_negative(convex_coefficient):
            raise ArithmeticError("convexity coefficient is not negative")

        endpoint_values = []
        raw_uppers = []
        for occupation in (lower, upper):
            counts = [arb(L - occupation), arb(occupation)]
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
                raise ArithmeticError(
                    f"interval {lower}..{upper} raw exponent is not negative"
                )
            value = (
                log_type_count
                + occupation * band_moment / order
                + conjugate * raw
            )
            endpoint_values.append(value)
            raw_uppers.append(str(raw.upper()))

        maximum = endpoint_values[0]
        if bool(endpoint_values[1].upper() > maximum.upper()):
            maximum = endpoint_values[1]
        log_contribution = arb(upper - lower + 1).log() + maximum
        if not definitely_negative(log_contribution + CELL_MARGIN_BITS * log2):
            raise ArithmeticError(
                f"interval {lower}..{upper} does not pass its cell threshold"
            )
        aggregate += log_contribution.exp()
        rows.append(
            {
                "occupation_interval": [lower, upper],
                "raw_endpoint_upper_bounds": raw_uppers,
                "log2_contribution_upper": str((log_contribution / log2).upper()),
                "margin_bits_lower_display": -float((log_contribution / log2).upper()),
            }
        )
        print(
            f"interval={index}/{len(source['intervals'])},q={lower}..{upper},"
            f"margin={rows[-1]['margin_bits_lower_display']:.9f}",
            flush=True,
        )

    if cursor != L + 1:
        raise ArithmeticError("interval cover does not reach the final occupation")
    if not definitely_positive(aggregate):
        raise ArithmeticError("aggregate interval mass is not positive")
    aggregate_log2 = aggregate.log() / log2
    if not definitely_negative(aggregate_log2 + 40):
        raise ArithmeticError("dense aggregate does not pass 40 bits")
    if not definitely_negative(aggregate_log2 + 468):
        raise ArithmeticError("dense aggregate does not pass 468 integer bits")

    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-b256-one-band-q-cover-d10-arb-v1",
        "status": "COMPLETE_OUTWARD_DENSE_CERTIFICATE",
        "claim": {
            "bad_weight_at_most": DISTANCE,
            "occupations": [MIN_OCCUPATION, L],
            "intervals_verified": len(rows),
            "gap_free_integer_cover": True,
            "fixed_witness_convexity_checked": True,
            "aggregate_log2_upper": str(aggregate_log2.upper()),
            "aggregate_margin_bits_lower_display": -float(aggregate_log2.upper()),
            "certified_integer_margin_bits": 468,
            "comparison_to_2^-40": True,
        },
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": 256,
            "witness_interpretation": "JSON binary64 values converted to exact rational numbers",
        },
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "outer_weight_band": list(BAND),
        },
        "rows": rows,
        "dependencies": [
            {"path": INPUT.name, "sha256": sha256(INPUT)},
            {"path": Path(__file__).name, "sha256": sha256(Path(__file__).resolve())},
            {
                "path": "certify_ebch32_parityfanout_dense_cells_outward.py",
                "sha256": sha256(
                    WORKSTREAM / "certify_ebch32_parityfanout_dense_cells_outward.py"
                ),
            },
            {
                "path": "ebch32_parityfanout31x33_ba3_B256_setup_outward.json",
                "sha256": sha256(
                    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
                ),
            },
        ],
        "scope": (
            "The receipt bounds the expected number of bad nonzero messages "
            "with 65 through 8192 active outer rows. It assumes the exact "
            "conditional row ensemble and the audited RM2Sub transfer interface."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
