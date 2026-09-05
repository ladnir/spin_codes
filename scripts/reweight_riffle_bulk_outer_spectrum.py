#!/usr/bin/env python3
"""Reweight a fixed-reference Riffle bulk receipt to a new outer spectrum.

The inner transfer and Chernoff optimization depend on the selected Bernoulli
reference, not on the outer spectrum.  When both receipts use the same fixed
reference probability, changing the spectrum only changes the per-active-block
likelihood envelope.  This script applies that exact additive correction and
recomputes the occupation and aggregate ledgers.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    load_outer_spectrum_logs,
    logsumexp,
    spectrum_bernoulli_envelope_log,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outer-spectrum", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    receipt = json.loads(args.input.read_text(encoding="utf-8"))
    parameters = receipt["parameters"]
    probabilities = {
        float(row["best_candidate_probability"])
        for row in receipt["occupation_rows"]
    }
    if len(probabilities) != 1:
        raise ValueError("reweighting currently requires one fixed reference probability")
    probability = probabilities.pop()
    old_envelopes = {
        float(row["probability"]): row
        for row in receipt["envelope"]["bernoulli_envelopes"]
    }
    if probability not in old_envelopes:
        raise ValueError("selected probability is absent from old envelope metadata")
    old_bits = float(old_envelopes[probability]["log2_mass"])

    outer_bits = int(parameters["outer_bits"])
    spectrum = load_outer_spectrum_logs(args.outer_spectrum, outer_bits)
    new_log, new_weight = spectrum_bernoulli_envelope_log(
        outer_bits,
        spectrum,
        probability,
        int(parameters["tail_endpoint_width"]),
    )
    new_bits = new_log / LOG2
    delta_bits = new_bits - old_bits

    for row in receipt["occupation_rows"]:
        occupation = int(row["active_regular_outer_blocks"])
        inner = float(row["inner_log2_upper"]) + occupation * delta_bits
        pointwise = float(row["outer_log2_envelope"]) + inner
        row["inner_log2_upper"] = inner
        row["pointwise_log2_upper"] = pointwise
        row["pointwise_margin_bits"] = -pointwise

    pointwise_logs = np.asarray(
        [float(row["pointwise_log2_upper"]) * LOG2 for row in receipt["occupation_rows"]]
    )
    partial_log2 = float(logsumexp(pointwise_logs) / LOG2)
    receipt["partial_log2_upper"] = partial_log2
    receipt["partial_margin_bits"] = -partial_log2
    receipt["dominant_rows"] = sorted(
        receipt["occupation_rows"],
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]
    receipt["parameters"]["outer_spectrum"] = str(args.outer_spectrum)
    receipt["envelope"]["bernoulli_envelopes"] = [
        {
            "probability": probability,
            "log2_mass": new_bits,
            "maximizing_weight": new_weight,
        }
    ]
    receipt["envelope"]["log2_eta"] = new_bits
    receipt["envelope"]["maximizing_regular_weight"] = new_weight
    receipt["reweighting"] = {
        "schema": "riffle-fixed-reference-outer-reweight-v1",
        "source_receipt": str(args.input),
        "old_envelope_bits_per_active_block": old_bits,
        "new_envelope_bits_per_active_block": new_bits,
        "delta_bits_per_active_block": delta_bits,
        "justification": (
            "The Bernoulli reference probability and every inner transfer are "
            "unchanged. The outer likelihood envelope is an additive constant "
            "per active block, independent of the Chernoff tilt."
        ),
    }
    receipt["scope"] += (
        " Reweighted exactly at the fixed Bernoulli reference to the outer "
        "spectrum named in parameters.outer_spectrum."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"old_envelope_bits,{old_bits:.12f}")
    print(f"new_envelope_bits,{new_bits:.12f}")
    print(f"delta_bits_per_active,{delta_bits:.12f}")
    print(f"partial_margin_bits,{-partial_log2:.9f}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
