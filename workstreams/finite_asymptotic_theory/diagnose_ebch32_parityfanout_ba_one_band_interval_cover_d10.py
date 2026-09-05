#!/usr/bin/env python3
"""One-band occupation-interval diagnostic at the 10% finite target."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

for key, value in (
    ("SPIN_EBCH_RELATIVE_DISTANCE", "0.10"),
    ("SPIN_EBCH_DISTANCE", "209715"),
    ("SPIN_EBCH_LOWER_WEIGHT", "24"),
    ("SPIN_EBCH_PARITY_FANOUT", "1"),
):
    if key in os.environ and os.environ[key] != value:
        raise RuntimeError(f"{key} must equal {value}")
    os.environ[key] = value

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_ebch32_ba_k20_three_band_q_probe import optimize  # noqa: E402
from diagnose_ebch32_ba_k20_three_band_qL import (  # noqa: E402
    B,
    EPOCHS,
    L,
    N,
    band_log_moment,
    conditioned_spectrum,
)
from diagnose_finite_k20_band_compositions_qL import log_multinomial  # noqa: E402
from diagnose_finite_k20_renyi_dense import high_precision_transfer_log_moment  # noqa: E402


DISTANCE = 209_715
MIN_OCCUPATION = 65
BAND = (24, 232)
CELL_MARGIN_BITS = 44.0
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10.json"


def evaluate_fixed(envelope, spectrum, occupation: int, witness: dict[str, object]):
    probabilities = np.asarray(witness["reference_type_probabilities"], dtype=np.float64)
    values = np.asarray(witness["reference_value_probabilities"], dtype=np.float64)
    surprisal = float(witness["surprisal"])
    order = float(witness["holder_order"])
    moment = band_log_moment(spectrum, BAND, float(values[0]), order)
    bit_probability = float(probabilities[1] * values[0])
    inner = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal, epochs=EPOCHS
    )
    counts = (L - occupation, occupation)
    log_type_count = log_multinomial(counts)
    log_conditioning = log_type_count + sum(
        count * math.log(float(probability))
        for count, probability in zip(counts, probabilities)
    )
    raw = inner + DISTANCE * surprisal - B * log_conditioning
    conjugate = (order - 1.0) / order
    value = log_type_count + occupation * moment / order + conjugate * min(0.0, raw)
    return raw, -value / math.log(2.0), moment / math.log(2.0)


def certify_interval(envelope, spectrum, lower: int, upper: int, inherited=()):
    for witness in inherited:
        endpoint_rows = [evaluate_fixed(envelope, spectrum, q, witness) for q in (lower, upper)]
        margin = min(row[1] for row in endpoint_rows)
        if all(row[0] < 0.0 for row in endpoint_rows) and (
            margin > CELL_MARGIN_BITS + math.log2(upper - lower + 1)
        ):
            selected = witness
            mode = "inherited-witness"
            break
    else:
        midpoint = (lower + upper) // 2
        selected = optimize(envelope, spectrum, (midpoint,), (BAND,))
        endpoint_rows = [evaluate_fixed(envelope, spectrum, q, selected) for q in (lower, upper)]
        margin = min(row[1] for row in endpoint_rows)
        if not all(row[0] < 0.0 for row in endpoint_rows) or (
            margin <= CELL_MARGIN_BITS + math.log2(upper - lower + 1)
        ):
            return None, (selected,)
        mode = "center-witness"
    return {
        "occupation_interval": [lower, upper],
        "optimizer_mode": mode,
        "reference_type_probabilities": selected["reference_type_probabilities"],
        "reference_value_probabilities": selected["reference_value_probabilities"],
        "surprisal": selected["surprisal"],
        "holder_order": selected["holder_order"],
        "endpoint_raw_log2_upper": [row[0] / math.log(2.0) for row in endpoint_rows],
        "endpoint_margin_bits": [row[1] for row in endpoint_rows],
        "band_log2_renyi_moment": endpoint_rows[0][2],
        "minimum_endpoint_margin_bits": margin,
        "log2_contribution_upper": math.log2(upper - lower + 1) - margin,
    }, ()


def main() -> None:
    if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") != "1":
        raise RuntimeError("set SPIN_EBCH_PARITY_FANOUT=1")
    if os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22") != "24":
        raise RuntimeError("set SPIN_EBCH_LOWER_WEIGHT=24")
    if int(os.environ.get("SPIN_EBCH_DISTANCE", str(DISTANCE))) != DISTANCE:
        raise RuntimeError("set SPIN_EBCH_DISTANCE=209715")

    spectrum, good = conditioned_spectrum()
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    pending = [(MIN_OCCUPATION, L, ())]
    accepted = []
    processed = 0
    while pending:
        lower, upper, seeds = pending.pop()
        processed += 1
        receipt, child_seeds = certify_interval(
            envelope, spectrum, lower, upper, seeds
        )
        if receipt is not None:
            accepted.append(receipt)
            continue
        if lower == upper:
            raise ArithmeticError(f"failed to certify occupation {lower}")
        midpoint = (lower + upper) // 2
        pending.append((midpoint + 1, upper, child_seeds))
        pending.append((lower, midpoint, child_seeds))
        if processed % 10 == 0:
            print(
                f"processed={processed},accepted={len(accepted)},pending={len(pending)}",
                flush=True,
            )

    accepted.sort(key=lambda row: int(row["occupation_interval"][0]))
    cursor = MIN_OCCUPATION
    for row in accepted:
        lower, upper = map(int, row["occupation_interval"])
        if lower != cursor:
            raise ArithmeticError("interval cover has a gap or overlap")
        cursor = upper + 1
    if cursor != L + 1:
        raise ArithmeticError("interval cover does not reach the final occupation")

    logs = [float(row["log2_contribution_upper"]) for row in accepted]
    largest = max(logs)
    aggregate_log2 = largest + math.log2(sum(2.0 ** (value - largest) for value in logs))
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-b256-one-band-q-cover-d10-v1",
        "status": "COMPLETE_BINARY64_DIAGNOSTIC_INTERVAL_COVER",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "distance": DISTANCE,
            "relative_bad_weight": DISTANCE / N,
            "occupation_interval": [MIN_OCCUPATION, L],
            "outer_weight_band": list(BAND),
            "conditioning_good_probability_diagnostic": good,
            "cell_margin_bits": CELL_MARGIN_BITS,
        },
        "convexity_rule": (
            "For a fixed witness with negative raw reference exponent, the raw "
            "exponent and final exponent are convex in the real occupation. "
            "Endpoint checks therefore cover every integer occupation in an interval."
        ),
        "processed_intervals": processed,
        "accepted_intervals": len(accepted),
        "aggregate_log2_upper": aggregate_log2,
        "aggregate_margin_bits": -aggregate_log2,
        "intervals": accepted,
        "limitations": [
            "Witness discovery uses nearest binary64 optimization.",
            "Final diagnostic transfer evaluations use non-interval high precision.",
            "An independent outward verifier must recompute every fixed witness."
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "processed_intervals": processed,
        "accepted_intervals": len(accepted),
        "aggregate_margin_bits": -aggregate_log2,
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
